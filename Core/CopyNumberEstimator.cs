using Tsumiki.Common;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// 各 unitig がゲノム中に何回現れるか(コピー数)を、カバレッジから推定する。
    ///
    /// 考え方: ゲノム中に1回しか現れない領域のカバレッジを基準値(baseline)とすると、
    /// n 回現れる反復配列にはリードが n 倍集まるので、カバレッジも約 n 倍になる。
    /// したがって「その unitig の平均カバレッジ / baseline」を丸めればコピー数になる。
    ///
    /// なぜ必要か:
    /// - 反復配列かどうかを、グラフの形(入次数・出次数)ではなく量的な根拠で判定できる。
    ///   入次数2・出次数2でも実は単一コピー(バブルの残骸)ということがありうるし、
    ///   逆に次数が1でも高カバレッジならタンデムリピートを1本に潰している疑いがある。
    /// - ゲノム全体の経路探索(GenomePathSearch)で「この unitig は2回使ってよい」という
    ///   予算になる。予算が無いと、探索は同じ反復を何度でも通れてしまい破綻する。
    ///
    /// baseline は全 unitig の平均カバレッジの「長さ加重中央値」を使う。
    /// 単純平均・単純中央値だと、本数として多い短い断片(エラー由来の残骸など)に
    /// 引きずられる。塩基数で重み付けすれば、ゲノムの大部分を占める単一コピー領域の
    /// 水準が選ばれる。
    /// </summary>
    internal static class CopyNumberEstimator
    {
        /// <summary>
        /// これを下回るカバレッジ比の unitig は、コピー数を推定できるだけの
        /// 根拠が無いとみなして 1 として扱う(0 コピーにはしない)。
        /// </summary>
        private const double MinimumRatioForMultiCopy = 1.5;

        /// <summary>
        /// コピー数の上限。これを超える比が出た場合、rRNA オペロンのような
        /// 高コピー反復か、あるいはカバレッジ異常のどちらかで区別がつかない。
        /// 経路探索の予算としては大きすぎると探索が発散するため頭打ちにする。
        /// </summary>
        private const int MaximumCopyNumber = 12;

        internal readonly record struct Result(
            double Baseline,
            IReadOnlyDictionary<int, double> Coverage,
            IReadOnlyDictionary<int, int> CopyNumber);

        /// <summary>
        /// unitig ID(1始まり) -> その unitig を構成する k-mer の平均カバレッジ、を計算する。
        /// </summary>
        public static Dictionary<int, double> ComputeCoverage(
            TrustedKmerIndex index,
            IReadOnlyDictionary<int, string> unitigSequences,
            int kmerLength)
        {
            Dictionary<int, double> coverage = [];
            foreach (var (id, sequence) in unitigSequences)
            {
                if (sequence.Length < kmerLength)
                {
                    coverage[id] = 0;
                    continue;
                }

                var bytes = new byte[sequence.Length];
                for (var i = 0; i < sequence.Length; i++)
                {
                    bytes[i] = Util.GetSimpleNucleotideID(sequence[i]);
                }

                ulong sum = 0;
                var count = 0;
                for (var i = 0; i + kmerLength <= bytes.Length; i++)
                {
                    sum += index.GetCoverage(bytes.AsSpan(i, kmerLength));
                    count++;
                }
                coverage[id] = count == 0 ? 0 : (double)sum / count;
            }
            return coverage;
        }

        /// <summary>
        /// カバレッジからコピー数を推定する。
        /// </summary>
        public static Result Estimate(
            IReadOnlyDictionary<int, double> coverage,
            IReadOnlyDictionary<int, int> unitigLengths)
        {
            var baseline = WeightedMedian(coverage, unitigLengths);

            Dictionary<int, int> copyNumber = [];
            foreach (var (id, cov) in coverage)
            {
                if (baseline <= 0)
                {
                    copyNumber[id] = 1;
                    continue;
                }

                var ratio = cov / baseline;
                if (ratio < MinimumRatioForMultiCopy)
                {
                    // 単一コピー(あるいは低カバレッジで判断できない)。
                    // 0 にはしない: 実際に配列は存在しており、経路から
                    // 締め出してしまうと組み立てられなくなる。
                    copyNumber[id] = 1;
                    continue;
                }

                copyNumber[id] = Math.Clamp((int)Math.Round(ratio), 1, MaximumCopyNumber);
            }

            return new Result(baseline, coverage, copyNumber);
        }

        /// <summary>
        /// 長さで重み付けしたカバレッジの中央値。ゲノムの大部分を占める
        /// 単一コピー領域の水準を推定するために使う。
        /// </summary>
        private static double WeightedMedian(
            IReadOnlyDictionary<int, double> coverage,
            IReadOnlyDictionary<int, int> unitigLengths)
        {
            var pairs = coverage
                .Where(kv => unitigLengths.ContainsKey(kv.Key) && kv.Value > 0)
                .Select(kv => (Length: (long)unitigLengths[kv.Key], Coverage: kv.Value))
                .OrderBy(p => p.Coverage)
                .ToList();
            if (pairs.Count == 0)
            {
                return 0;
            }

            var totalLength = pairs.Sum(p => p.Length);
            if (totalLength == 0)
            {
                return 0;
            }

            var half = totalLength / 2.0;
            long cumulative = 0;
            foreach (var (length, cov) in pairs)
            {
                cumulative += length;
                if (cumulative >= half)
                {
                    return cov;
                }
            }
            return pairs[^1].Coverage;
        }

        /// <summary>
        /// 推定結果の要約をコンソールへ出力する。
        /// 「単一コピーが何本・何bp、2コピー以上が何本・何bp」が分かると、
        /// 反復配列がアセンブリのどれだけを占めているかが把握できる。
        /// </summary>
        public static void Report(Result result, IReadOnlyDictionary<int, int> unitigLengths)
        {
            Console.WriteLine($"[Info] Single-copy coverage baseline estimated as {result.Baseline:0.#} (length-weighted median).");

            var byCopy = result.CopyNumber
                .GroupBy(kv => kv.Value)
                .OrderBy(g => g.Key)
                .Select(g => (Copy: g.Key,
                              Count: g.Count(),
                              Bases: g.Sum(kv => (long)unitigLengths.GetValueOrDefault(kv.Key, 0))))
                .ToList();

            var summary = string.Join(", ", byCopy.Select(b => $"x{b.Copy}: {b.Count} unitig(s)/{b.Bases:N0}bp"));
            Console.WriteLine($"[Info] Estimated copy numbers -- {summary}");

            var repeatBases = byCopy.Where(b => b.Copy >= 2).Sum(b => b.Bases);
            var totalBases = byCopy.Sum(b => b.Bases);
            if (totalBases > 0)
            {
                Console.WriteLine(
                    $"[Info] Multi-copy (repeat) content: {repeatBases:N0}bp of {totalBases:N0}bp " +
                    $"({100.0 * repeatBases / totalBases:0.0}% of the assembly is sequence that occurs more than once).");
            }
        }
    }
}
