using Tsumiki.IO;

namespace Tsumiki.Core
{
    /// <summary>
    /// アセンブリ結果(unitigs/contigs/scaffolds)の品質を大まかに把握するための
    /// 基本統計量(N50/L50/総延長/配列数/GC含量)。リファレンスなしで計算できる
    /// 範囲の指標のみを対象とする。
    /// </summary>
    internal readonly struct AssemblyStats(
        int sequenceCount,
        long totalLength,
        int maxLength,
        int minLength,
        int n50,
        int l50,
        double gcPercent)
    {
        public readonly int SequenceCount = sequenceCount;
        public readonly long TotalLength = totalLength;
        public readonly int MaxLength = maxLength;
        public readonly int MinLength = minLength;

        /// <summary>
        /// 長い順に並べて累積長が全長の50%に達した時点の配列長。
        /// </summary>
        public readonly int N50 = n50;

        /// <summary>
        /// N50 に達するまでに必要だった配列の本数(1始まり)。
        /// </summary>
        public readonly int L50 = l50;

        public readonly double GcPercent = gcPercent;

        public override string ToString()
        {
            return $"count={this.SequenceCount}, total_length={this.TotalLength}, " +
                   $"N50={this.N50}, L50={this.L50}, max={this.MaxLength}, min={this.MinLength}, " +
                   $"GC%={this.GcPercent:0.00}";
        }
    }

    internal static class AssemblyStatsReporter
    {
        public static AssemblyStats Compute(IEnumerable<string> sequences)
        {
            var lengths = new List<int>();
            long totalLength = 0;
            long gcCount = 0;
            long baseCount = 0;

            foreach (var seq in sequences)
            {
                lengths.Add(seq.Length);
                totalLength += seq.Length;
                foreach (var c in seq)
                {
                    if (c is 'N' or 'n')
                    {
                        continue;
                    }
                    baseCount++;
                    if (c is 'G' or 'g' or 'C' or 'c')
                    {
                        gcCount++;
                    }
                }
            }

            if (lengths.Count == 0)
            {
                return new AssemblyStats(0, 0, 0, 0, 0, 0, 0);
            }

            lengths.Sort();
            lengths.Reverse();

            var half = totalLength / 2.0;
            long cumulative = 0;
            var n50 = lengths[^1];
            var l50 = lengths.Count;
            for (var i = 0; i < lengths.Count; i++)
            {
                cumulative += lengths[i];
                if (cumulative >= half)
                {
                    n50 = lengths[i];
                    l50 = i + 1;
                    break;
                }
            }

            var gcPercent = baseCount == 0 ? 0.0 : (100.0 * gcCount / baseCount);

            return new AssemblyStats(
                sequenceCount: lengths.Count,
                totalLength: totalLength,
                maxLength: lengths[0],
                minLength: lengths[^1],
                n50: n50,
                l50: l50,
                gcPercent: gcPercent);
        }

        public static AssemblyStats ComputeFromFasta(string fastaPath)
        {
            return Compute(ReadSequences(fastaPath));
        }

        private static IEnumerable<string> ReadSequences(string fastaPath)
        {
            using var reader = new FastaReader(fastaPath);
            while (reader.HasNext())
            {
                yield return reader.NextSequence().Seq;
            }
        }

        /// <summary>
        /// 他アセンブラとの比較で慣習的に使われる最小長。ABySS 付属の
        /// abyss-fac が既定でこの長さ以上の配列だけを集計するため、
        /// 公表されている N50 等はほぼこの条件で計算されている。
        /// 全件の統計だけを出していると、数十bpの断片が大量に混じった
        /// こちらの数字と比較して不当に悪く見える(あるいはその逆になる)。
        /// </summary>
        public const int ComparableMinLength = 500;

        /// <summary>
        /// fastaPath の統計量を計算し、"[Stats] label: ..." の形式でコンソールへ出力する。
        /// 全配列を対象とした統計に加えて、他アセンブラの公表値と直接比較できるよう
        /// ComparableMinLength 以上の配列だけに絞った統計も併記する。
        /// </summary>
        public static void Report(string label, string fastaPath)
        {
            if (!File.Exists(fastaPath))
            {
                Console.WriteLine($"[Stats] {label}: (file not found: {fastaPath})");
                return;
            }
            var stats = ComputeFromFasta(fastaPath);
            Console.WriteLine($"[Stats] {label}: {stats}");

            var filtered = Compute(ReadSequences(fastaPath).Where(s => s.Length >= ComparableMinLength));
            Console.WriteLine($"[Stats] {label} (>={ComparableMinLength}bp, comparable to abyss-fac): {filtered}");
        }
    }
}
