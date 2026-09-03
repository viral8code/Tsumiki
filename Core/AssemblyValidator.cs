using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// 出来上がったアセンブリが、観測された k-mer とその出現回数に対して
    /// 辻褄が合っているかを検証する。リファレンス配列を使わずに
    /// 「自分の出した答えが入力データと整合しているか」を確かめる自己検査。
    ///
    /// 見るのは次の2つのずれ:
    ///
    /// 1. 取りこぼし(missing): 信頼できる k-mer 集合には存在するのに、
    ///    アセンブリのどこにも現れない k-mer。グラフ簡略化で削りすぎたか、
    ///    分岐が解けずに経路から漏れた領域を意味する。
    ///
    /// 2. 出しすぎ(over-represented): カバレッジから推定されるコピー数よりも
    ///    多くアセンブリ中に現れる k-mer。反復配列の複製(repeat resolution)が
    ///    行き過ぎた場合や、同じ領域を2通りに組み立ててしまった場合に出る。
    ///    総延長が実際のゲノムサイズより大きくなる原因は、ほぼこれ。
    ///
    /// どちらも「多少はある」のが普通で、ゼロにはならない(エラー由来の
    /// k-mer が信頼集合に残っていれば missing 側に出る)。絶対値ではなく、
    /// 変更の前後でこの割合がどう動いたかを見るための指標として使う。
    /// </summary>
    internal static class AssemblyValidator
    {
        internal readonly record struct Result(
            long TrustedKmers,
            long AssemblyKmerInstances,
            long DistinctKmersInAssembly,
            long MissingKmers,
            long OverRepresentedKmers,
            long ExcessInstances)
        {
            /// <summary>信頼できる k-mer のうち、アセンブリに現れなかった割合。</summary>
            public double MissingPercent => this.TrustedKmers == 0 ? 0 : 100.0 * this.MissingKmers / this.TrustedKmers;

            /// <summary>
            /// アセンブリ中の k-mer 延べ数のうち、コピー数の推定を超えて
            /// 余分に現れている分の割合。総延長の水増し量にほぼ対応する。
            /// </summary>
            public double ExcessPercent => this.AssemblyKmerInstances == 0 ? 0 : 100.0 * this.ExcessInstances / this.AssemblyKmerInstances;
        }

        /// <summary>
        /// fastaPath のアセンブリを、index が保持する信頼できる k-mer 集合と
        /// 突き合わせる。expectedCopies は「その k-mer が何回現れてよいか」を
        /// カバレッジから見積もるための単一コピー基準値。
        /// </summary>
        public static Result Validate(string fastaPath, TrustedKmerIndex index, int kmerLength, double singleCopyBaseline)
        {
            // アセンブリ中に各 k-mer が何回現れるかを数える。
            // 正規化(canonical)して数えるため、逆相補は同じものとして扱う。
            Dictionary<string, int> observed = [];
            long instances = 0;

            using (var reader = new FastaReader(fastaPath))
            {
                while (reader.HasNext())
                {
                    var seq = reader.NextSequence().Seq;
                    for (var i = 0; i + kmerLength <= seq.Length; i++)
                    {
                        var kmer = seq.Substring(i, kmerLength);
                        if (kmer.Contains('N'))
                        {
                            continue;
                        }
                        var canonical = Canonical(kmer);
                        observed[canonical] = observed.GetValueOrDefault(canonical) + 1;
                        instances++;
                    }
                }
            }

            long trusted = 0;
            long missing = 0;
            long overRepresented = 0;
            long excessInstances = 0;

            foreach (var kmerBytes in index.EnumerateTrustedKmers())
            {
                trusted++;
                var text = string.Concat(kmerBytes.Select(Util.ByteToBaseString));
                var canonical = Canonical(text);

                var seen = observed.GetValueOrDefault(canonical);
                if (seen == 0)
                {
                    missing++;
                    continue;
                }

                // カバレッジから期待されるコピー数。基準値が取れていない場合は
                // 判定を諦める(1コピー扱いにすると全部を過剰と誤判定してしまう)。
                if (singleCopyBaseline <= 0)
                {
                    continue;
                }
                var coverage = index.GetCoverage(kmerBytes);
                var expected = Math.Max(1, (int)Math.Round(coverage / singleCopyBaseline));

                if (seen > expected)
                {
                    overRepresented++;
                    excessInstances += seen - expected;
                }
            }

            return new Result(trusted, instances, observed.Count, missing, overRepresented, excessInstances);
        }

        private static string Canonical(string kmer)
        {
            var revComp = Util.ReverseComprement(kmer);
            return string.CompareOrdinal(kmer, revComp) <= 0 ? kmer : revComp;
        }

        public static void Report(string label, Result result)
        {
            Console.WriteLine(
                $"[Check] {label}: {result.TrustedKmers:N0} trusted k-mer(s); " +
                $"{result.MissingKmers:N0} ({result.MissingPercent:0.00}%) do not appear in the assembly at all " +
                "(sequence that was trimmed away or never reached by any path).");
            Console.WriteLine(
                $"[Check] {label}: {result.AssemblyKmerInstances:N0} k-mer instance(s) in the assembly; " +
                $"{result.ExcessInstances:N0} ({result.ExcessPercent:0.00}%) are more copies than the coverage supports " +
                $"across {result.OverRepresentedKmers:N0} distinct k-mer(s) -- this is the part of the total length that is inflated.");
        }
    }
}
