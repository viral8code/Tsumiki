using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// アセンブリが観測された k-mer とその出現回数に対して辻褄が合っているかの自己検査。
    ///
    /// 取りこぼし(信頼集合にあるのにアセンブリに現れない)は削りすぎか
    /// 経路から漏れた領域を、出しすぎ(コピー数を超えて現れる)は
    /// 反復の複製過多か同じ領域の二重組み立てを意味する。
    ///
    /// どちらもゼロにはならない(エラー由来の k-mer が信頼集合に残れば
    /// 取りこぼし側に出る)。絶対値ではなく変更前後の増減を見る指標。
    /// </summary>
    internal static class AssemblyValidator
    {
        /// <summary>
        /// FASTA のアセンブリを、k-mer インデックスが保持する信頼できる k-mer 集合と
        /// 突き合わせる。p_単一コピー基準値 は「その k-mer が何回現れてよいか」を
        /// カバレッジから見積もるために使う。
        /// </summary>
        public static 整合性検査結果 Get_検査結果(
            string p_FASTAパス, TrustedKmerIndex p_kmerインデックス, int p_k長, double p_単一コピー基準値)
        {
            // 逆相補は同一視して数える。キーは 2bit パックした UInt128 で、
            // 文字列キーだとアセンブリ規模で 1GB を超え、同時に生きている
            // k-mer インデックスと合わせてメモリが厳しくなる。
            if (p_k長 > 64)
            {
                Console.WriteLine($"[Check] {Path.GetFileName(p_FASTAパス)}: skipped (self-check currently supports k <= 64 only).");
                return default;
            }

            Dictionary<UInt128, int> l_観測 = [];
            long l_延べ数 = 0;

            using (var l_読み込み = new FastaReader(p_FASTAパス))
            {
                while (l_読み込み.Get_続きがあるか())
                {
                    var l_配列 = l_読み込み.Get_次の配列().A_配列;
                    for (var i = 0; i + p_k長 <= l_配列.Length; i++)
                    {
                        if (!KmerPacking.Get_正規化パック(l_配列, i, p_k長, out var l_正規形))
                        {
                            continue;
                        }
                        l_観測[l_正規形] = l_観測.GetValueOrDefault(l_正規形) + 1;
                        l_延べ数++;
                    }
                }
            }

            long l_信頼kmer数 = 0;
            long l_取りこぼし数 = 0;
            long l_出しすぎ種類数 = 0;
            long l_余分な延べ数 = 0;

            foreach (var l_kmer in p_kmerインデックス.Get_信頼kmer一覧())
            {
                l_信頼kmer数++;
                var l_正規形 = KmerPacking.Get_正規化パック(l_kmer);

                var l_出現数 = l_観測.GetValueOrDefault(l_正規形);
                if (l_出現数 == 0)
                {
                    l_取りこぼし数++;
                    continue;
                }

                // カバレッジから期待されるコピー数。基準値が取れていない場合は
                // 判定を諦める(1コピー扱いにすると全部を過剰と誤判定してしまう)。
                if (p_単一コピー基準値 <= 0)
                {
                    continue;
                }
                var l_カバレッジ = p_kmerインデックス.Get_カバレッジ(l_kmer);
                var l_期待コピー数 = Math.Max(1, (int)Math.Round(l_カバレッジ / p_単一コピー基準値));

                if (l_出現数 > l_期待コピー数)
                {
                    l_出しすぎ種類数++;
                    l_余分な延べ数 += l_出現数 - l_期待コピー数;
                }
            }

            return new 整合性検査結果(l_信頼kmer数, l_延べ数, l_観測.Count, l_取りこぼし数, l_出しすぎ種類数, l_余分な延べ数);
        }

        public static void V_出力_検査結果(string p_ラベル, 整合性検査結果 p_結果)
        {
            Console.WriteLine(
                $"[Check] {p_ラベル}: {p_結果.A_信頼kmer数:N0} trusted k-mer(s); " +
                $"{p_結果.A_取りこぼし数:N0} ({p_結果.A_取りこぼし率:0.00}%) do not appear in the assembly at all " +
                "(sequence that was trimmed away or never reached by any path).");
            Console.WriteLine(
                $"[Check] {p_ラベル}: {p_結果.A_アセンブリ内の延べ数:N0} k-mer instance(s) in the assembly; " +
                $"{p_結果.A_余分な延べ数:N0} ({p_結果.A_出しすぎ率:0.00}%) are more copies than the coverage supports " +
                $"across {p_結果.A_出しすぎkmer種類数:N0} distinct k-mer(s) -- this is the part of the total length that is inflated.");
        }
    }
}
