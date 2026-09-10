using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// アセンブリが観測された k-mer とその出現回数に対して辻褄が合っているかの自己検査
    /// </summary>
    /// <remarks>
    /// 取りこぼし (信頼集合にあるのにアセンブリに現れない) は削りすぎか経路から漏れた領域を、
    /// 出しすぎ (コピー数を超えて現れる) は反復の複製過多か同じ領域の二重組み立てを意味する<br/>
    /// どちらもゼロにはならない (エラー由来の k-mer が信頼集合に残れば取りこぼし側に出る) ため、
    /// 絶対値ではなく変更前後の増減を見る指標になる
    /// </remarks>
    internal static class AssemblyValidator
    {
        #region 公開メソッド

        /// <summary>
        /// FASTA のアセンブリを、k-mer インデックスが保持する信頼できる k-mer 集合と
        /// 突き合わせる
        /// </summary>
        /// <param name="p_FASTAパス">突き合わせる FASTA のパス</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_単一コピー基準値">その k-mer が何回現れてよいかをカバレッジから見積もるための基準値</param>
        /// <returns>自己検査の結果</returns>
        public static 整合性検査結果? Get_検査結果(string p_FASTAパス, TrustedKmerIndex p_kmerインデックス, int p_k長, double p_単一コピー基準値)
        {
            // 逆相補は同一視して数える
            // キーは 2 bit パックした UInt128 で、k が 64 を超えるとパックが収まらないので正規形のハッシュに切り替える
            // 数えるだけで配列を戻さないので、ハッシュで足りる
            Dictionary<UInt128, int> l_観測 = [];
            var l_延べ数 = 0L;

            using (var l_読み込み = new FastaReader(p_FASTAパス))
            {
                while (l_読み込み.Get_続きがあるか())
                {
                    var l_配列 = l_読み込み.Get_次の配列().A_配列;
                    for (var i = 0; i + p_k長 <= l_配列.Length; i++)
                    {
                        if (!KmerPacking.Get_正規化キー(l_配列, i, p_k長, out var l_正規形))
                        {
                            continue;
                        }
                        l_観測[l_正規形] = l_観測.GetValueOrDefault(l_正規形) + 1;
                        l_延べ数++;
                    }
                }
            }

            var l_信頼kmer数 = 0L;
            var l_取りこぼし数 = 0L;
            var l_出しすぎ種類数 = 0L;
            var l_余分な延べ数 = 0L;

            foreach (var l_kmer in p_kmerインデックス.Get_信頼kmer一覧())
            {
                l_信頼kmer数++;
                var l_正規形 = KmerPacking.Get_正規化キー(l_kmer);

                var l_出現数 = l_観測.GetValueOrDefault(l_正規形);
                if (l_出現数 == 0)
                {
                    l_取りこぼし数++;
                    continue;
                }

                // カバレッジから期待されるコピー数
                // 基準値が取れていない場合は判定を諦める (1 コピー扱いにすると全部を過剰と誤判定してしまう)
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

        /// <summary>
        /// 自己検査の結果をログへ出力する
        /// </summary>
        /// <param name="p_ラベル">どの成果物に対する検査かを表す名前</param>
        /// <param name="p_検査結果">自己検査の結果、調べられなかった場合は null</param>
        public static void V_出力_検査結果(string p_ラベル, 整合性検査結果? p_検査結果)
        {
            if (p_検査結果 is not { } p_結果)
            {
                Logger.V_出力(メッセージID.検査_対象外のk長, p_ラベル);
                return;
            }

            Logger.V_出力(メッセージID.検査_取りこぼし, p_ラベル, p_結果.A_信頼kmer数, p_結果.A_取りこぼし数, p_結果.A_取りこぼし率);
            Logger.V_出力(メッセージID.検査_出しすぎ, p_ラベル, p_結果.A_アセンブリ内の延べ数, p_結果.A_余分な延べ数, p_結果.A_出しすぎ率, p_結果.A_出しすぎkmer種類数);
        }

        #endregion
    }
}
