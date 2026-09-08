using System.Collections.Concurrent;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// カットオフで落ちた k-mer のうち、リードの中で信頼できる k-mer に
    /// 挟まれているものを救い上げる。
    ///
    /// カットオフは出現回数だけを見るため、カバレッジがたまたま薄い領域の
    /// 真の k-mer もエラーと一緒に落ちる。落ちた場所ではグラフが千切れ、
    /// その領域は以降どの工程からも見えなくなる。
    ///
    /// 「前後が信頼できる k-mer で、同じリードの中で連続している」という
    /// 条件は出現回数とは独立した証拠で、単独の低頻度 k-mer とは区別できる。
    /// ただし通常の信頼 k-mer より根拠は弱いため、観測した回数をそのまま
    /// カバレッジとして与え、名目値で水増ししない。
    /// </summary>
    internal static class MercyKmerRescuer
    {
        /// <summary>
        /// 救済の対象とする、信頼できない窓の連続長の上限。
        /// 長く途切れている箇所は、カバレッジが薄いのではなく
        /// そもそも別の配列を読んでいる可能性が高くなる。
        /// </summary>
        private const int 救済する連の上限 = 8;

        /// <summary>
        /// 救済に必要な観測回数。1回しか見ていない k-mer は、
        /// 挟まれていてもエラーと区別できない。
        /// </summary>
        private const int 救済に必要な観測数 = 2;

        /// <summary>
        /// 落ちた k-mer を救済し、実際に足した件数を返す。
        /// k が 64 を超える場合は 2bit パックが UInt128 に収まらないため
        /// 何もせず 0 を返す。
        /// </summary>
        public static int Get_救済数(
            Parameters p_引数, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            if (p_k長 > 64)
            {
                Logger.V_出力(メッセージID.救済_対象外のk長, p_k長);
                return 0;
            }

            // 候補は数百万件になりうる。ワーカーごとに辞書を持つと
            // その本数だけ複製することになるため、1つを共有する。
            ConcurrentDictionary<UInt128, int> l_候補 = [];

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            ReadPipeline.V_実行(
                l_スレッド数,
                l_スレッド数 * 256,
                FastqReader.Get_生リード列(p_引数.A_リード1のパス, p_引数.A_リード2のパス),
                (l_リード, _) => V_集める_1リード(l_リード, p_kmerインデックス, p_k長, l_候補));

            var l_追加数 = 0;
            var l_kmer = new byte[p_k長];
            foreach (var (l_正規形, l_観測数) in l_候補)
            {
                if (l_観測数 < 救済に必要な観測数)
                {
                    continue;
                }
                V_復元_塩基列(l_正規形, p_k長, l_kmer);
                if (p_kmerインデックス.V_追加_信頼kmer(l_kmer, (ulong)l_観測数))
                {
                    l_追加数++;
                }
            }

            Logger.V_出力(メッセージID.救済したkmer数, l_追加数, l_候補.Count);
            return l_追加数;
        }

        /// <summary>
        /// 1本のリードから救済候補を集める。信頼できない窓の連なりが
        /// 両側を信頼できる窓に挟まれている場合だけを候補にする。
        /// </summary>
        private static void V_集める_1リード(
            string p_リード, TrustedKmerIndex p_kmerインデックス, int p_k長,
            ConcurrentDictionary<UInt128, int> p_候補)
        {
            if (p_リード.Length < p_k長 + 2)
            {
                // 両側に信頼できる窓を要求する以上、窓が3つ取れなければ意味がない。
                return;
            }

            var l_塩基列 = Util.V_変換_塩基列(p_リード);
            var l_窓数 = l_塩基列.Length - p_k長 + 1;

            // 曖昧塩基を含む窓は候補にできない(パックも登録もできない)。
            // 信頼できるかどうかも見ないまま、連を切る壁として扱う。
            var l_有効 = new bool[l_窓数];
            var l_信頼 = new bool[l_窓数];
            var l_曖昧数 = 0;
            for (var i = 0; i < p_k長; i++)
            {
                if (l_塩基列[i] is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    l_曖昧数++;
                }
            }
            for (var i = 0; i < l_窓数; i++)
            {
                if (i > 0)
                {
                    if (l_塩基列[i - 1] is < Consts.塩基ID.A or > Consts.塩基ID.T)
                    {
                        l_曖昧数--;
                    }
                    if (l_塩基列[i + p_k長 - 1] is < Consts.塩基ID.A or > Consts.塩基ID.T)
                    {
                        l_曖昧数++;
                    }
                }
                l_有効[i] = l_曖昧数 == 0;
                l_信頼[i] = l_有効[i] && p_kmerインデックス.Get_含まれるか(l_塩基列.AsSpan(i, p_k長));
            }

            for (var i = 1; i < l_窓数 - 1; i++)
            {
                if (l_信頼[i] || !l_有効[i] || !l_信頼[i - 1])
                {
                    continue;
                }

                var l_終わり = i;
                while (l_終わり < l_窓数 && !l_信頼[l_終わり] && l_有効[l_終わり])
                {
                    l_終わり++;
                }
                var l_連の長さ = l_終わり - i;
                if (l_終わり < l_窓数 && l_信頼[l_終わり] && l_連の長さ <= 救済する連の上限)
                {
                    for (var j = i; j < l_終わり; j++)
                    {
                        var l_正規形 = KmerPacking.Get_正規化パック(l_塩基列.AsSpan(j, p_k長));
                        _ = p_候補.AddOrUpdate(l_正規形, 1, (_, l_既存) => l_既存 + 1);
                    }
                }
                i = l_終わり;
            }
        }

        /// <summary>パック済みの正規形を塩基ID列へ戻す。</summary>
        private static void V_復元_塩基列(UInt128 p_パック済み, int p_k長, byte[] p_出力)
        {
            for (var i = 0; i < p_k長; i++)
            {
                var l_コドン = (int)((p_パック済み >> (2 * (p_k長 - 1 - i))) & 3);
                p_出力[i] = (byte)(l_コドン + 1);
            }
        }
    }
}
