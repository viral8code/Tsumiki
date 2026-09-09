using System.Collections.Concurrent;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Core.Preprocessing
{
    /// <summary>
    /// カットオフで落ちた k-mer のうち、リードの中で信頼できる k-mer に挟まれているものを救い上げる
    /// </summary>
    /// <remarks>
    /// カットオフは出現回数だけを見るため、カバレッジがたまたま薄い領域の真の k-mer もエラーと一緒に落ちる<br/>
    /// 落ちた場所ではグラフが千切れ、その領域は以降どの工程からも見えなくなる<br/>
    /// 前後が信頼できる k-mer で同じリードの中で連続しているという条件は
    /// 出現回数とは独立した証拠で、単独の低頻度 k-mer とは区別できる<br/>
    /// ただし通常の信頼 k-mer より根拠は弱いため、観測した回数をそのままカバレッジとして与え、名目値で水増ししない
    /// </remarks>
    internal static class MercyKmerRescuer
    {
        /// <summary>
        /// 救済の対象とする、信頼できない窓の連続長の上限<br/>
        /// 長く途切れている箇所は、カバレッジが薄いのではなく
        /// そもそも別の配列を読んでいる可能性が高くなる
        /// </summary>
        private const int 救済する連の上限 = 8;

        /// <summary>
        /// 救済に必要な観測回数
        /// </summary>
        /// <remarks>
        /// 1 回しか見ていない k-mer は、挟まれていてもエラーと区別できない
        /// </remarks>
        private const int 救済に必要な観測数 = 2;

        /// <summary>
        /// 落ちた k-mer を救済し、実際に足した件数を返す
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <returns>救済して足した k-mer の件数</returns>
        public static int Get_救済数(
            Parameters p_引数, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            // 候補は数百万件になりうるので、ワーカーごとに辞書を持つとその本数だけ複製することになるため、1 つを共有する
            // 値に塩基列そのものを持つのは、キーが k > 64 でハッシュになり配列を戻せなくなるためで、
            // 救済は集合へ足す処理なので実体が要る
            ConcurrentDictionary<UInt128, (int A_観測数, byte[] A_kmer)> l_候補 = [];

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            ReadPipeline.V_実行(
                l_スレッド数,
                l_スレッド数 * 256,
                FastqReader.Get_生リード列(p_引数.A_リード1のパス, p_引数.A_リード2のパス),
                (l_リード, _) => V_集める_1リード(l_リード, p_kmerインデックス, p_k長, l_候補));

            var l_追加数 = 0;
            foreach (var (_, l_候補中身) in l_候補)
            {
                if (l_候補中身.A_観測数 < 救済に必要な観測数)
                {
                    continue;
                }
                if (p_kmerインデックス.V_追加_信頼kmer(l_候補中身.A_kmer, (ulong)l_候補中身.A_観測数))
                {
                    l_追加数++;
                }
            }

            Logger.V_出力(メッセージID.救済したkmer数, l_追加数, l_候補.Count);
            return l_追加数;
        }

        /// <summary>
        /// 1 本のリードから救済候補を集める
        /// </summary>
        /// <remarks>
        /// 信頼できない窓の連なりが、両側を信頼できる窓に挟まれている場合だけを候補にする
        /// </remarks>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_候補">集めた救済候補</param>
        private static void V_集める_1リード(
            string p_リード, TrustedKmerIndex p_kmerインデックス, int p_k長,
            ConcurrentDictionary<UInt128, (int A_観測数, byte[] A_kmer)> p_候補)
        {
            if (p_リード.Length < p_k長 + 2)
            {
                // 両側に信頼できる窓を要求する以上、窓が 3 つ取れなければ意味がない
                return;
            }

            var l_塩基列 = Util.V_変換_塩基列(p_リード);
            var l_窓数 = l_塩基列.Length - p_k長 + 1;

            // 曖昧塩基を含む窓は候補にできない (パックも登録もできない) ので、
            // 信頼できるかどうかも見ないまま、連を切る壁として扱う
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
                        var l_窓 = l_塩基列.AsSpan(j, p_k長);
                        var l_キー = KmerPacking.Get_正規化キー(l_窓);
                        var l_控え = l_窓.ToArray();
                        _ = p_候補.AddOrUpdate(
                            l_キー,
                            _ => (1, l_控え),
                            (_, l_既存) => (l_既存.A_観測数 + 1, l_既存.A_kmer));
                    }
                }
                i = l_終わり;
            }
        }
    }
}
