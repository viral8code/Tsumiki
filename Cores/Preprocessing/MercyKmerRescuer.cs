using System.Collections.Concurrent;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Preprocessing
{
    /// <summary>
    /// カットオフで落ちた k-mer のうち、リードの中で信頼できる k-mer に挟まれているものを救い上げる
    /// </summary>
    internal static class MercyKmerRescuer
    {
        #region 定数

        /// <summary>
        /// 救済の対象とする、信頼できない窓の連続長の上限
        /// </summary>
        private const int 救済する連の上限 = 8;

        /// <summary>
        /// 救済に必要な観測回数
        /// </summary>
        private const int 救済に必要な観測数 = 2;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 落ちた k-mer を救済し、実際に足した件数を返す
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <returns>救済して足した k-mer の件数</returns>
        public static int Get_救済数(Parameters p_引数, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            using var l_計測 = new StageTimer($"mercy k={p_k長}");
            ConcurrentDictionary<UInt128, (int A_観測数, byte[] A_kmer)> l_候補 = [];

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, FastqReader.Get_生リード列([.. p_引数.A_ライブラリ群.SelectMany(x => new[] { x.A_リード1, x.A_リード2 })]), (l_リード, _) => V_集める_1リード(l_リード, p_kmerインデックス, p_k長, l_候補));

            var l_追加数 = 0;
            foreach (var (_, l_候補中身) in l_候補)
            {
                if (l_候補中身.A_観測数 < 救済に必要な観測数)
                {
                    continue;
                }
                if (p_kmerインデックス.Try追加_信頼kmer(l_候補中身.A_kmer, (ulong)l_候補中身.A_観測数))
                {
                    l_追加数++;
                }
            }

            Logger.V_出力(メッセージID.救済したkmer数, l_追加数, l_候補.Count);
            return l_追加数;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 本のリードから救済候補を集める
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_候補">集めた救済候補</param>
        private static void V_集める_1リード(string p_リード, TrustedKmerIndex p_kmerインデックス, int p_k長, ConcurrentDictionary<UInt128, (int A_観測数, byte[] A_kmer)> p_候補)
        {
            if (p_リード.Length < p_k長 + 2)
            {
                return;
            }

            var l_塩基列 = Util.V_変換_塩基列(p_リード);
            var l_窓数 = l_塩基列.Length - p_k長 + 1;

            var l_有効 = new bool[l_窓数];
            var l_信頼 = new bool[l_窓数];
            if (p_k長 <= TrustedKmerIndex.パック値のk上限)
            {
                V_判定_窓_パック(p_リード, p_kmerインデックス, p_k長, l_有効, l_信頼);
            }
            else
            {
                V_判定_窓(l_塩基列, p_kmerインデックス, p_k長, l_有効, l_信頼);
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
                if (l_終わり < l_窓数 && l_信頼[l_終わり] && l_連の長さ <= 救済する連の上限 && Is行き止まり同士(p_kmerインデックス, l_塩基列.AsSpan(i - 1, p_k長), l_塩基列.AsSpan(l_終わり, p_k長)))
                {
                    for (var j = i; j < l_終わり; j++)
                    {
                        var l_窓 = l_塩基列.AsSpan(j, p_k長);
                        var l_キー = KmerPacking.TryGet_正規化キー(l_窓);
                        var l_控え = l_窓.ToArray();
                        _ = p_候補.AddOrUpdate(l_キー, _ => (1, l_控え), (_, l_既存) => (l_既存.A_観測数 + 1, l_既存.A_kmer));
                    }
                }
                i = l_終わり;
            }
        }

        /// <summary>
        /// 救済する連の左の信頼 k-mer に他の出口が無く、右の信頼 k-mer に他の入口が無いか
        /// </summary>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_左">連の直前の信頼できる k-mer</param>
        /// <param name="p_右">連の直後の信頼できる k-mer</param>
        /// <returns>両側とも行き止まりなら true</returns>
        internal static bool Is行き止まり同士(TrustedKmerIndex p_kmerインデックス, Span<byte> p_左, Span<byte> p_右)
        {
            return p_kmerインデックス.Get_出次数(p_左) == 0 && p_kmerインデックス.Get_入次数(p_右) == 0;
        }

        /// <summary>
        /// 各窓が曖昧塩基を含まないか、信頼できる k-mer かを、転がしたパック値で判定する (k &lt;= 128)
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_有効">曖昧塩基を含まない窓</param>
        /// <param name="p_信頼">信頼できる k-mer の窓</param>
        private static void V_判定_窓_パック(string p_リード, TrustedKmerIndex p_kmerインデックス, int p_k長, bool[] p_有効, bool[] p_信頼)
        {
            var l_窓 = new RollingKmer(p_k長);
            for (var i = 0; i < p_リード.Length; i++)
            {
                if (!l_窓.Try追加(p_リード[i], out var l_キー))
                {
                    continue;
                }

                var l_開始 = i - p_k長 + 1;
                p_有効[l_開始] = true;
                p_信頼[l_開始] = p_kmerインデックス.Haskmer_正規形(l_キー.A_上位, l_キー.A_下位);
            }
        }

        /// <summary>
        /// 各窓が曖昧塩基を含まないか、信頼できる k-mer かを判定する
        /// </summary>
        /// <param name="p_塩基列">リードの塩基 ID 列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_有効">曖昧塩基を含まない窓</param>
        /// <param name="p_信頼">信頼できる k-mer の窓</param>
        private static void V_判定_窓(byte[] p_塩基列, TrustedKmerIndex p_kmerインデックス, int p_k長, bool[] p_有効, bool[] p_信頼)
        {
            var l_曖昧数 = 0;
            for (var i = 0; i < p_k長; i++)
            {
                if (p_塩基列[i] is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    l_曖昧数++;
                }
            }
            for (var i = 0; i < p_有効.Length; i++)
            {
                if (i > 0)
                {
                    if (p_塩基列[i - 1] is < Consts.塩基ID.A or > Consts.塩基ID.T)
                    {
                        l_曖昧数--;
                    }

                    if (p_塩基列[i + p_k長 - 1] is < Consts.塩基ID.A or > Consts.塩基ID.T)
                    {
                        l_曖昧数++;
                    }
                }
                p_有効[i] = l_曖昧数 == 0;
                p_信頼[i] = p_有効[i] && p_kmerインデックス.Haskmer(p_塩基列.AsSpan(i, p_k長));
            }
        }

        #endregion
    }
}
