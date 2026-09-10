using System.Text;
using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Reporting;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// スキャフォールドの N を、グラフ上で両側を繋ぐ経路を探して実配列に置き換える
    /// </summary>
    /// <remarks>
    /// contig が途切れる原因は配列の不在ではなく分岐の未解決であることが多く、
    /// その場合ギャップを埋める配列は k-mer 集合の中に実在する<br/>
    /// 経路がちょうど 1 本に定まったときだけ埋める<br/>
    /// 複数見つかった場合は
    /// どれが正しいか決められないため N のまま残す<br/>
    /// 誤った配列で埋めるより、
    /// 分からないことが分かる状態のほうが下流の解析にとって安全
    /// </remarks>
    internal static class GapFiller
    {
        #region 公開メソッド

        /// <summary>
        /// スキャフォールドを読み込み、埋められるギャップを埋めて同じパスへ書き戻す
        /// </summary>
        /// <param name="p_スキャフォールドパス"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        public static ギャップ充填統計 V_充填_ギャップ(string p_スキャフォールドパス, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_スキャフォールド群 = FastaReader.Get_全エントリ(p_スキャフォールドパス);

            var l_総ギャップ数 = 0;
            var l_埋めたギャップ数 = 0;
            var l_埋めた塩基数 = 0;
            var l_一意でない数 = 0;
            var l_到達不能数 = 0;

            List<(string A_ID, string A_配列)> l_結果 = [];
            foreach (var (l_ID, l_配列) in l_スキャフォールド群)
            {
                var l_出力 = new StringBuilder();
                var l_位置 = 0;
                while (l_位置 < l_配列.Length)
                {
                    if (l_配列[l_位置] != 'N')
                    {
                        _ = l_出力.Append(l_配列[l_位置]);
                        l_位置++;
                        continue;
                    }

                    // N の連続区間 = 1 つのギャップ
                    var l_ギャップ開始 = l_位置;
                    while (l_位置 < l_配列.Length && l_配列[l_位置] == 'N')
                    {
                        l_位置++;
                    }
                    var l_ギャップ長 = l_位置 - l_ギャップ開始;
                    l_総ギャップ数++;

                    var l_埋めた配列 = Get_ギャップを埋める配列(l_出力, l_配列, l_ギャップ長, l_位置, p_kmerインデックス, p_k長, out var l_判定);
                    if (l_埋めた配列 != null)
                    {
                        _ = l_出力.Append(l_埋めた配列);
                        l_埋めたギャップ数++;
                        l_埋めた塩基数 += l_埋めた配列.Length;
                    }
                    else
                    {
                        if (l_判定 == ギャップ充填判定.一意でない)
                        {
                            l_一意でない数++;
                        }
                        else
                        {
                            l_到達不能数++;
                        }
                        AmbiguityRecorder.V_記録(
                            l_判定 == ギャップ充填判定.一意でない
                                ? 曖昧箇所の種別.経路が一意でない
                                : 曖昧箇所の種別.到達不能,
                            $"{l_ID.TrimStart('>')}:{l_ギャップ開始}-{l_位置}");
                        _ = l_出力.Append('N', l_ギャップ長);
                    }
                }
                l_結果.Add((l_ID, l_出力.ToString()));
            }

            using (var l_書き込み = new FastaWriter(p_スキャフォールドパス))
            {
                foreach (var (l_ID, l_配列) in l_結果)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }

            return new ギャップ充填統計(l_総ギャップ数, l_埋めたギャップ数, l_埋めた塩基数, l_一意でない数, l_到達不能数);
        }

        /// <summary>
        /// ギャップ充填の結果をログへ出力する
        /// </summary>
        /// <param name="p_統計">ギャップ充填の結果</param>
        public static void V_出力_充填統計(ギャップ充填統計 p_統計)
        {
            if (p_統計.A_総ギャップ数 == 0)
            {
                Logger.V_出力(メッセージID.ギャップ充填_対象なし);
                return;
            }
            Logger.V_出力(メッセージID.ギャップ充填統計, p_統計.A_埋めたギャップ数, p_統計.A_総ギャップ数, p_統計.A_埋めた塩基数, p_統計.A_一意に定まらなかった数, p_統計.A_到達できなかった数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// ギャップの左右の足場から、その間を埋める配列を探す
        /// </summary>
        /// <param name="p_左側の出力"></param>
        /// <param name="p_配列"></param>
        /// <param name="p_ギャップ長"></param>
        /// <param name="p_ギャップ終端"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_判定"></param>
        /// <returns></returns>
        /// <remarks>
        /// 見つからない/一意に定まらない場合は null を返す
        /// </remarks>
        private static string? Get_ギャップを埋める配列(StringBuilder p_左側の出力, string p_配列, int p_ギャップ長, int p_ギャップ終端, TrustedKmerIndex p_kmerインデックス, int p_k長, out ギャップ充填判定 p_判定)
        {
            p_判定 = ギャップ充填判定.到達不能;

            if (p_ギャップ長 > Consts.ギャップ充填のギャップ長上限 || p_左側の出力.Length < p_k長)
            {
                return null;
            }

            if (p_ギャップ終端 + p_k長 > p_配列.Length)
            {
                return null;
            }

            // 左側の足場: 既に書き出した配列の末尾 k-mer
            var l_左のkmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                l_左のkmer[i] = Util.Get_塩基ID(p_左側の出力[p_左側の出力.Length - p_k長 + i]);
            }

            // 右側の足場: ギャップ直後の k-mer
            // ここへ到達できれば繋がったことになる
            var l_目標kmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                l_目標kmer[i] = Util.Get_塩基ID(p_配列[p_ギャップ終端 + i]);
            }

            if (Array.IndexOf(l_左のkmer, Consts.無効な塩基) >= 0 || Array.IndexOf(l_目標kmer, Consts.無効な塩基) >= 0)
            {
                return null;
            }

            if (!p_kmerインデックス.Get_含まれるか(l_左のkmer) || !p_kmerインデックス.Get_含まれるか(l_目標kmer))
            {
                // 足場そのものが信頼できる k-mer 集合に無いなら探索しても意味がない
                return null;
            }

            var l_最小長 = Math.Max(0, p_ギャップ長 - Consts.ギャップ充填の長さの余裕幅);
            var l_最大長 = p_ギャップ長 + Consts.ギャップ充填の長さの余裕幅;

            (var l_経路, p_判定) = ConstrainedPathFinder.Get_経路(l_左のkmer, l_目標kmer, l_最小長, l_最大長, p_kmerインデックス, p_k長);
            return l_経路;
        }

        #endregion
    }
}
