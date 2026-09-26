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
    /// scaffold の N を、グラフ上で両側を繋ぐ経路を探して実配列に置き換える
    /// </summary>
    internal static class GapFiller
    {
        #region 定数

        /// <summary>
        /// 経路上の各 k-mer に求める最小カバレッジ
        /// </summary>
        private const ulong C_経路の最小カバレッジ = 2UL;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// scaffold を読み込み、埋められるギャップを埋めて同じパスへ書き戻す
        /// </summary>
        /// <param name="p_scaffoldパス"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        public static ギャップ充填統計 V_充填_ギャップ(string p_scaffoldパス, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_scaffold群 = FastaReader.Get_全エントリ(p_scaffoldパス);

            var l_総ギャップ数 = 0;
            var l_埋めたギャップ数 = 0;
            var l_埋めた塩基数 = 0;
            var l_一意でない数 = 0;
            var l_到達不能数 = 0;

            List<(string A_ID, string A_配列)> l_結果 = [];
            foreach (var (l_ID, l_配列) in l_scaffold群)
            {
                var l_出力 = new StringBuilder();
                var l_位置 = 0;
                while (l_位置 < l_配列.Length)
                {
                    if (!Util.Isギャップ文字(l_配列[l_位置]))
                    {
                        _ = l_出力.Append(l_配列[l_位置]);
                        l_位置++;
                        continue;
                    }

                    var l_ギャップ開始 = l_位置;
                    l_位置 = Util.Get_ギャップの終わり(l_配列, l_位置);

                    var l_ギャップ長 = l_位置 - l_ギャップ開始;
                    l_総ギャップ数++;

                    var l_埋めた配列 = Get_ギャップ充填配列(l_出力, l_配列, l_ギャップ長, l_位置, p_kmerインデックス, p_k長, out var l_判定, out var l_Isアンカー不足, out var l_Is支持不足, out var l_安定ID);
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

                        var l_種別 = l_Isアンカー不足 ? 曖昧箇所の種別.アンカー不足
                            : l_Is支持不足 ? 曖昧箇所の種別.支持なし
                            : l_判定 switch
                            {
                                ギャップ充填判定.一意でない => 曖昧箇所の種別.経路が一意でない,
                                ギャップ充填判定.探索打切り => 曖昧箇所の種別.探索打切り,
                                _ => 曖昧箇所の種別.到達不能,
                            };
                        AmbiguityRecorder.V_記録(l_種別, $"{l_ID.TrimStart('>')}:{l_ギャップ開始}-{l_位置}", p_安定ID: l_安定ID);
                        _ = l_出力.Append(l_配列, l_ギャップ開始, l_ギャップ長);
                    }
                }

                l_結果.Add((l_ID, l_出力.ToString()));
            }

            using (var l_書き込み = new FastaWriter(p_scaffoldパス))
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
        /// <param name="p_Isアンカー不足">足場そのものが信頼できる k-mer 集合に無く、探索を始められなかった場合 true</param>
        /// <param name="p_Is支持不足">経路自体は一意に見つかったが、経路上の k-mer カバレッジが薄く採用を見送った場合 true</param>
        /// <param name="p_安定ID">この箇所を再実行をまたいで追跡するための安定ID</param>
        /// <returns></returns>
        private static string? Get_ギャップ充填配列(StringBuilder p_左側の出力, string p_配列, int p_ギャップ長, int p_ギャップ終端, TrustedKmerIndex p_kmerインデックス, int p_k長, out ギャップ充填判定 p_判定, out bool p_Isアンカー不足, out bool p_Is支持不足, out string p_安定ID)
        {
            p_判定 = ギャップ充填判定.到達不能;
            p_Isアンカー不足 = false;
            p_Is支持不足 = false;
            p_安定ID = string.Empty;

            if (p_ギャップ長 > Consts.ギャップ充填のギャップ長上限 || p_左側の出力.Length < p_k長)
            {
                return null;
            }

            if (p_ギャップ終端 + p_k長 > p_配列.Length)
            {
                return null;
            }

            var l_左のkmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                l_左のkmer[i] = Util.Get_塩基ID(p_左側の出力[p_左側の出力.Length - p_k長 + i]);
            }

            var l_目標kmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                l_目標kmer[i] = Util.Get_塩基ID(p_配列[p_ギャップ終端 + i]);
            }

            if (Array.IndexOf(l_左のkmer, Consts.無効な塩基) >= 0 || Array.IndexOf(l_目標kmer, Consts.無効な塩基) >= 0)
            {
                return null;
            }

            p_安定ID = AmbiguityRecorder.Get_安定ID(new string([.. l_左のkmer.Select(Util.Get_塩基文字)]), new string([.. l_目標kmer.Select(Util.Get_塩基文字)]));

            if (!p_kmerインデックス.Haskmer(l_左のkmer) || !p_kmerインデックス.Haskmer(l_目標kmer))
            {
                p_Isアンカー不足 = true;
                return null;
            }

            var l_最小長 = Math.Max(0, p_ギャップ長 - Consts.ギャップ充填の長さの余裕幅);
            var l_最大長 = p_ギャップ長 + Consts.ギャップ充填の長さの余裕幅;

            (var l_経路, p_判定) = ConstrainedPathFinder.Get_経路(l_左のkmer, l_目標kmer, l_最小長, l_最大長, p_kmerインデックス, p_k長);
            if (l_経路 is null)
            {
                return null;
            }

            if (!Has十分なカバレッジ支持(p_左側の出力, l_経路, p_配列, p_ギャップ終端, p_kmerインデックス, p_k長))
            {
                p_Is支持不足 = true;
                return null;
            }

            return l_経路;
        }

        /// <summary>
        /// 経路上のすべての k-mer が、最低限のカバレッジを持つことを確かめる
        /// </summary>
        /// <param name="p_左側の出力">既に書き出した配列</param>
        /// <param name="p_経路">見つかった充填配列</param>
        /// <param name="p_配列">元の scaffold 配列</param>
        /// <param name="p_ギャップ終端">ギャップの終端位置</param>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>すべての k-mer が最小カバレッジ以上なら true</returns>
        private static bool Has十分なカバレッジ支持(StringBuilder p_左側の出力, string p_経路, string p_配列, int p_ギャップ終端, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_接続 = p_左側の出力.ToString(p_左側の出力.Length - p_k長, p_k長) + p_経路 + p_配列.Substring(p_ギャップ終端, p_k長);
            for (var i = 0; i + p_k長 <= l_接続.Length; i++)
            {
                var l_kmer = Get_kmerバイト列(l_接続, i, p_k長);
                if (l_kmer is null || p_kmerインデックス.Get_カバレッジ(l_kmer) < C_経路の最小カバレッジ)
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// 指定した位置の k-mer を塩基 ID 列にして返す
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_開始">k-mer の開始位置</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>塩基 ID 列、曖昧塩基を含む場合は null</returns>
        private static byte[]? Get_kmerバイト列(string p_配列, int p_開始, int p_k長)
        {
            var l_kmer = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_配列[p_開始 + i]);
                if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    return null;
                }

                l_kmer[i] = l_塩基ID;
            }

            return l_kmer;
        }

        #endregion
    }
}
