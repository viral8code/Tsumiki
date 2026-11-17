using System.Text;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 最終配列に残った N の列を、ほかの k の contig のうち両脇の錨を同じ向きで順に含むものの、錨の間の配列で埋める<br/>
    /// 埋める配列は、両脇の片と繋いでリードに続けて出てくる (PairBridgeFiller.Is写しがリードに続けて出る) ものがただ 1 つに決まるときだけ使う
    /// </summary>
    internal static class OtherKGapFiller
    {
        #region 定数

        /// <summary>
        /// N の列の両脇から取る錨の長さ
        /// </summary>
        private const int C_錨長 = 60;

        /// <summary>
        /// 埋める配列の長さの上限
        /// </summary>
        private const int C_埋める長さの上限 = 2000;

        /// <summary>
        /// 埋める配列をリードで確かめる窓の長さの、リード長に対する比の分子<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分子 = 7;

        /// <summary>
        /// 埋める配列をリードで確かめる窓の長さの、リード長に対する比の分母<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分母 = 10;

        /// <summary>
        /// 連結した配列の区切りの文字 (塩基にもギャップ文字にもない文字)
        /// </summary>
        private const char C_区切り = '#';

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群の N の列を、ほかの k の contig の配列で埋める
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_ほかのkの配列群">ほかの k の contig の配列</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_統計">埋めた結果の統計</param>
        /// <param name="p_埋めた番号">配列 ID ごとの、埋めた N の列の番号 (埋める前の配列での 0 始まり)</param>
        /// <returns>埋めた後の配列の ID と配列 (順と ID は入力のまま)</returns>
        public static List<(string A_ID, string A_配列)> Get_埋めた配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IReadOnlyList<string> p_ほかのkの配列群, ReadMinimizerIndex p_リード索引, int p_リード長, out ほかのkで埋めた統計 p_統計, out Dictionary<string, HashSet<int>> p_埋めた番号)
        {
            var l_全体 = string.Join(C_区切り, p_配列群.Select(x => x.A_配列));
            var l_ほかのk = string.Join(C_区切り, p_ほかのkの配列群);
            var l_窓長 = Math.Max(ReadMinimizerIndex.C_最短の問い合わせ長, p_リード長 * C_確かめる窓の分子 / C_確かめる窓の分母);
            p_埋めた番号 = [];
            int l_対象数 = 0, l_埋めた数 = 0, l_錨が取れない数 = 0, l_反復の端の数 = 0, l_無い数 = 0, l_リードに無い数 = 0;
            List<(string A_ID, string A_配列)> l_結果 = [];
            foreach (var (l_ID, l_配列) in p_配列群)
            {
                var l_出力 = new StringBuilder(l_配列.Length);
                var l_番号 = 0;
                var l_位置 = 0;
                while (l_位置 < l_配列.Length)
                {
                    var l_始め = l_配列.IndexOf(Consts.塩基コード.N, l_位置);
                    if (l_始め < 0)
                    {
                        _ = l_出力.Append(l_配列.AsSpan(l_位置));
                        break;
                    }

                    var l_終わり = l_始め;
                    while (l_終わり < l_配列.Length && l_配列[l_終わり] == Consts.塩基コード.N)
                    {
                        l_終わり++;
                    }

                    _ = l_出力.Append(l_配列.AsSpan(l_位置, l_始め - l_位置));
                    l_対象数++;
                    var l_埋める配列 = Get_埋める配列(l_配列, l_始め, l_終わり, l_全体, l_ほかのk, p_リード索引, l_窓長, ref l_錨が取れない数, ref l_反復の端の数, ref l_無い数, ref l_リードに無い数);
                    if (l_埋める配列 is not null)
                    {
                        l_埋めた数++;
                        Logger.V_出力(メッセージID.ほかのkで埋めたNの列の内訳, l_ID, l_始め, l_終わり - l_始め, l_埋める配列.Length, l_配列.Substring(l_始め - C_錨長, C_錨長), l_配列.Substring(l_終わり, C_錨長));
                        _ = l_出力.Append(l_埋める配列);
                        if (!p_埋めた番号.TryGetValue(l_ID, out var l_番号群))
                        {
                            l_番号群 = [];
                            p_埋めた番号[l_ID] = l_番号群;
                        }

                        _ = l_番号群.Add(l_番号);
                    }
                    else
                    {
                        _ = l_出力.Append(l_配列.AsSpan(l_始め, l_終わり - l_始め));
                    }

                    l_番号++;
                    l_位置 = l_終わり;
                }

                l_結果.Add((l_ID, l_出力.ToString()));
            }

            p_統計 = new ほかのkで埋めた統計(l_対象数, l_埋めた数, l_錨が取れない数, l_反復の端の数, l_無い数, l_リードに無い数);
            return l_結果;
        }

        /// <summary>
        /// 一部の N の列を埋めたあとの、長さ不明のギャップの N の列の番号を求める<br/>
        /// 埋めた番号を除き、それより後ろの番号は手前で埋めた数だけ詰める
        /// </summary>
        /// <param name="p_長さ不明の番号">埋める前の、長さ不明のギャップの N の列の番号</param>
        /// <param name="p_埋めた番号">埋めた N の列の番号 (埋める前の配列での番号)</param>
        /// <returns>埋めた後の番号</returns>
        public static HashSet<int> Get_埋めた後の番号(IReadOnlySet<int> p_長さ不明の番号, IReadOnlySet<int> p_埋めた番号)
        {
            return [.. p_長さ不明の番号.Where(x => !p_埋めた番号.Contains(x)).Select(x => x - p_埋めた番号.Count(y => y < x))];
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 つの N の列を埋める配列を求める<br/>
        /// 両脇の錨 (C_錨長) が A・C・G・T だけでアセンブリの中で 1 回だけで、<br/>
        /// 両脇の片 (隣の N の列か配列の端まで) がどちらも PairBridgeFiller.C_反復を調べる長さ 以上あり、PairBridgeFiller.Is反復の端 でないときだけ候補を探す<br/>
        /// ほかの k の連結配列から、左の錨のあと C_埋める長さの上限 以内に右の錨が続く所 (逆相補の向きも) を集め、錨の間の配列を候補にする<br/>
        /// 候補のうち、両脇の片と繋いでリードに続けて出てくるものがただ 1 つなら、それを返す
        /// </summary>
        /// <param name="p_配列">N の列を含む配列</param>
        /// <param name="p_始め">N の列の始まり</param>
        /// <param name="p_終わり">N の列の終わり (含まない)</param>
        /// <param name="p_全体">アセンブリの全配列を区切りで連結した配列</param>
        /// <param name="p_ほかのk">ほかの k の contig を区切りで連結した配列</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_窓長">リードで確かめる窓の長さ</param>
        /// <param name="p_錨が取れない数">錨が取れないか一意でない数 (加算する)</param>
        /// <param name="p_反復の端の数">両脇の片が短いか、反復の端である数 (加算する)</param>
        /// <param name="p_無い数">ほかの k に候補が無い数 (加算する)</param>
        /// <param name="p_リードに無い数">候補がリードでただ 1 つに決まらない数 (加算する)</param>
        /// <returns>埋める配列、埋めなければ null</returns>
        private static string? Get_埋める配列(string p_配列, int p_始め, int p_終わり, string p_全体, string p_ほかのk, ReadMinimizerIndex p_リード索引, int p_窓長, ref int p_錨が取れない数, ref int p_反復の端の数, ref int p_無い数, ref int p_リードに無い数)
        {
            if (p_始め < C_錨長 || p_終わり + C_錨長 > p_配列.Length)
            {
                p_錨が取れない数++;
                return null;
            }

            var l_左の錨 = p_配列.Substring(p_始め - C_錨長, C_錨長);
            var l_右の錨 = p_配列.Substring(p_終わり, C_錨長);
            if (!Is一意な錨(l_左の錨, p_全体) || !Is一意な錨(l_右の錨, p_全体))
            {
                p_錨が取れない数++;
                return null;
            }

            var l_窓 = PairBridgeFiller.C_反復を調べる長さ;
            var l_左の片の長さ = p_始め - (p_配列.LastIndexOf(Consts.塩基コード.N, p_始め - 1) + 1);
            var l_右の片の終わり = p_配列.IndexOf(Consts.塩基コード.N, p_終わり);
            var l_右の片の長さ = (l_右の片の終わり < 0 ? p_配列.Length : l_右の片の終わり) - p_終わり;
            if (l_左の片の長さ < l_窓 || l_右の片の長さ < l_窓 || PairBridgeFiller.Is反復の端(p_配列.Substring(p_始め - l_窓, l_窓), p_配列.Substring(p_終わり, l_窓), null))
            {
                p_反復の端の数++;
                return null;
            }

            HashSet<string> l_候補 = new(StringComparer.Ordinal);
            V_集める_錨の間(p_ほかのk, l_左の錨, l_右の錨, false, l_候補);
            V_集める_錨の間(p_ほかのk, Util.V_逆相補(l_右の錨), Util.V_逆相補(l_左の錨), true, l_候補);
            if (l_候補.Count == 0)
            {
                p_無い数++;
                return null;
            }

            var l_片の長さ = p_窓長 - 1;
            var l_左 = p_配列.Substring(Math.Max(0, p_始め - l_片の長さ), Math.Min(p_始め, l_片の長さ));
            var l_右 = p_配列.Substring(p_終わり, Math.Min(p_配列.Length - p_終わり, l_片の長さ));
            var l_通った = l_候補.Where(x => PairBridgeFiller.Is写しがリードに続けて出る(p_リード索引, l_左, x, l_右, p_窓長)).Take(2).ToList();
            if (l_通った.Count != 1)
            {
                p_リードに無い数++;
                return null;
            }

            return l_通った[0];
        }

        /// <summary>
        /// 錨が A・C・G・T だけで、連結配列の中に両方の向きを合わせて 1 回だけ現れるか
        /// </summary>
        /// <param name="p_錨">錨</param>
        /// <param name="p_全体">アセンブリの全配列を区切りで連結した配列</param>
        /// <returns>一意な錨なら真</returns>
        private static bool Is一意な錨(string p_錨, string p_全体)
        {
            if (p_錨.AsSpan().IndexOfAnyExcept(Consts.塩基文字) >= 0)
            {
                return false;
            }

            var l_逆 = Util.V_逆相補(p_錨);
            return Get_出現数(p_全体, p_錨, 2) + (l_逆 == p_錨 ? 0 : Get_出現数(p_全体, l_逆, 2)) == 1;
        }

        /// <summary>
        /// 連結配列の中の、部分配列の出現数を上限まで数える
        /// </summary>
        /// <param name="p_全体">連結配列</param>
        /// <param name="p_部分">部分配列</param>
        /// <param name="p_上限">ここまで数えたら打ち切る</param>
        /// <returns>出現数 (上限で頭打ち)</returns>
        private static int Get_出現数(string p_全体, string p_部分, int p_上限)
        {
            var l_数 = 0;
            for (var l_位置 = p_全体.IndexOf(p_部分, StringComparison.Ordinal); l_位置 >= 0 && l_数 < p_上限; l_位置 = p_全体.IndexOf(p_部分, l_位置 + 1, StringComparison.Ordinal))
            {
                l_数++;
            }

            return l_数;
        }

        /// <summary>
        /// ほかの k の連結配列で、左の錨のあと C_埋める長さの上限 以内に右の錨が続く所ごとに、錨の間の配列を候補に足す<br/>
        /// 間に区切りか A・C・G・T 以外を含むものは除く。逆相補の向きで探したときは、間の配列を逆相補にして足す
        /// </summary>
        /// <param name="p_ほかのk">ほかの k の contig を区切りで連結した配列</param>
        /// <param name="p_左">探す左の錨</param>
        /// <param name="p_右">探す右の錨</param>
        /// <param name="p_Is逆相補">逆相補の向きで探したか</param>
        /// <param name="p_候補">候補を足す先</param>
        internal static void V_集める_錨の間(string p_ほかのk, string p_左, string p_右, bool p_Is逆相補, HashSet<string> p_候補)
        {
            for (var l_左の位置 = p_ほかのk.IndexOf(p_左, StringComparison.Ordinal); l_左の位置 >= 0; l_左の位置 = p_ほかのk.IndexOf(p_左, l_左の位置 + 1, StringComparison.Ordinal))
            {
                var l_間の始め = l_左の位置 + p_左.Length;
                var l_探す長さ = Math.Min(p_ほかのk.Length - l_間の始め, C_埋める長さの上限 + p_右.Length);
                var l_右の位置 = p_ほかのk.IndexOf(p_右, l_間の始め, l_探す長さ, StringComparison.Ordinal);
                if (l_右の位置 < 0)
                {
                    continue;
                }

                var l_間 = p_ほかのk.AsSpan(l_間の始め, l_右の位置 - l_間の始め);
                if (l_間.IndexOfAnyExcept(Consts.塩基文字) >= 0)
                {
                    continue;
                }

                _ = p_候補.Add(p_Is逆相補 ? Util.V_逆相補(l_間.ToString()) : l_間.ToString());
            }
        }

        #endregion
    }
}
