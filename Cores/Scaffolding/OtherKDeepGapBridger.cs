using System.Text;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 最終配列に残った N の列を、両脇の片の端から離れた深い錨と、ほかの k の contig の錨の間の配列で繋ぐ<br/>
    /// 左の錨から右の錨までがリード 1 本に収まり、その配列がリードに出てくるときだけ繋ぐ (両脇の片が重なる形なら畳む)
    /// </summary>
    internal static class OtherKDeepGapBridger
    {
        #region 定数

        /// <summary>
        /// 繋ぐ形をリードで確かめる窓の長さの、リード長に対する比の分子<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分子 = 7;

        /// <summary>
        /// 繋ぐ形をリードで確かめる窓の長さの、リード長に対する比の分母<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分母 = 10;

        /// <summary>
        /// 左の錨から右の錨までをリードで確かめるときに要る、リードの出現数の下限<br/>
        /// PairBridgeFiller.Is写しがリードに続けて出る の窓の支持数と同じ
        /// </summary>
        private const int C_支持数の下限 = 2;

        /// <summary>
        /// 繋いだ所のログに出す、左の錨までと右の錨からの配列の長さ<br/>
        /// 評価のあとで、最終配列の中の場所を探すのに使う
        /// </summary>
        private const int C_ログの脇の長さ = 60;

        /// <summary>
        /// ほかの k の contig を連結するときの区切りの文字 (塩基にもギャップ文字にもない文字)
        /// </summary>
        private const char C_区切り = '#';

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群の N の列を、深い錨とほかの k の contig の配列で繋ぐ<br/>
        /// N の列ごとに、左の片 (手前の N の列か配列の端から。手前の N の列を繋いだときは、そこで使った右の錨から) と右の片 (次の N の列か配列の端まで) を OtherKEdgeBridger.Get_繋ぐ形 に渡し、<br/>
        /// 左の錨から右の錨までが Is錨から錨までがリードにある を満たすときだけ、左の片を左の錨まで、右の片を右の錨から残して、錨の間の配列で繋ぐ
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_ほかのkの配列群">ほかの k の contig の配列</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_統計">繋いだ結果の統計</param>
        /// <param name="p_繋いだ番号">配列 ID ごとの、繋いだ N の列の番号 (繋ぐ前の配列での 0 始まり)</param>
        /// <returns>繋いだ後の配列の ID と配列 (順と ID は入力のまま)</returns>
        public static List<(string A_ID, string A_配列)> Get_繋いだ配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IReadOnlyList<string> p_ほかのkの配列群, ReadMinimizerIndex p_リード索引, int p_リード長, out 深い錨で繋いだ統計 p_統計, out Dictionary<string, HashSet<int>> p_繋いだ番号)
        {
            var l_窓長 = Math.Max(ReadMinimizerIndex.C_最短の問い合わせ長, p_リード長 * C_確かめる窓の分子 / C_確かめる窓の分母);
            var l_錨表 = OtherKRepeatReplacer.Get_一意な錨の表(p_配列群);
            var l_ほかのk = string.Join(C_区切り, p_ほかのkの配列群);
            var l_錨長 = OtherKRepeatReplacer.C_錨長;
            p_繋いだ番号 = [];
            int l_対象数 = 0, l_繋いだ数 = 0, l_畳んだ数 = 0, l_決まらない数 = 0, l_確かめられない数 = 0;
            List<(string A_ID, string A_配列)> l_結果 = [];
            foreach (var (l_ID, l_配列) in p_配列群)
            {
                var l_列群 = OverlapGapFolder.Get_Nの列群(l_配列);
                var l_出力 = new StringBuilder(l_配列.Length);
                var l_次の始め = 0;
                for (var l_番号 = 0; l_番号 < l_列群.Count; l_番号++)
                {
                    var (l_始め, l_終わり) = l_列群[l_番号];
                    var l_左 = l_配列[l_次の始め..l_始め];
                    var l_右 = l_配列[l_終わり..(l_番号 + 1 < l_列群.Count ? l_列群[l_番号 + 1].A_始め : l_配列.Length)];
                    l_対象数++;
                    var l_形 = OtherKEdgeBridger.Get_繋ぐ形(l_左, l_右, l_錨表, l_ほかのk, p_リード索引, p_リード長);
                    if (l_形 is not { } l_決まった形)
                    {
                        l_決まらない数++;
                        _ = l_出力.Append(l_配列.AsSpan(l_次の始め, l_終わり - l_次の始め));
                        l_次の始め = l_終わり;
                        continue;
                    }

                    var l_左の片 = l_左[..(l_決まった形.A_左の錨 + l_錨長)];
                    var l_錨から錨まで = string.Concat(l_左.AsSpan(l_決まった形.A_左の錨, l_錨長), l_決まった形.A_埋める配列, l_右.AsSpan(l_決まった形.A_右の錨, l_錨長));
                    if (!Is錨から錨までがリードにある(l_錨から錨まで, p_リード索引, p_リード長, l_窓長))
                    {
                        l_確かめられない数++;
                        _ = l_出力.Append(l_配列.AsSpan(l_次の始め, l_終わり - l_次の始め));
                        l_次の始め = l_終わり;
                        continue;
                    }

                    var l_差 = l_決まった形.A_埋める配列.Length - (l_決まった形.A_u側 + l_決まった形.A_x側);
                    l_繋いだ数++;
                    l_畳んだ数 += l_差 < 0 ? 1 : 0;
                    Logger.V_出力(メッセージID.深い錨で繋いだNの列の内訳, l_ID, l_始め, l_終わり - l_始め, l_決まった形.A_埋める配列.Length, l_差, l_左の片[^Math.Min(l_左の片.Length, C_ログの脇の長さ)..], l_右.Substring(l_決まった形.A_右の錨, Math.Min(l_右.Length - l_決まった形.A_右の錨, C_ログの脇の長さ)));
                    _ = l_出力.Append(l_左の片).Append(l_決まった形.A_埋める配列);
                    if (!p_繋いだ番号.TryGetValue(l_ID, out var l_番号群))
                    {
                        l_番号群 = [];
                        p_繋いだ番号[l_ID] = l_番号群;
                    }

                    _ = l_番号群.Add(l_番号);
                    l_次の始め = l_終わり + l_決まった形.A_右の錨;
                }

                _ = l_出力.Append(l_配列.AsSpan(l_次の始め));
                l_結果.Add((l_ID, l_出力.ToString()));
            }

            p_統計 = new 深い錨で繋いだ統計(l_対象数, l_繋いだ数, l_畳んだ数, l_決まらない数, l_確かめられない数);
            return l_結果;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 左の錨から右の錨まで (左の錨 + 錨の間の配列 + 右の錨) をリードで確かめられるか<br/>
        /// 窓長以下なら、OtherKEdgeBridger.Get_繋ぐ形 の窓の確かめで済んでいるので真。リード長以下ならリードに C_支持数の下限 か所以上あれば真、リード長より長ければ偽<br/>
        /// 錨の間が窓より長いと、縦の反復の単位の数や、参照で 2 か所にある錨のどちらに続くかを、窓では見分けられない
        /// </summary>
        /// <param name="p_錨から錨まで">左の錨 + 錨の間の配列 + 右の錨</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_窓長">リードで確かめる窓の長さ</param>
        /// <returns>確かめられれば真</returns>
        internal static bool Is錨から錨までがリードにある(string p_錨から錨まで, ReadMinimizerIndex p_リード索引, int p_リード長, int p_窓長)
        {
            return p_錨から錨まで.Length <= p_窓長
                || (p_錨から錨まで.Length <= p_リード長 && p_リード索引.Get_出現数(p_錨から錨まで, C_支持数の下限) >= C_支持数の下限);
        }

        #endregion
    }
}
