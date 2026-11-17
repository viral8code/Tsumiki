using System.Text;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 最終配列に残った N の列のうち、左の片の末尾と右の片の先頭が完全に一致する所を、リードに確かめて畳む<br/>
    /// 畳んだ形 (左の片 + 右の片から重なりを除いたもの) がリードに続けて出てくる重なりの長さが、ただ 1 つに決まるときだけ畳む
    /// </summary>
    internal static class OverlapGapFolder
    {
        #region 定数

        /// <summary>
        /// 畳む重なりの長さの下限<br/>
        /// これより短い一致は、偶然の一致と区別しにくい
        /// </summary>
        private const int C_重なりの下限 = 20;

        /// <summary>
        /// 畳む重なりの長さの上限
        /// </summary>
        private const int C_重なりの上限 = 2000;

        /// <summary>
        /// 畳んだ形をリードで確かめる窓の長さの、リード長に対する比の分子<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分子 = 7;

        /// <summary>
        /// 畳んだ形をリードで確かめる窓の長さの、リード長に対する比の分母<br/>
        /// 窓はリード長 × C_確かめる窓の分子 / C_確かめる窓の分母
        /// </summary>
        private const int C_確かめる窓の分母 = 10;

        /// <summary>
        /// 重なりの両脇をリードで確かめるときに要る、リードの出現数の下限<br/>
        /// PairBridgeFiller.Is写しがリードに続けて出る の窓の支持数と同じ
        /// </summary>
        private const int C_支持数の下限 = 2;

        /// <summary>
        /// 畳んだ所のログに出す、重なりの外の両脇の配列の長さ<br/>
        /// 評価のあとで、最終配列の中の場所を探すのに使う
        /// </summary>
        private const int C_ログの脇の長さ = 60;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群の N の列のうち、両脇の片が完全に重なる所を畳む<br/>
        /// 両脇の片 (隣の N の列か配列の端まで。手前の N の列を畳んだときは、その重なりを除く) がどちらも PairBridgeFiller.C_反復を調べる長さ 以上あり、<br/>
        /// 重なりを除いて反復の端でなく、畳んだ形が PairBridgeFiller.Is写しがリードに続けて出る を満たす重なりの長さが、ただ 1 つのときだけ畳む
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_統計">畳んだ結果の統計</param>
        /// <param name="p_畳んだ番号">配列 ID ごとの、畳んだ N の列の番号 (畳む前の配列での 0 始まり)</param>
        /// <returns>畳んだ後の配列の ID と配列 (順と ID は入力のまま)</returns>
        public static List<(string A_ID, string A_配列)> Get_畳んだ配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, ReadMinimizerIndex p_リード索引, int p_リード長, out Nの列を畳んだ統計 p_統計, out Dictionary<string, HashSet<int>> p_畳んだ番号)
        {
            var l_窓長 = Math.Max(ReadMinimizerIndex.C_最短の問い合わせ長, p_リード長 * C_確かめる窓の分子 / C_確かめる窓の分母);
            p_畳んだ番号 = [];
            int l_対象数 = 0, l_畳んだ数 = 0, l_片が短い数 = 0, l_重なりが無い数 = 0, l_反復の端の数 = 0, l_リードに無い数 = 0, l_決まらない数 = 0;
            List<(string A_ID, string A_配列)> l_結果 = [];
            foreach (var (l_ID, l_配列) in p_配列群)
            {
                var l_列群 = Get_Nの列群(l_配列);
                var l_出力 = new StringBuilder(l_配列.Length);
                var l_次の始め = 0;
                for (var l_番号 = 0; l_番号 < l_列群.Count; l_番号++)
                {
                    var (l_始め, l_終わり) = l_列群[l_番号];
                    var l_左 = l_配列[Math.Min(l_次の始め, l_始め)..l_始め];
                    var l_右 = l_配列[l_終わり..(l_番号 + 1 < l_列群.Count ? l_列群[l_番号 + 1].A_始め : l_配列.Length)];
                    l_対象数++;
                    var l_重なり = Get_畳む重なり(l_左, l_右, p_リード索引, p_リード長, l_窓長, ref l_片が短い数, ref l_重なりが無い数, ref l_反復の端の数, ref l_リードに無い数, ref l_決まらない数);
                    _ = l_出力.Append(l_配列.AsSpan(l_次の始め, l_始め - l_次の始め));
                    if (l_重なり is not { } l_o)
                    {
                        _ = l_出力.Append(l_配列.AsSpan(l_始め, l_終わり - l_始め));
                        l_次の始め = l_終わり;
                        continue;
                    }

                    l_畳んだ数++;
                    Logger.V_出力(メッセージID.重なりで畳んだNの列の内訳, l_ID, l_始め, l_終わり - l_始め, l_o, l_左[^Math.Min(l_左.Length, l_o + C_ログの脇の長さ)..^l_o], l_右[l_o..Math.Min(l_右.Length, l_o + C_ログの脇の長さ)]);
                    if (!p_畳んだ番号.TryGetValue(l_ID, out var l_番号群))
                    {
                        l_番号群 = [];
                        p_畳んだ番号[l_ID] = l_番号群;
                    }

                    _ = l_番号群.Add(l_番号);
                    l_次の始め = l_終わり + l_o;
                }

                _ = l_出力.Append(l_配列.AsSpan(l_次の始め));
                l_結果.Add((l_ID, l_出力.ToString()));
            }

            p_統計 = new Nの列を畳んだ統計(l_対象数, l_畳んだ数, l_片が短い数, l_重なりが無い数, l_反復の端の数, l_リードに無い数, l_決まらない数);
            return l_結果;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列の N の列の始めと終わり (含まない) を前から並べる
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <returns>N の列の始めと終わり</returns>
        internal static List<(int A_始め, int A_終わり)> Get_Nの列群(string p_配列)
        {
            List<(int A_始め, int A_終わり)> l_列群 = [];
            for (var l_始め = p_配列.IndexOf(Consts.塩基コード.N); l_始め >= 0; l_始め = p_配列.IndexOf(Consts.塩基コード.N, l_列群[^1].A_終わり))
            {
                var l_終わり = l_始め;
                while (l_終わり < p_配列.Length && p_配列[l_終わり] == Consts.塩基コード.N)
                {
                    l_終わり++;
                }

                l_列群.Add((l_始め, l_終わり));
            }

            return l_列群;
        }

        /// <summary>
        /// 1 つの N の列を畳む重なりの長さを求める<br/>
        /// 両脇の片がどちらも PairBridgeFiller.C_反復を調べる長さ 以上あるとき、左の片の末尾と右の片の先頭が完全に一致する長さ (C_重なりの下限 から C_重なりの上限) のうち、<br/>
        /// 重なりを除いて PairBridgeFiller.Is反復の端 でなく、畳んだ形が PairBridgeFiller.Is写しがリードに続けて出る と Is重なりの両脇がリードで続く を満たすものが、ただ 1 つならそれを返す
        /// </summary>
        /// <param name="p_左">左の片 (隣の N の列か配列の端から、N の列の前まで)</param>
        /// <param name="p_右">右の片 (N の列のあとから、隣の N の列か配列の端まで)</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_窓長">リードで確かめる窓の長さ</param>
        /// <param name="p_片が短い数">両脇の片のどちらかが反復を調べる長さより短い数 (加算する)</param>
        /// <param name="p_重なりが無い数">両脇の片が完全に一致する重なりが無い数 (加算する)</param>
        /// <param name="p_反復の端の数">一致する重なりはあるが、どれも重なりを除いて反復の端である数 (加算する)</param>
        /// <param name="p_リードに無い数">反復の端でない重なりはあるが、畳んだ形がどれもリードに続けて出てこない数 (加算する)</param>
        /// <param name="p_決まらない数">畳んだ形がリードに続けて出てくる重なりが 2 つ以上ある数 (加算する)</param>
        /// <returns>畳む重なりの長さ、畳まなければ null</returns>
        private static int? Get_畳む重なり(string p_左, string p_右, ReadMinimizerIndex p_リード索引, int p_リード長, int p_窓長, ref int p_片が短い数, ref int p_重なりが無い数, ref int p_反復の端の数, ref int p_リードに無い数, ref int p_決まらない数)
        {
            if (p_左.Length < PairBridgeFiller.C_反復を調べる長さ || p_右.Length < PairBridgeFiller.C_反復を調べる長さ)
            {
                p_片が短い数++;
                return null;
            }

            var l_一致群 = Get_一致する重なり群(p_左, p_右);
            if (l_一致群.Count == 0)
            {
                p_重なりが無い数++;
                return null;
            }

            var l_反復でない群 = l_一致群.Where(x => !PairBridgeFiller.Is反復の端(p_左, p_右, x)).ToList();
            if (l_反復でない群.Count == 0)
            {
                p_反復の端の数++;
                return null;
            }

            var l_通った = l_反復でない群.Where(x => PairBridgeFiller.Is写しがリードに続けて出る(p_リード索引, p_左[..^x], p_左[^x..], p_右[x..], p_窓長) && Is重なりの両脇がリードで続く(p_左, p_右, x, p_リード索引, p_リード長, p_窓長)).Take(2).ToList();
            if (l_通った.Count == 0)
            {
                p_リードに無い数++;
                return null;
            }

            if (l_通った.Count > 1)
            {
                p_決まらない数++;
                return null;
            }

            return l_通った[0];
        }

        /// <summary>
        /// 重なりと、その外の左右 1 塩基ずつを繋いだ配列 (左の片の最後の 1 塩基 + 重なり + 右の片の重なりの次の 1 塩基) が、リードに C_支持数の下限 か所以上あるか<br/>
        /// この配列が窓長以下なら、窓の確かめで済んでいるので真。リード長より長ければ確かめられないので偽<br/>
        /// 重なりが窓より長いと、窓は重なりの両脇を一度に含まないので、縦の反復の単位を抜いて畳んだ形でも窓がどれもリードに出てくる
        /// </summary>
        /// <param name="p_左">左の片</param>
        /// <param name="p_右">右の片</param>
        /// <param name="p_重なり">重なりの長さ</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_窓長">リードで確かめる窓の長さ</param>
        /// <returns>確かめられれば真</returns>
        internal static bool Is重なりの両脇がリードで続く(string p_左, string p_右, int p_重なり, ReadMinimizerIndex p_リード索引, int p_リード長, int p_窓長)
        {
            var l_長さ = p_重なり + 2;
            return l_長さ <= p_窓長
                || (l_長さ <= p_リード長 && p_リード索引.Get_出現数(string.Concat(p_左.AsSpan(p_左.Length - p_重なり - 1), p_右.AsSpan(p_重なり, 1)), C_支持数の下限) >= C_支持数の下限);
        }

        /// <summary>
        /// 左の片の末尾と右の片の先頭が完全に一致する長さを、C_重なりの下限 から短い順に集める<br/>
        /// 一致する所は A・C・G・T だけのものに限る。長さの上限は C_重なりの上限 と両脇の片の長さの短いほう
        /// </summary>
        /// <param name="p_左">左の片</param>
        /// <param name="p_右">右の片</param>
        /// <returns>一致する重なりの長さ</returns>
        internal static List<int> Get_一致する重なり群(string p_左, string p_右)
        {
            List<int> l_一致群 = [];
            var l_上限 = Math.Min(C_重なりの上限, Math.Min(p_左.Length, p_右.Length));
            for (var l_長さ = C_重なりの下限; l_長さ <= l_上限; l_長さ++)
            {
                var l_重なり = p_右.AsSpan(0, l_長さ);
                if (p_左.AsSpan(p_左.Length - l_長さ).SequenceEqual(l_重なり) && l_重なり.IndexOfAnyExcept(Consts.塩基文字) < 0)
                {
                    l_一致群.Add(l_長さ);
                }
            }

            return l_一致群;
        }

        #endregion
    }
}
