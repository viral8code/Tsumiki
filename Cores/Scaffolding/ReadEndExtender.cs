using System.Text;
using Tsumiki.Commons;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 配列の端を、端の手前の錨を含むリードの続きの多数決で伸ばす<br/>
    /// 錨は端から C_置き直す長さ だけ手前で終わる配列。錨の続きを 1 塩基ずつ、リードの合意がはっきりした所まで取り、端の末尾 C_置き直す長さ 塩基を除いてそれを足す<br/>
    /// 配列の右端を伸ばし、左端は逆相補の右端として同じ処理をして戻す
    /// </summary>
    internal static class ReadEndExtender
    {
        #region 定数

        /// <summary>
        /// 端の錨の長さ (ReadMinimizerIndex の最短の問い合わせ長)
        /// </summary>
        private const int C_錨長 = ReadMinimizerIndex.C_最短の問い合わせ長;

        /// <summary>
        /// 伸ばす対象の配列の最短の長さ
        /// </summary>
        private const int C_対象の最短長 = 500;

        /// <summary>
        /// 端の片がこれより短ければ伸ばさない
        /// </summary>
        private const int C_端の片の最短長 = 200;

        /// <summary>
        /// 端から手前にずらした錨の位置、および置き直す塩基数 (末尾からこの長さを除いて合意を足す)
        /// </summary>
        private const int C_置き直す長さ = 30;

        /// <summary>
        /// Get_前後群 に渡す、続きを見るリード数の上限
        /// </summary>
        private const int C_続きを見るリード数の上限 = 64;

        /// <summary>
        /// 錨の手前を比べる長さ (リードの手前が端の片と同じ場所か)
        /// </summary>
        private const int C_手前を比べる長さ = 20;

        /// <summary>
        /// リードの続きの末尾から使わない塩基数 (リードの末尾は読み違いが多いので、多数決に入れない)
        /// </summary>
        private const int C_使わない末尾の長さ = 20;

        /// <summary>
        /// 続きの位置ごとに、その位置まで続きがあるリードの数の下限
        /// </summary>
        private const int C_要るリード数 = 8;

        /// <summary>
        /// 最多の塩基が、その位置まで続きがあるリードのうちで占める割合の下限 (これ以上なら合意とみなす)
        /// </summary>
        private const double C_合意とみなす割合 = 0.95D;

        /// <summary>
        /// 片側の端を伸ばす処理を繰り返す回数の上限 (伸ばした新しい端から、同じ条件でもう一度伸ばす)
        /// </summary>
        private const int C_繰り返す回数 = 3;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 長さが足りる配列の両端を、リードの続きの合意で伸ばす
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安 (Scaffolder.Get_一意の出現数) 、分からなければ 0</param>
        /// <param name="p_伸ばした端の数">伸ばした端の数</param>
        /// <param name="p_伸ばした塩基数">伸ばした塩基数 (合意の長さ − 置き直す長さ の合計)</param>
        /// <returns>伸ばした後の配列の ID と配列 (順と ID は入力のまま)</returns>
        public static List<(string A_ID, string A_配列)> Get_伸ばした配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, ReadMinimizerIndex p_リード索引, int p_リード長, int p_一意の出現数, out int p_伸ばした端の数, out int p_伸ばした塩基数)
        {
            p_伸ばした端の数 = 0;
            p_伸ばした塩基数 = 0;
            List<(string A_ID, string A_配列)> l_結果 = [];
            foreach (var (l_ID, l_配列) in p_配列群)
            {
                if (l_配列.Length < C_対象の最短長 || l_ID.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase))
                {
                    l_結果.Add((l_ID, l_配列));
                    continue;
                }

                var l_右 = Get_右端を繰り返し伸ばした(l_配列, p_リード索引, p_リード長, p_一意の出現数, out var l_右の塩基数);
                var l_左 = Get_右端を繰り返し伸ばした(Util.V_逆相補_曖昧塩基あり(l_右), p_リード索引, p_リード長, p_一意の出現数, out var l_左の塩基数);
                var l_伸ばした配列 = l_左の塩基数 > 0 ? Util.V_逆相補_曖昧塩基あり(l_左) : l_右;
                p_伸ばした端の数 += (l_右の塩基数 > 0 ? 1 : 0) + (l_左の塩基数 > 0 ? 1 : 0);
                p_伸ばした塩基数 += l_右の塩基数 + l_左の塩基数;
                l_結果.Add((l_ID, l_伸ばした配列));
            }

            return l_結果;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列の右端を、伸びなくなるか C_繰り返す回数 回に達するまで、伸ばした配列に対して繰り返し伸ばす
        /// </summary>
        /// <param name="p_配列">向き付きの配列</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安、分からなければ 0</param>
        /// <param name="p_伸ばした塩基数">すべての回で伸ばした塩基数の合計、伸ばさなければ 0</param>
        /// <returns>伸ばした後の配列</returns>
        private static string Get_右端を繰り返し伸ばした(string p_配列, ReadMinimizerIndex p_リード索引, int p_リード長, int p_一意の出現数, out int p_伸ばした塩基数)
        {
            p_伸ばした塩基数 = 0;
            var l_配列 = p_配列;
            for (var i = 0; i < C_繰り返す回数; i++)
            {
                l_配列 = Get_右端を伸ばした(l_配列, p_リード索引, p_リード長, p_一意の出現数, out var l_今回の塩基数);
                if (l_今回の塩基数 == 0)
                {
                    break;
                }

                p_伸ばした塩基数 += l_今回の塩基数;
            }

            return l_配列;
        }

        /// <summary>
        /// 配列の右端 (最後のギャップ文字より後ろの片) を、錨の続きの合意で伸ばす<br/>
        /// 片が短い、錨が反復の中、合意が置き直す長さ以下なら、伸ばさずに元の配列を返す
        /// </summary>
        /// <param name="p_配列">向き付きの配列</param>
        /// <param name="p_リード索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安、分からなければ 0</param>
        /// <param name="p_伸ばした塩基数">伸ばした塩基数、伸ばさなければ 0</param>
        /// <returns>伸ばした後の配列</returns>
        private static string Get_右端を伸ばした(string p_配列, ReadMinimizerIndex p_リード索引, int p_リード長, int p_一意の出現数, out int p_伸ばした塩基数)
        {
            p_伸ばした塩基数 = 0;
            var l_片の始め = 0;
            for (var i = p_配列.Length - 1; i >= 0; i--)
            {
                if (Util.Isギャップ文字(p_配列[i]))
                {
                    l_片の始め = i + 1;
                    break;
                }
            }

            var l_片 = p_配列[l_片の始め..];
            if (l_片.Length < C_端の片の最短長)
            {
                return p_配列;
            }

            var l_錨の始め = l_片.Length - C_置き直す長さ - C_錨長;
            var l_錨 = l_片.Substring(l_錨の始め, C_錨長);
            if (l_錨.AsSpan().IndexOfAnyExcept(Consts.塩基文字) >= 0 || Scaffolder.Is反復の中(p_リード索引, l_片, l_錨の始め, p_一意の出現数))
            {
                return p_配列;
            }

            var l_手前 = l_片[..l_錨の始め];
            var l_続き群 = p_リード索引.Get_前後群(l_錨, p_リード長, C_続きを見るリード数の上限).Where(x => Scaffolder.Is手前が同じ場所(l_手前, x.A_前, C_手前を比べる長さ)).Select(x => x.A_続き).Where(x => x.Length > C_使わない末尾の長さ).Select(x => x[..^C_使わない末尾の長さ]).ToList();
            var l_合意 = Get_合意(l_続き群);
            if (l_合意.Length <= C_置き直す長さ)
            {
                return p_配列;
            }

            p_伸ばした塩基数 = l_合意.Length - C_置き直す長さ;
            return string.Concat(p_配列.AsSpan(0, p_配列.Length - C_置き直す長さ), l_合意);
        }

        /// <summary>
        /// 錨の続きを先頭から 1 塩基ずつ見て、その位置まで続きがあるリードが 要るリード数 本以上あり、最多の塩基の割合が 合意とみなす割合 以上の間だけ、その塩基を取る
        /// </summary>
        /// <param name="p_続き群">錨の直後の塩基列 (リードの端で止まる)</param>
        /// <returns>合意の塩基列</returns>
        private static string Get_合意(List<string> p_続き群)
        {
            var l_合意 = new StringBuilder();
            for (var j = 0; ; j++)
            {
                var l_数 = new int[Consts.塩基文字.Length];
                var l_候補数 = 0;
                foreach (var l_続き in p_続き群)
                {
                    if (l_続き.Length > j)
                    {
                        l_候補数++;
                        var l_番号 = Consts.塩基文字.IndexOf(l_続き[j]);
                        if (l_番号 >= 0)
                        {
                            l_数[l_番号]++;
                        }
                    }
                }

                if (l_候補数 < C_要るリード数)
                {
                    break;
                }

                var l_最多 = 0;
                for (var k = 1; k < l_数.Length; k++)
                {
                    if (l_数[k] > l_数[l_最多])
                    {
                        l_最多 = k;
                    }
                }

                if ((double)l_数[l_最多] / l_候補数 < C_合意とみなす割合)
                {
                    break;
                }

                l_合意.Append(Consts.塩基文字[l_最多]);
            }

            return l_合意.ToString();
        }

        #endregion
    }
}
