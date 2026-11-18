using Tsumiki.Commons;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 同じ所の候補の形どうしを、リード長より少し短い窓で比べて、リードで説明できる形をただ 1 つ選ぶ<br/>
    /// ReadGraphGapFiller の道の選び方 (ID13) と、OtherKRepeatReplacer の最終配列の形とほかの k の形の比べ (ID14) で使う
    /// </summary>
    internal static class ReadWindowSelector
    {
        #region 定数

        /// <summary>
        /// 形どうしを比べる窓の長さの、リード長に対する百分率<br/>
        /// リード長と同じだと窓を丸ごと含むリードが少なく、窓にリードが無いことが当てにならない
        /// </summary>
        private const int C_窓の百分率 = 85;

        /// <summary>
        /// 百分率の分母
        /// </summary>
        private const int C_百分率の分母 = 100;

        /// <summary>
        /// ほかの形にだけある窓にリードが無いことで、その形を外すのに要る窓の数の下限
        /// </summary>
        private const int C_窓の下限数 = 10;

        /// <summary>
        /// 窓に要る、リードの出現数の下限
        /// </summary>
        private const int C_支持数の下限 = 2;

        /// <summary>
        /// 窓がほかの形にあるかを調べるとき、形の両側に付ける片の長さ
        /// </summary>
        internal const int C_広い形の脇の長さ = 2000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 形どうしを比べる窓の長さ (リード長 × C_窓の百分率 / 100。リードの索引に問い合わせられる最短の長さ以上)
        /// </summary>
        /// <param name="p_リード長">代表リード長</param>
        /// <returns>窓の長さ</returns>
        internal static int Get_窓長(int p_リード長)
        {
            return Math.Max(ReadMinimizerIndex.C_最短の問い合わせ長, p_リード長 * C_窓の百分率 / C_百分率の分母);
        }

        /// <summary>
        /// 形の群を窓で比べ、条件を満たす形がただ 1 つならその番号を返す。条件は、その形の窓がどれもリードに C_支持数の下限 か所以上あり、<br/>
        /// ほかのどの形についても、ほかの形の窓のうちその形の広い形に無いもの (ほかの形にだけある窓) がリードに 1 か所も無く、<br/>
        /// その窓が C_窓の下限数 個以上あるか、その形の窓のうちほかの形の広い形に無いものがあること
        /// </summary>
        /// <param name="p_形群">形 (片の端 リード長 bp + 候補の配列 + 片の端 リード長 bp)</param>
        /// <param name="p_広い形群">広い形 (片を C_広い形の脇の長さ まで付けた形。p_形群 と同じ順)</param>
        /// <param name="p_窓長">窓の長さ</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_条件を満たす数">条件を満たした形の数</param>
        /// <param name="p_相手の窓の最少">選んだ形に対して、ほかの形にだけある窓の数の最少 (選ばなければ 0)</param>
        /// <returns>選んだ形の番号 (選ばなければ null)</returns>
        internal static int? Get_選んだ番号(IReadOnlyList<string> p_形群, IReadOnlyList<string> p_広い形群, int p_窓長, ReadMinimizerIndex p_索引, out int p_条件を満たす数, out int p_相手の窓の最少)
        {
            Dictionary<string, int> l_数の表 = [];
            List<(int A_番号, int A_最少)> l_満たす形群 = [];
            for (var i = 0; i < p_形群.Count; i++)
            {
                if (Get_窓群(p_形群[i], p_窓長).Any(x => Get_窓の出現数(x, l_数の表, p_索引) < C_支持数の下限))
                {
                    continue;
                }

                var l_Is満たす = true;
                var l_最少 = int.MaxValue;
                for (var j = 0; j < p_形群.Count && l_Is満たす; j++)
                {
                    if (j == i)
                    {
                        continue;
                    }

                    string[] l_相手だけの窓群 = [.. Get_窓群(p_形群[j], p_窓長).Where(x => !p_広い形群[i].Contains(x, StringComparison.Ordinal))];
                    l_Is満たす = !l_相手だけの窓群.Any(x => Get_窓の出現数(x, l_数の表, p_索引) > 0)
                        && (l_相手だけの窓群.Length >= C_窓の下限数 || Get_窓群(p_形群[i], p_窓長).Any(x => !p_広い形群[j].Contains(x, StringComparison.Ordinal)));
                    l_最少 = Math.Min(l_最少, l_相手だけの窓群.Length);
                }

                if (l_Is満たす)
                {
                    l_満たす形群.Add((i, l_最少));
                }
            }

            p_条件を満たす数 = l_満たす形群.Count;
            p_相手の窓の最少 = l_満たす形群.Count == 1 ? l_満たす形群[0].A_最少 : 0;
            return l_満たす形群.Count == 1 ? l_満たす形群[0].A_番号 : null;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列の長さ p_窓長 の窓を、先頭から 1 塩基ずつずらして並べる
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_窓長">窓の長さ</param>
        /// <returns>窓の配列</returns>
        private static IEnumerable<string> Get_窓群(string p_配列, int p_窓長)
        {
            for (var l_始め = 0; l_始め + p_窓長 <= p_配列.Length; l_始め++)
            {
                yield return p_配列.Substring(l_始め, p_窓長);
            }
        }

        /// <summary>
        /// 窓がリードかその逆相補に出てくる箇所の数を、C_支持数の下限 まで数える (表に覚えて、同じ窓は 1 回だけ問い合わせる)。A・C・G・T 以外を含む窓は 0
        /// </summary>
        /// <param name="p_窓">窓</param>
        /// <param name="p_数の表">数えた窓の表</param>
        /// <param name="p_索引">リードの索引</param>
        /// <returns>出てくる箇所の数 (C_支持数の下限 で頭打ち)</returns>
        private static int Get_窓の出現数(string p_窓, Dictionary<string, int> p_数の表, ReadMinimizerIndex p_索引)
        {
            if (!p_数の表.TryGetValue(p_窓, out var l_数))
            {
                l_数 = p_窓.AsSpan().IndexOfAnyExcept(Consts.塩基文字) < 0 ? p_索引.Get_出現数(p_窓, C_支持数の下限) : 0;
                p_数の表[p_窓] = l_数;
            }

            return l_数;
        }

        #endregion
    }
}
