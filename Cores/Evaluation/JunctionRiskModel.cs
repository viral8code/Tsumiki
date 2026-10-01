using Tsumiki.Models.Evaluation;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 継ぎ目の証拠から、繋ぐ相手を誤っている確率を出す
    /// </summary>
    internal static class JunctionRiskModel
    {
        #region 定数

        /// <summary>
        /// ペアで跨げるとみなす期待の組の数
        /// </summary>
        private const double C_跨げるとみなす期待 = 1D;

        /// <summary>
        /// 組の比を信じる期待の組の数
        /// </summary>
        private const double C_比を信じる期待 = 3D;

        /// <summary>
        /// 得点の説明変数ごとの係数 (切片、log 反復長、コピー数、log 跨ぐ読み、log 跨ぐ組、log 期待の組、組の比、期待なし、外れ割合、log 外れ錨、切れ端、低い側の深さ比、ギャップ)
        /// </summary>
        private static readonly double[] C_得点の係数 = [-4.352129400596203, 0.5317252598043638, -0.3835578043493378, 0.002719622406046047, -0.42392366115898206, 0.34530838038785383, -0.7026936261060341, -0.3314233256216074, 0.507846892920463, 0.5120702419783515, -0.07885861690564737, -0.37559067759075626, 0.024224181417808632];

        /// <summary>
        /// 得点の説明変数ごとの標準化の平均
        /// </summary>
        private static readonly double[] C_得点の平均 = [0D, 4.798165102230567, 1.3770635816403263, 2.224909065530946, 3.210547747773368, 3.235847617462876, 0.8140429834462061, 0.2137885628291949, 0.20263026711813173, 1.5837253908127171, 0.03386192626034558, 1.1195739277652421, 0.0035741158765989467];

        /// <summary>
        /// 得点の説明変数ごとの標準化の幅
        /// </summary>
        private static readonly double[] C_得点の幅 = [1D, 1.2382849707315227, 0.5989939148506482, 1.9421077467545627, 1.8045768799211892, 1.744288246714602, 0.48723385814245174, 0.409979283906658, 0.308402145166564, 1.0305702532096517, 0.12464616938405142, 0.346636093137907, 0.05967697690316483];

        /// <summary>
        /// ギャップの組の較正 (切片、得点、低い側の深さ比、切れ端)
        /// </summary>
        private static readonly double[] C_ギャップの較正 = [-3.442442144518912, 1.3226755816286004, -0.6038858069649776, 0.712143698115858];

        /// <summary>
        /// 組で跨げる組の較正 (切片、得点)
        /// </summary>
        private static readonly double[] C_跨げる組の較正 = [0.6126689139472675, 2.617854620679609];

        /// <summary>
        /// 組で跨げない組の較正 (切片、得点、低い側の深さ比、切れ端)
        /// </summary>
        private static readonly double[] C_跨げない組の較正 = [0.3030979086174536, 0.928841530928149, -2.7693864888039, 2.177682353092452];

        /// <summary>
        /// コピー数の上限
        /// </summary>
        private const double C_コピー数の上限 = 5D;

        /// <summary>
        /// 組の比と低い側の深さ比の上限
        /// </summary>
        private const double C_比の上限 = 2D;

        /// <summary>
        /// 切れ端の上限
        /// </summary>
        private const double C_切れ端の上限 = 3D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 候補を組に分け、繋ぐ相手を誤っている確率を付ける
        /// </summary>
        /// <param name="p_候補"></param>
        /// <returns>評価</returns>
        public static 継ぎ目の評価 Get_評価(継ぎ目候補 p_候補)
        {
            var l_得点 = Get_得点(p_候補);
            var l_深さ = Math.Min(p_候補.A_低い側の深さ比, C_比の上限);
            var l_切れ端 = Math.Min(p_候補.A_切れ端, C_切れ端の上限);
            var l_組 = Get_組(p_候補);
            var l_線形和 = l_組 switch
            {
                継ぎ目の組.ギャップ => C_ギャップの較正[0] + (C_ギャップの較正[1] * l_得点) + (C_ギャップの較正[2] * l_深さ) + (C_ギャップの較正[3] * l_切れ端),
                継ぎ目の組.組で跨げる => C_跨げる組の較正[0] + (C_跨げる組の較正[1] * l_得点),
                _ => C_跨げない組の較正[0] + (C_跨げない組の較正[1] * l_得点) + (C_跨げない組の較正[2] * l_深さ) + (C_跨げない組の較正[3] * l_切れ端),
            };
            return new 継ぎ目の評価(p_候補, l_組, Get_シグモイド(l_線形和));
        }

        /// <summary>
        /// 候補の組
        /// </summary>
        /// <param name="p_候補"></param>
        /// <returns></returns>
        public static 継ぎ目の組 Get_組(継ぎ目候補 p_候補)
        {
            return p_候補.A_Isギャップ ? 継ぎ目の組.ギャップ : p_候補.A_期待の組 >= C_跨げるとみなす期待 ? 継ぎ目の組.組で跨げる : 継ぎ目の組.組で跨げない;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 人工の誤結合で学習した得点 (ロジット)
        /// </summary>
        /// <param name="p_候補"></param>
        /// <returns></returns>
        internal static double Get_得点(継ぎ目候補 p_候補)
        {
            var l_Is比を信じる = p_候補.A_期待の組 >= C_比を信じる期待;
            ReadOnlySpan<double> l_説明変数 =
            [
                1D,
                Math.Log(1D + p_候補.A_反復長),
                Math.Min(p_候補.A_コピー数, C_コピー数の上限),
                Math.Log(1D + p_候補.A_跨ぐ読み),
                Math.Log(1D + p_候補.A_跨ぐ組),
                Math.Log(1D + p_候補.A_期待の組),
                l_Is比を信じる ? Math.Min(p_候補.A_組の比, C_比の上限) : 0D,
                l_Is比を信じる ? 0D : 1D,
                p_候補.A_外れ割合,
                Math.Log(1D + p_候補.A_外れ錨),
                Math.Min(p_候補.A_切れ端, C_切れ端の上限),
                Math.Min(p_候補.A_低い側の深さ比, C_比の上限),
                p_候補.A_Isギャップ ? 1D : 0D,
            ];
            var l_得点 = 0D;
            for (var i = 0; i < l_説明変数.Length; i++)
            {
                l_得点 += C_得点の係数[i] * (l_説明変数[i] - C_得点の平均[i]) / C_得点の幅[i];
            }

            return l_得点;
        }

        /// <summary>
        /// ロジットを確率にする
        /// </summary>
        /// <param name="p_ロジット"></param>
        /// <returns></returns>
        private static double Get_シグモイド(double p_ロジット)
        {
            return 1D / (1D + Math.Exp(-Math.Clamp(p_ロジット, -30D, 30D)));
        }

        #endregion
    }
}
