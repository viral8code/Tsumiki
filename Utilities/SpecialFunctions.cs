namespace Tsumiki.Utilities
{
    /// <summary>
    /// 負の二項分布の当てはめに要る、.NET に無いガンマ関数まわりの特殊関数
    /// </summary>
    internal static class SpecialFunctions
    {
        #region 定数

        /// <summary>
        /// 漸近展開へ移る前に漸化式で引き上げる引数の下限
        /// </summary>
        private const double C_漸近展開の下限 = 6D;

        /// <summary>
        /// Lanczos 近似の係数 (g = 7、n = 9)
        /// </summary>
        private static readonly double[] C_Lanczos係数 =
        [
            0.99999999999980993D,
            676.5203681218851D,
            -1259.1392167224028D,
            771.32342877765313D,
            -176.61502916214059D,
            12.507343278686905D,
            -0.13857109526572012D,
            9.9843695780195716e-6D,
            1.5056327351493116e-7D,
        ];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// log ガンマ関数 (引数は正)
        /// </summary>
        /// <param name="p_引数"></param>
        /// <returns></returns>
        public static double Get_対数ガンマ(double p_引数)
        {
            var l_x = p_引数 - 1D;
            var l_和 = C_Lanczos係数[0];
            for (var i = 1; i < C_Lanczos係数.Length; i++)
            {
                l_和 += C_Lanczos係数[i] / (l_x + i);
            }

            var l_t = l_x + 7.5D;
            return (0.5D * Math.Log(2D * Math.PI)) + ((l_x + 0.5D) * Math.Log(l_t)) - l_t + Math.Log(l_和);
        }

        /// <summary>
        /// ディガンマ関数 (log ガンマ関数の 1 階微分、引数は正)
        /// </summary>
        /// <param name="p_引数"></param>
        /// <returns></returns>
        public static double Get_ディガンマ(double p_引数)
        {
            var l_x = p_引数;
            var l_補正 = 0D;
            while (l_x < C_漸近展開の下限)
            {
                l_補正 -= 1D / l_x;
                l_x += 1D;
            }

            var l_逆2乗 = 1D / (l_x * l_x);
            return l_補正 + Math.Log(l_x) - (0.5D / l_x)
                + (l_逆2乗 * (-(1D / 12D) + (l_逆2乗 * ((1D / 120D) + (l_逆2乗 * (-(1D / 252D) + (l_逆2乗 / 240D)))))));
        }

        /// <summary>
        /// トリガンマ関数 (log ガンマ関数の 2 階微分、引数は正)
        /// </summary>
        /// <param name="p_引数"></param>
        /// <returns></returns>
        public static double Get_トリガンマ(double p_引数)
        {
            var l_x = p_引数;
            var l_補正 = 0D;
            while (l_x < C_漸近展開の下限)
            {
                l_補正 += 1D / (l_x * l_x);
                l_x += 1D;
            }

            var l_逆2乗 = 1D / (l_x * l_x);
            return l_補正 + (1D / l_x) + (0.5D * l_逆2乗)
                + (l_逆2乗 / l_x * ((1D / 6D) - (l_逆2乗 * ((1D / 30D) - (l_逆2乗 * ((1D / 42D) - (l_逆2乗 / 30D)))))));
        }

        #endregion
    }
}
