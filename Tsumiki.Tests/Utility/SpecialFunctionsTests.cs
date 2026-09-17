using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 負の二項分布の当てはめに使う特殊関数を、閉じた形で値が分かる点で固定する
    /// </summary>
    public class SpecialFunctionsTests
    {
        #region 公開メソッド

        /// <summary>
        /// log ガンマ関数が既知の値を返す
        /// </summary>
        [Fact]
        public void V_対数ガンマが既知の値を返す()
        {
            Assert.Equal(0D, SpecialFunctions.Get_対数ガンマ(1D), 10);
            Assert.Equal(0D, SpecialFunctions.Get_対数ガンマ(2D), 10);
            Assert.Equal(Math.Log(2D), SpecialFunctions.Get_対数ガンマ(3D), 10);
            Assert.Equal(Math.Log(Math.Sqrt(Math.PI)), SpecialFunctions.Get_対数ガンマ(0.5D), 10);
            Assert.Equal(Math.Log(3_628_800D), SpecialFunctions.Get_対数ガンマ(11D), 8);
        }

        /// <summary>
        /// ディガンマ関数が既知の値を返す
        /// </summary>
        [Fact]
        public void V_ディガンマが既知の値を返す()
        {
            const double l_オイラーの定数 = 0.5772156649015329D;
            Assert.Equal(-l_オイラーの定数, SpecialFunctions.Get_ディガンマ(1D), 9);
            Assert.Equal(1D - l_オイラーの定数, SpecialFunctions.Get_ディガンマ(2D), 9);
            Assert.Equal(1.5D - l_オイラーの定数, SpecialFunctions.Get_ディガンマ(3D), 9);
        }

        /// <summary>
        /// トリガンマ関数が既知の値を返す
        /// </summary>
        [Fact]
        public void V_トリガンマが既知の値を返す()
        {
            var l_π二乗の6分の1 = Math.PI * Math.PI / 6D;
            Assert.Equal(l_π二乗の6分の1, SpecialFunctions.Get_トリガンマ(1D), 9);
            Assert.Equal(l_π二乗の6分の1 - 1D, SpecialFunctions.Get_トリガンマ(2D), 9);
            Assert.Equal(l_π二乗の6分の1 - 1.25D, SpecialFunctions.Get_トリガンマ(3D), 9);
        }

        /// <summary>
        /// ディガンマ・トリガンマが log ガンマの差分と整合する
        /// </summary>
        [Fact]
        public void V_ディガンマとトリガンマが対数ガンマの微分と整合する()
        {
            const double l_刻み = 1e-5D;
            foreach (var l_点 in new[] { 0.7D, 2.3D, 9.5D, 120D })
            {
                var l_数値微分 = (SpecialFunctions.Get_対数ガンマ(l_点 + l_刻み) - SpecialFunctions.Get_対数ガンマ(l_点 - l_刻み)) / (2D * l_刻み);
                Assert.Equal(l_数値微分, SpecialFunctions.Get_ディガンマ(l_点), 5);

                var l_二階 = (SpecialFunctions.Get_ディガンマ(l_点 + l_刻み) - SpecialFunctions.Get_ディガンマ(l_点 - l_刻み)) / (2D * l_刻み);
                Assert.Equal(l_二階, SpecialFunctions.Get_トリガンマ(l_点), 5);
            }
        }

        #endregion
    }
}
