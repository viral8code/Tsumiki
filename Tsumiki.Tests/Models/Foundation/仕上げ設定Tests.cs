using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Models.Foundation
{
    /// <summary>
    /// v0.2 仕上げ設定の契約を検証する
    /// </summary>
    public class 仕上げ設定Tests
    {
        #region 公開メソッド

        /// <summary>
        /// 既定では新しい経路を有効にしない
        /// </summary>
        [Fact]
        public void V_既定値_v02を無効にする()
        {
            var l_設定 = new 仕上げ設定();

            Assert.False(l_設定.A_Is有効);
            l_設定.V_検証();
        }

        /// <summary>
        /// 未知の schema を明示的に拒否する
        /// </summary>
        [Fact]
        public void V_検証_未知schemaを拒否する()
        {
            var l_設定 = new 仕上げ設定 { A_schemaバージョン = 3 };

            _ = Assert.Throws<NotSupportedException>(l_設定.V_検証);
        }

        /// <summary>
        /// 複製後の仕上げ設定を元設定と共有しない
        /// </summary>
        [Fact]
        public void V_複製_仕上げ設定を共有しない()
        {
            var l_元 = new Parameters { A_仕上げ設定 = new 仕上げ設定 { A_Is有効 = true } };
            var l_複製 = l_元.Get_複製();

            l_複製.A_仕上げ設定 = new 仕上げ設定();

            Assert.True(l_元.A_仕上げ設定.A_Is有効);
            Assert.False(l_複製.A_仕上げ設定.A_Is有効);
        }

        #endregion
    }
}
