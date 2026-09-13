using Tsumiki.Models.Reporting;

namespace Tsumiki.Tests.Models.Reporting
{
    /// <summary>
    /// v0.2 仕上げ実行状態の契約を検証する
    /// </summary>
    public class 仕上げ実行状態Tests
    {
        #region 公開メソッド

        /// <summary>
        /// draft 状態には理由を要求する
        /// </summary>
        [Fact]
        public void V_検証_draftに理由を要求する()
        {
            var l_結果 = new 仕上げ実行結果(仕上げ実行状態.証拠不足, null);

            _ = Assert.Throws<ArgumentException>(l_結果.V_検証);
        }

        /// <summary>
        /// 完了状態では理由を持たない
        /// </summary>
        [Fact]
        public void V_検証_完了を受理する()
        {
            var l_結果 = new 仕上げ実行結果(仕上げ実行状態.完了, null);

            l_結果.V_検証();
        }

        #endregion
    }
}
