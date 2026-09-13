using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Models.Foundation
{
    /// <summary>
    /// グラフ revision の契約を検証する
    /// </summary>
    public class グラフrevisionTests
    {
        #region 公開メソッド

        /// <summary>
        /// SHA-256 形式の revision を受理する
        /// </summary>
        [Fact]
        public void V_検証_SHA256を受理する()
        {
            var l_revision = new グラフrevision(0L, new string('a', 64));

            l_revision.V_検証();
        }

        /// <summary>
        /// 空のハッシュを明示的に拒否する
        /// </summary>
        [Fact]
        public void V_検証_空ハッシュを拒否する()
        {
            var l_revision = new グラフrevision(0L, string.Empty);

            _ = Assert.Throws<ArgumentException>(l_revision.V_検証);
        }

        #endregion
    }
}
