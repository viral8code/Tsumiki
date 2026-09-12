using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 局所経路探索の状態識別を検証する
    /// </summary>
    public class ConstrainedPathFinderTests
    {
        #region 公開メソッド

        /// <summary>
        /// CT と GT を同じ探索状態にまとめない
        /// </summary>
        [Fact]
        public void V_状態キー_四番目の塩基を二ビットに収める()
        {
            var l_集合 = new LocalKmerSet(2);
            l_集合.V_登録(new byte[] { 1, 1 });
            l_集合.V_登録(new byte[] { 1, 2 });
            l_集合.V_登録(new byte[] { 1, 3 });
            l_集合.V_登録(new byte[] { 2, 4 });
            l_集合.V_登録(new byte[] { 3, 4 });
            var l_結果 = ConstrainedPathFinder.Get_経路([1, 1], [3, 4], 0, 0, l_集合, 2);
            Assert.Equal(ギャップ充填判定.充填済み, l_結果.A_判定);
            Assert.Equal(string.Empty, l_結果.A_経路);
        }

        #endregion
    }
}
