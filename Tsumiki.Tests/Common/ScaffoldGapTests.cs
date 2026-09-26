using Tsumiki.Commons;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// scaffold のギャップ (N の連続と、未確認の繋ぎ目の印) を見分ける処理の検証
    /// </summary>
    public class ScaffoldGapTests
    {
        #region 公開メソッド

        /// <summary>
        /// N は続く限り 1 つのギャップ、未確認の繋ぎ目の印は 1 文字で 1 つのギャップ
        /// </summary>
        /// <param name="p_配列">scaffold 配列</param>
        /// <param name="p_開始">ギャップの文字がある位置</param>
        /// <param name="p_期待">ギャップの直後の位置</param>
        [Theory]
        [InlineData("ACNNNGT", 2, 5)]
        [InlineData("ACnGT", 2, 3)]
        [InlineData("ACnNNGT", 2, 3)]
        [InlineData("ACNNnGT", 2, 4)]
        [InlineData("ACNN", 2, 4)]
        public void V_ギャップの終わり_Nの連続と印を分けて返す(string p_配列, int p_開始, int p_期待)
        {
            Assert.True(Util.Isギャップ文字(p_配列[p_開始]));
            Assert.Equal(p_期待, Util.Get_ギャップの終わり(p_配列, p_開始));
        }

        #endregion
    }
}
