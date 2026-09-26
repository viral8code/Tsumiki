using Tsumiki.Commons;
using Tsumiki.Cores.Pipeline;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終成果物で、未確認の繋ぎ目を長さ不明のギャップにし、AGP に書き分ける処理の検証
    /// </summary>
    public class UnknownGapTests
    {
        #region 公開メソッド

        /// <summary>
        /// 未確認の繋ぎ目の印だけを N 100 個にし、推定長の N はそのまま残す
        /// </summary>
        /// <param name="p_配列">置き換える前の配列</param>
        /// <param name="p_長さ不明の番号">長さ不明になるはずの N の連続の番号</param>
        [Theory]
        [InlineData("ACnGT", new[] { 0 })]
        [InlineData("ACNGTnA", new[] { 1 })]
        [InlineData("ACnGTNNNAnC", new[] { 0, 2 })]
        [InlineData("ACNGT", new int[0])]
        public void Get_長さ不明のギャップへ置換_未確認の繋ぎ目だけを置き換える(string p_配列, int[] p_長さ不明の番号)
        {
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(p_配列, l_番号群);

            Assert.Equal(p_長さ不明の番号, l_番号群.Order());
            Assert.DoesNotContain(Consts.未確認の繋ぎ目, l_結果);
            Assert.Equal(p_配列.Length + (p_長さ不明の番号.Length * (FinalAssemblyPipeline.C_長さ不明のギャップ長 - 1)), l_結果.Length);
            Assert.Equal(p_配列.Replace("N", string.Empty).Replace("n", string.Empty), l_結果.Replace("N", string.Empty));
        }

        /// <summary>
        /// AGP では、長さ不明のギャップを U、推定長のギャップを N、配列の片を W として書く
        /// </summary>
        [Fact]
        public void Get_AGP行_ギャップの種類を書き分ける()
        {
            var l_番号群 = new HashSet<int>();
            var l_配列 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換("ACGTnACNNNGT", l_番号群);

            var l_行群 = FinalAssemblyPipeline.Get_AGP行("S1", l_配列, l_番号群).Select(x => x.Split('\t')).ToList();

            Assert.Equal(["W", "U", "W", "N", "W"], l_行群.Select(x => x[4]));
            Assert.Equal("100", l_行群[1][5]);
            Assert.Equal("3", l_行群[3][5]);
            Assert.Equal(["S1_1", "S1_2", "S1_3"], l_行群.Where(x => x[4] == "W").Select(x => x[5]));
            Assert.Equal(l_配列.Length.ToString(), l_行群[^1][2]);
        }

        #endregion
    }
}
