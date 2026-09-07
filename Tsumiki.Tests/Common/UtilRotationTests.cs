using Tsumiki.Common;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// 環状配列の開始位置を辞書式順序で最小の回転へ正規化する処理(Booth の
    /// アルゴリズム)の検証。環状 contig は開始位置が walk の起点という
    /// 偶然の産物でしかないため、同じ環状配列ならどの回転から出発しても
    /// 同じ正規化結果になることが下流の比較・再現性の前提になる。
    /// </summary>
    public class UtilRotationTests
    {
        [Theory]
        [InlineData("")]
        [InlineData("A")]
        public void Get_最小回転_TrivialInputs_ReturnsAsIs(string input)
        {
            Assert.Equal(input, Util.Get_最小回転(input));
        }

        [Fact]
        public void Get_最小回転_FindsTheLexicographicallySmallestRotation()
        {
            // "BAAB" の回転は BAAB, AABB, ABBA, BBAA。辞書式最小は AABB。
            Assert.Equal("AABB", Util.Get_最小回転("BAAB"));
        }

        /// <summary>
        /// 核心の性質: 同じ環状配列のどの回転から始めても、正規化結果は同一。
        /// </summary>
        [Theory]
        [InlineData(0)]
        [InlineData(1)]
        [InlineData(5)]
        [InlineData(10)]
        [InlineData(19)]
        public void Get_最小回転_IsInvariantAcrossAllRotationsOfTheSameSequence(int rotateBy)
        {
            const string original = "ACGTTGCAACGTAGGCTTAA"; // 20bp、非反復的
            var rotated = original[rotateBy..] + original[..rotateBy];

            Assert.Equal(Util.Get_最小回転(original), Util.Get_最小回転(rotated));
        }

        [Fact]
        public void Get_最小回転_HandlesRepetitiveSequences()
        {
            // 全て同じ文字なら、どの回転でも結果は同じ文字列になる。
            const string repetitive = "AAAAAA";
            Assert.Equal(repetitive, Util.Get_最小回転(repetitive));
        }

        [Fact]
        public void Get_最小回転_ResultIsAlwaysAValidRotationOfTheInput()
        {
            const string original = "TGGCAAGTCACTCTCGACCGA";
            var result = Util.Get_最小回転(original);

            Assert.Equal(original.Length, result.Length);
            Assert.Contains(result, original + original);
        }
    }
}
