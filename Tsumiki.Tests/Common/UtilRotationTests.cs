using Tsumiki.Commons;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// 環状配列の開始位置を辞書式順序で最小の回転へ正規化する処理 (Booth のアルゴリズム) の検証
    /// </summary>
    public class UtilRotationTests
    {
        #region 定数

        /// <summary>
        /// 項目 AABB
        /// </summary>
        private const string C_項目_AABB = "AABB";

        /// <summary>
        /// 項目 BAAB
        /// </summary>
        private const string C_項目_BAAB = "BAAB";

        /// <summary>
        /// 塩基配列 ACGTTGCAACGTAGGCTTAA
        /// </summary>
        private const string C_塩基配列_ACGTTGCAACGTAGGCTTAA = "ACGTTGCAACGTAGGCTTAA";

        /// <summary>
        /// 塩基配列 AAAAAA
        /// </summary>
        private const string C_塩基配列_AAAAAA = "AAAAAA";

        /// <summary>
        /// 塩基配列 TGGCAAGTCACTCTCGACCGA
        /// </summary>
        private const string C_塩基配列_TGGCAAGTCACTCTCGACCGA = "TGGCAAGTCACTCTCGACCGA";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 空文字と 1 文字はそのまま返す
        /// </summary>
        /// <param name="p_文字列"></param>
        [Theory]
        [InlineData("")]
        [InlineData("A")]
        public void V_自明な入力はそのまま返す(string p_文字列)
        {
            Assert.Equal(p_文字列, Util.Get_最小回転(p_文字列));
        }

        /// <summary>
        /// 辞書式最小の回転を求める
        /// </summary>
        [Fact]
        public void V_辞書式最小の回転を求める()
        {
            Assert.Equal(C_項目_AABB, Util.Get_最小回転(C_項目_BAAB));
        }

        /// <summary>
        /// 核心の性質: 同じ環状配列のどの回転から始めても、正規化結果は同一
        /// </summary>
        /// <param name="p_回転量"></param>
        [Theory]
        [InlineData(0)]
        [InlineData(1)]
        [InlineData(5)]
        [InlineData(10)]
        [InlineData(19)]
        public void V_どの回転から始めても結果は同じ(int p_回転量)
        {
            const string l_元 = C_塩基配列_ACGTTGCAACGTAGGCTTAA;
            var l_回転後 = l_元[p_回転量..] + l_元[..p_回転量];
            Assert.Equal(Util.Get_最小回転(l_元), Util.Get_最小回転(l_回転後));
        }

        /// <summary>
        /// 同じ文字の繰り返しでも結果を返す
        /// </summary>
        [Fact]
        public void V_反復配列も扱える()
        {
            const string l_反復配列 = C_塩基配列_AAAAAA;
            Assert.Equal(l_反復配列, Util.Get_最小回転(l_反復配列));
        }

        /// <summary>
        /// 結果は入力の回転のいずれかになっている
        /// </summary>
        [Fact]
        public void V_結果は入力の有効な回転になっている()
        {
            const string l_元 = C_塩基配列_TGGCAAGTCACTCTCGACCGA;
            var l_結果 = Util.Get_最小回転(l_元);
            Assert.Equal(l_元.Length, l_結果.Length);
            Assert.Contains(l_結果, l_元 + l_元);
        }

        #endregion
    }
}
