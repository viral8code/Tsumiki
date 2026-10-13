using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Models.Foundation
{
    /// <summary>
    /// 間引きの記録 (間引く前のライブラリ群、採用割合) の複製と表示の検証
    /// </summary>
    public class 間引き設定Tests
    {
        #region 公開メソッド

        /// <summary>
        /// Get_複製 で、間引く前のライブラリ群と採用割合が写り、元と独立している
        /// </summary>
        [Fact]
        public void Get_複製で間引き前ライブラリ群と採用割合が写る()
        {
            var l_元 = new Parameters
            {
                A_間引き前ライブラリ群 = [("a.fq", "b.fq")],
                A_間引きの採用割合 = 0.25
            };
            var l_複製 = l_元.Get_複製();
            Assert.Equal(l_元.A_間引き前ライブラリ群, l_複製.A_間引き前ライブラリ群);
            Assert.Equal(0.25, l_複製.A_間引きの採用割合);

            l_複製.A_間引きの採用割合 = 0.5;
            Assert.Equal(0.25, l_元.A_間引きの採用割合);
        }

        /// <summary>
        /// 既定値は、間引き前ライブラリ群が空、採用割合が 1.0
        /// </summary>
        [Fact]
        public void 既定値は間引き前ライブラリ群が空で採用割合が1_0()
        {
            var l_設定 = new Parameters();
            Assert.Empty(l_設定.A_間引き前ライブラリ群);
            Assert.Equal(1.0, l_設定.A_間引きの採用割合);
        }

        /// <summary>
        /// ToString には間引きの記録が現れない (再開の署名を変えないため)
        /// </summary>
        [Fact]
        public void ToStringには間引きの記録が現れない()
        {
            var l_設定 = new Parameters
            {
                A_間引き前ライブラリ群 = [("a.fq", "b.fq")],
                A_間引きの採用割合 = 0.3456789
            };
            var l_文字列 = l_設定.ToString();
            Assert.DoesNotContain("0.3456789", l_文字列);
            Assert.DoesNotContain("a.fq", l_文字列);
        }

        #endregion
    }
}
