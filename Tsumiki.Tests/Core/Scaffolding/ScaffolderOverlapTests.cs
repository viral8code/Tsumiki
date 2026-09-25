using System.Text;
using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// scaffold の連結で k-1 の重なりを畳む条件
    /// </summary>
    /// <remarks>
    /// de Bruijn グラフ上で隣り合う contig は k-1 だけ重なるので、畳まずに N で繋ぐとその k-1 塩基が二重に出る
    /// </remarks>
    public class ScaffolderOverlapTests
    {
        #region 公開メソッド

        /// <summary>
        /// k-1 の重なりが一致すれば畳む
        /// </summary>
        [Fact]
        public void Get_畳める重なり長_kマイナス1の重なりが一致すれば畳む()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };
            var l_出力 = new StringBuilder("ACGGATCTGACCTTAGG");

            Assert.Equal(7, Scaffolder.Get_畳める重なり長(l_出力, "CCTTAGGTTACGACT"));
        }

        /// <summary>
        /// 重なりが一致しなければ畳まない
        /// </summary>
        [Fact]
        public void Get_畳める重なり長_重なりが一致しなければ畳まない()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };
            var l_出力 = new StringBuilder("ACGGATCTGACCTTAGG");

            Assert.Equal(0, Scaffolder.Get_畳める重なり長(l_出力, "CCTTAGATTACGACT"));
        }

        /// <summary>
        /// k-1 より短い配列は畳まない
        /// </summary>
        [Fact]
        public void Get_畳める重なり長_kマイナス1より短い配列は畳まない()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            Assert.Equal(0, Scaffolder.Get_畳める重なり長(new StringBuilder("CCTTAGG"), "CCTTAG"));
        }

        /// <summary>
        /// k-1 が一致しなくても、それより短く 15 塩基以上で完全に一致する重なりは、最長のものを畳む
        /// </summary>
        [Theory]
        [InlineData(88, 88)]
        [InlineData(40, 40)]
        [InlineData(15, 15)]
        [InlineData(14, 0)]
        [InlineData(6, 0)]
        public void Get_畳める重なり長_短い完全一致の重なりも畳む(int p_真の重なり, int p_期待)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 89, A_スレッド数 = 1 };
            var l_乱数 = new Random(p_真の重なり);
            string Get_乱配列(int p_長さ) => new([.. Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)])]);
            var l_共通 = Get_乱配列(p_真の重なり);
            var l_左 = Get_乱配列(200) + l_共通;
            var l_右 = l_共通 + Get_乱配列(200);

            Assert.Equal(p_期待, Scaffolder.Get_畳める重なり長(new StringBuilder(l_左), l_右));
        }

        #endregion
    }
}
