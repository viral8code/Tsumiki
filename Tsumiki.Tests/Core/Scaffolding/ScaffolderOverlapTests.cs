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

        #endregion
    }
}
