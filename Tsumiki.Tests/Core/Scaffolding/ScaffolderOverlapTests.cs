using System.Text;
using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

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
        /// k-1 より短い重なりは、繋いだ配列がリードに 2 か所以上あれば畳む
        /// </summary>
        [Theory]
        [InlineData(40, 3, 40)]
        [InlineData(5, 2, 5)]
        [InlineData(0, 2, 0)]
        [InlineData(40, 1, null)]
        [InlineData(40, 0, null)]
        public void Get_リードで確かめた重なり長_リードに繋いだ配列があれば畳む(int p_真の重なり, int p_リード数, int? p_期待)
        {
            var l_乱数 = new Random(p_真の重なり * 10 + p_リード数);
            var l_共通 = Get_乱配列(l_乱数, p_真の重なり);
            var l_左 = Get_乱配列(l_乱数, 200) + l_共通;
            var l_右 = l_共通 + Get_乱配列(l_乱数, 200);
            var l_繋いだ配列 = l_左 + l_右[p_真の重なり..];
            var l_リード群 = Enumerable.Range(0, p_リード数).Select(i => l_繋いだ配列.Substring(150 + (p_真の重なり / 2) + i, 100)).Append(Get_乱配列(l_乱数, 100)).ToList();
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);

            Assert.Equal(p_期待, Scaffolder.Get_リードで確かめた重なり長(l_索引, new StringBuilder(l_左), l_右, 89, 100));
        }

        /// <summary>
        /// 繋ぎ方が 2 通りともリードにあるときは、どちらとも決めずに畳まない
        /// </summary>
        [Fact]
        public void Get_リードで確かめた重なり長_繋ぎ方が2通りあれば畳まない()
        {
            var l_乱数 = new Random(7);
            var l_単位 = Get_乱配列(l_乱数, 20);
            var l_左 = Get_乱配列(l_乱数, 200) + l_単位 + l_単位;
            var l_右 = l_単位 + l_単位 + Get_乱配列(l_乱数, 200);
            List<string> l_リード群 = [];
            foreach (var l_重なり in new[] { 40, 20 })
            {
                var l_繋いだ配列 = l_左 + l_右[l_重なり..];
                l_リード群.Add(l_繋いだ配列.Substring(170, 100));
                l_リード群.Add(l_繋いだ配列.Substring(171, 100));
            }

            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);

            Assert.Null(Scaffolder.Get_リードで確かめた重なり長(l_索引, new StringBuilder(l_左), l_右, 89, 100));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 乱数で配列を作る
        /// </summary>
        /// <param name="p_乱数">乱数</param>
        /// <param name="p_長さ">配列の長さ</param>
        /// <returns>A・C・G・T からなる配列</returns>
        private static string Get_乱配列(Random p_乱数, int p_長さ)
        {
            return new string([.. Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[p_乱数.Next(4)])]);
        }

        #endregion

    }
}
