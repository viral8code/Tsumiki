using Tsumiki.Core.Evaluation;
using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// アセンブリの長さの統計を求める処理の検証
    /// </summary>
    public class AssemblyStatsReporterTests
    {
        /// <summary>
        /// 入力が空ならすべて 0 を返すことを確かめる
        /// </summary>
        [Fact]
        public void 入力が空ならすべて0を返す()
        {
            var l_stats = AssemblyStatsReporter.Get_統計([]);

            Assert.Equal(0, l_stats.A_配列数);
            Assert.Equal(0, l_stats.A_総延長);
            Assert.Equal(0, l_stats.A_N50);
            Assert.Equal(0, l_stats.A_L50);
            Assert.Equal(0, l_stats.A_GC率);
        }

        /// <summary>
        /// 配列が 1 本なら N50 はその長さと等しいことを確かめる
        /// </summary>
        [Fact]
        public void 配列が1本ならN50はその長さと等しい()
        {
            var l_stats = AssemblyStatsReporter.Get_統計(["ACGTACGTAC"]);

            Assert.Equal(1, l_stats.A_配列数);
            Assert.Equal(10, l_stats.A_総延長);
            Assert.Equal(10, l_stats.A_N50);
            Assert.Equal(1, l_stats.A_L50);
            Assert.Equal(10, l_stats.A_最大長);
            Assert.Equal(10, l_stats.A_最小長);
        }

        /// <summary>
        /// 既知の例で N50 と L50 を正しく求めることを確かめる
        /// </summary>
        [Fact]
        public void 既知の例でN50とL50を求める()
        {
            // 長さ: 100, 90, 80, 70, 60, 50, 40, 30, 20, 10 (合計 550)
            // 半分 (275) に達するのは 100+90+80+70=340 の時点 (4 本目) なので N50=70, L50=4
            List<string> l_sequences =
            [
                new string('A', 100),
                new string('A', 90),
                new string('A', 80),
                new string('A', 70),
                new string('A', 60),
                new string('A', 50),
                new string('A', 40),
                new string('A', 30),
                new string('A', 20),
                new string('A', 10),
            ];

            var l_stats = AssemblyStatsReporter.Get_統計(l_sequences);

            Assert.Equal(10, l_stats.A_配列数);
            Assert.Equal(550, l_stats.A_総延長);
            Assert.Equal(70, l_stats.A_N50);
            Assert.Equal(4, l_stats.A_L50);
            Assert.Equal(100, l_stats.A_最大長);
            Assert.Equal(10, l_stats.A_最小長);
        }

        /// <summary>
        /// GC 率は N の連続を無視し、大文字小文字を区別しないことを確かめる
        /// </summary>
        [Fact]
        public void GC率はNの連続を無視し大文字小文字を区別しない()
        {
            // G/C: 4, A/T: 4, N: 2 -> GC% は N を除いた 8 塩基中 4 塩基 = 50%
            var l_stats = AssemblyStatsReporter.Get_統計(["ggccaattNN"]);

            Assert.Equal(50.0, l_stats.A_GC率, precision: 6);
        }

        /// <summary>
        /// FASTA ファイルから配列を読み込んで統計を求めることを確かめる
        /// </summary>
        [Fact]
        public void FASTAファイルから配列を読み込んで統計を求める()
        {
            var l_path = Path.GetTempFileName();
            try
            {
                File.WriteAllText(l_path, ">seq1\nACGTACGTAC\n>seq2\nACGT\n");

                var l_stats = AssemblyStatsReporter.Get_統計_FASTA(l_path);

                Assert.Equal(2, l_stats.A_配列数);
                Assert.Equal(14, l_stats.A_総延長);
            }
            finally
            {
                File.Delete(l_path);
            }
        }
    }
}
