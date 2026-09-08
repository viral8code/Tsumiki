using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    public class AssemblyStatsReporterTests
    {
        [Fact]
        public void Compute_EmptyInput_ReturnsAllZeros()
        {
            var stats = AssemblyStatsReporter.Get_統計([]);

            Assert.Equal(0, stats.A_配列数);
            Assert.Equal(0, stats.A_総延長);
            Assert.Equal(0, stats.A_N50);
            Assert.Equal(0, stats.A_L50);
            Assert.Equal(0, stats.A_GC率);
        }

        [Fact]
        public void Compute_SingleSequence_N50EqualsItsLength()
        {
            var stats = AssemblyStatsReporter.Get_統計(["ACGTACGTAC"]);

            Assert.Equal(1, stats.A_配列数);
            Assert.Equal(10, stats.A_総延長);
            Assert.Equal(10, stats.A_N50);
            Assert.Equal(1, stats.A_L50);
            Assert.Equal(10, stats.A_最大長);
            Assert.Equal(10, stats.A_最小長);
        }

        [Fact]
        public void Compute_KnownN50Example()
        {
            // 長さ: 100, 90, 80, 70, 60, 50, 40, 30, 20, 10 (合計550)。
            // 半分(275)に達するのは 100+90+80+70=340 の時点(4本目)なので N50=70, L50=4。
            List<string> sequences =
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

            var stats = AssemblyStatsReporter.Get_統計(sequences);

            Assert.Equal(10, stats.A_配列数);
            Assert.Equal(550, stats.A_総延長);
            Assert.Equal(70, stats.A_N50);
            Assert.Equal(4, stats.A_L50);
            Assert.Equal(100, stats.A_最大長);
            Assert.Equal(10, stats.A_最小長);
        }

        [Fact]
        public void Compute_GcPercent_IgnoresNRunsAndIsCaseInsensitive()
        {
            // G/C: 4, A/T: 4, N: 2 -> GC% は N を除いた8塩基中4塩基 = 50%。
            var stats = AssemblyStatsReporter.Get_統計(["ggccaattNN"]);

            Assert.Equal(50.0, stats.A_GC率, precision: 6);
        }

        [Fact]
        public void ComputeFromFasta_ReadsSequencesFromFile()
        {
            var path = Path.GetTempFileName();
            try
            {
                File.WriteAllText(path, ">seq1\nACGTACGTAC\n>seq2\nACGT\n");

                var stats = AssemblyStatsReporter.Get_統計_FASTA(path);

                Assert.Equal(2, stats.A_配列数);
                Assert.Equal(14, stats.A_総延長);
            }
            finally
            {
                File.Delete(path);
            }
        }
    }
}
