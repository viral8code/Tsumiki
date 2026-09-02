using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    public class AssemblyStatsReporterTests
    {
        [Fact]
        public void Compute_EmptyInput_ReturnsAllZeros()
        {
            var stats = AssemblyStatsReporter.Compute([]);

            Assert.Equal(0, stats.SequenceCount);
            Assert.Equal(0, stats.TotalLength);
            Assert.Equal(0, stats.N50);
            Assert.Equal(0, stats.L50);
            Assert.Equal(0, stats.GcPercent);
        }

        [Fact]
        public void Compute_SingleSequence_N50EqualsItsLength()
        {
            var stats = AssemblyStatsReporter.Compute(["ACGTACGTAC"]);

            Assert.Equal(1, stats.SequenceCount);
            Assert.Equal(10, stats.TotalLength);
            Assert.Equal(10, stats.N50);
            Assert.Equal(1, stats.L50);
            Assert.Equal(10, stats.MaxLength);
            Assert.Equal(10, stats.MinLength);
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

            var stats = AssemblyStatsReporter.Compute(sequences);

            Assert.Equal(10, stats.SequenceCount);
            Assert.Equal(550, stats.TotalLength);
            Assert.Equal(70, stats.N50);
            Assert.Equal(4, stats.L50);
            Assert.Equal(100, stats.MaxLength);
            Assert.Equal(10, stats.MinLength);
        }

        [Fact]
        public void Compute_GcPercent_IgnoresNRunsAndIsCaseInsensitive()
        {
            // G/C: 4, A/T: 4, N: 2 -> GC% は N を除いた8塩基中4塩基 = 50%。
            var stats = AssemblyStatsReporter.Compute(["ggccaattNN"]);

            Assert.Equal(50.0, stats.GcPercent, precision: 6);
        }

        [Fact]
        public void ComputeFromFasta_ReadsSequencesFromFile()
        {
            var path = Path.GetTempFileName();
            try
            {
                File.WriteAllText(path, ">seq1\nACGTACGTAC\n>seq2\nACGT\n");

                var stats = AssemblyStatsReporter.ComputeFromFasta(path);

                Assert.Equal(2, stats.SequenceCount);
                Assert.Equal(14, stats.TotalLength);
            }
            finally
            {
                File.Delete(path);
            }
        }
    }
}
