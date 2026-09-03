using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// スキャフォールドのギャップ(N の連続)を、de Bruijn グラフ上で
    /// 両端を繋ぐ経路を探して実配列に置き換える処理の検証。
    ///
    /// contig が途切れるのは配列が存在しないからではなく、分岐でどちらへ
    /// 進むか決められなかったからであることが多い。その場合ギャップを埋める
    /// 配列は k-mer 集合の中に実在しており、両端から辿れば復元できる。
    /// </summary>
    public class GapFillerTests : IDisposable
    {
        private readonly string _tempDir;

        public GapFillerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_gapfiller_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        private TrustedKmerIndex BuildIndex(int kmerLength, params string[] sequences)
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = kmerLength, ThreadCount = 1 };
            var index = new TrustedKmerIndex(this._tempDir);
            foreach (var seq in sequences)
            {
                var bytes = seq.Select(Util.GetSimpleNucleotideID).ToArray();
                for (var i = 0; i + kmerLength <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < 3; rep++)
                    {
                        index.Add(bytes.AsSpan(i, kmerLength), workerIndex: 0);
                    }
                }
            }
            _ = index.Cutoff(bounds: 2);
            return index;
        }

        private string WriteScaffold(string name, string sequence)
        {
            var path = Path.Combine(this._tempDir, name);
            using (var writer = new FastaWriter(path))
            {
                writer.Write("SCAFFOLD1", sequence);
            }
            return path;
        }

        private static string ReadSingleSequence(string path)
        {
            using var reader = new FastaReader(path);
            Assert.True(reader.HasNext());
            return reader.NextSequence().Seq;
        }

        [Fact]
        public void Run_UniquePathThroughTheGraph_RestoresTheTrueSequence()
        {
            const int k = 21;
            // 200bp の非反復的な配列。k=21 なので偶然の重複はまず起きない。
            var truth = RandomSequence(200, seed: 20260903);

            using var index = this.BuildIndex(k, truth);

            // 真ん中 40bp を N に置き換えたスキャフォールドを作る。
            const int gapStart = 80;
            const int gapLength = 40;
            var withGap = truth[..gapStart] + new string('N', gapLength) + truth[(gapStart + gapLength)..];
            var path = this.WriteScaffold("scaffolds.fasta", withGap);

            var stats = GapFiller.Run(path, index, k);

            Assert.Equal(1, stats.TotalGaps);
            Assert.Equal(1, stats.FilledGaps);
            Assert.Equal(gapLength, stats.FilledBases);

            // 埋めた結果は元の配列そのものに戻っていなければならない。
            Assert.Equal(truth, ReadSingleSequence(path));
        }

        [Fact]
        public void Run_GapLengthEstimateSlightlyOff_StillFillsUsingTheMargin()
        {
            const int k = 21;
            var truth = RandomSequence(200, seed: 7);

            using var index = this.BuildIndex(k, truth);

            // 実際の欠損は 40bp だが、推定を誤って 30 個の N になっている状況。
            // ギャップ長推定はインサートサイズ推定のばらつきを引き継ぐため、
            // ぴったりの長さしか探さないと現実にはまず埋まらない。
            const int gapStart = 80;
            const int actualMissing = 40;
            var withGap = truth[..gapStart] + new string('N', 30) + truth[(gapStart + actualMissing)..];
            var path = this.WriteScaffold("scaffolds_off.fasta", withGap);

            var stats = GapFiller.Run(path, index, k);

            Assert.Equal(1, stats.FilledGaps);
            Assert.Equal(truth, ReadSingleSequence(path));
        }

        /// <summary>
        /// ギャップを埋める経路が複数ある場合、どれが正しいか決められない。
        /// 誤った配列で埋めるより N のまま残すほうが下流の解析にとって安全。
        /// </summary>
        [Fact]
        public void Run_MultiplePathsFitTheGap_LeavesItAsNRatherThanGuessing()
        {
            const int k = 21;
            var prefix = RandomSequence(80, seed: 11);
            var suffix = RandomSequence(80, seed: 12);
            // 同じ長さで中身だけ違う2通りの中間配列を、どちらも k-mer 集合に入れる。
            var middleA = RandomSequence(40, seed: 13);
            var middleB = RandomSequence(40, seed: 14);

            using var index = this.BuildIndex(k, prefix + middleA + suffix, prefix + middleB + suffix);

            var withGap = prefix + new string('N', 40) + suffix;
            var path = this.WriteScaffold("scaffolds_ambiguous.fasta", withGap);

            var stats = GapFiller.Run(path, index, k);

            Assert.Equal(1, stats.TotalGaps);
            Assert.Equal(0, stats.FilledGaps);
            Assert.Equal(1, stats.AmbiguousGaps);
            // N はそのまま残っていること。
            Assert.Contains('N', ReadSingleSequence(path));
        }

        /// <summary>
        /// 両端を繋ぐ経路がグラフ上に存在しない(本当に配列が無い)場合は、
        /// 当然埋められない。
        /// </summary>
        [Fact]
        public void Run_NoPathConnectsTheTwoSides_LeavesItAsN()
        {
            const int k = 21;
            var left = RandomSequence(80, seed: 21);
            var right = RandomSequence(80, seed: 22);

            // 左右それぞれの k-mer は入れるが、両者を繋ぐ配列は入れない。
            using var index = this.BuildIndex(k, left, right);

            var withGap = left + new string('N', 40) + right;
            var path = this.WriteScaffold("scaffolds_unreachable.fasta", withGap);

            var stats = GapFiller.Run(path, index, k);

            Assert.Equal(1, stats.TotalGaps);
            Assert.Equal(0, stats.FilledGaps);
            Assert.Equal(1, stats.UnreachableGaps);
            Assert.Contains('N', ReadSingleSequence(path));
        }

        [Fact]
        public void Run_NoGaps_LeavesTheSequenceUntouched()
        {
            const int k = 21;
            var truth = RandomSequence(150, seed: 31);
            using var index = this.BuildIndex(k, truth);

            var path = this.WriteScaffold("scaffolds_nogap.fasta", truth);
            var stats = GapFiller.Run(path, index, k);

            Assert.Equal(0, stats.TotalGaps);
            Assert.Equal(truth, ReadSingleSequence(path));
        }
    }
}
