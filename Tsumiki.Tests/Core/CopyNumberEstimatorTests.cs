using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// カバレッジから unitig のコピー数を推定する処理の検証。
    ///
    /// ゲノム中に1回しか現れない領域のカバレッジを基準値とすると、n 回現れる
    /// 反復配列にはリードが n 倍集まる。したがってカバレッジ比を丸めれば
    /// コピー数になる。これが分かると、反復配列かどうかをグラフの形ではなく
    /// 量的な根拠で判定でき、経路探索では「何回まで使ってよいか」の予算になる。
    /// </summary>
    public class CopyNumberEstimatorTests : IDisposable
    {
        private readonly string _tempDir;

        public CopyNumberEstimatorTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_copynumber_tests_" + Guid.NewGuid().ToString("N"));
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

        [Fact]
        public void Estimate_SeparatesSingleCopyFromTwoCopyAndFourCopySequences()
        {
            const int k = 21;
            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = 1 };

            // 単一コピー相当を2本(長さで基準値を支配させる)、
            // 2倍・4倍のカバレッジで登録する配列を1本ずつ用意する。
            var single1 = RandomSequence(400, seed: 1);
            var single2 = RandomSequence(400, seed: 2);
            var doubled = RandomSequence(120, seed: 3);
            var quadrupled = RandomSequence(120, seed: 4);

            using var index = new TrustedKmerIndex(this._tempDir);

            void Add(string seq, int depth)
            {
                var bytes = seq.Select(Util.GetSimpleNucleotideID).ToArray();
                for (var i = 0; i + k <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < depth; rep++)
                    {
                        index.Add(bytes.AsSpan(i, k), workerIndex: 0);
                    }
                }
            }

            Add(single1, 20);
            Add(single2, 20);
            Add(doubled, 40);
            Add(quadrupled, 80);

            _ = index.Cutoff(bounds: 2);

            Dictionary<int, string> unitigs = new()
            {
                [1] = single1,
                [2] = single2,
                [3] = doubled,
                [4] = quadrupled,
            };
            var lengths = unitigs.ToDictionary(kv => kv.Key, kv => kv.Value.Length);

            var coverage = CopyNumberEstimator.ComputeCoverage(index, unitigs, k);
            var result = CopyNumberEstimator.Estimate(coverage, lengths);

            // 基準値は長さ加重中央値なので、長い単一コピー配列の水準になるはず。
            Assert.InRange(result.Baseline, 15, 25);

            Assert.Equal(1, result.CopyNumber[1]);
            Assert.Equal(1, result.CopyNumber[2]);
            Assert.Equal(2, result.CopyNumber[3]);
            Assert.Equal(4, result.CopyNumber[4]);
        }

        /// <summary>
        /// カバレッジがわずかに高いだけの配列を反復と誤判定してはいけない。
        /// 実データのカバレッジは領域ごとにかなりばらつくため、
        /// 1.5倍未満は単一コピーとして扱う。
        /// </summary>
        [Fact]
        public void Estimate_TreatsMildlyElevatedCoverageAsSingleCopy()
        {
            const int k = 21;
            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = 1 };

            var baselineSeq = RandomSequence(400, seed: 5);
            var slightlyHigher = RandomSequence(120, seed: 6);

            using var index = new TrustedKmerIndex(this._tempDir);

            void Add(string seq, int depth)
            {
                var bytes = seq.Select(Util.GetSimpleNucleotideID).ToArray();
                for (var i = 0; i + k <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < depth; rep++)
                    {
                        index.Add(bytes.AsSpan(i, k), workerIndex: 0);
                    }
                }
            }

            Add(baselineSeq, 20);
            Add(slightlyHigher, 26); // 1.3倍

            _ = index.Cutoff(bounds: 2);

            Dictionary<int, string> unitigs = new() { [1] = baselineSeq, [2] = slightlyHigher };
            var lengths = unitigs.ToDictionary(kv => kv.Key, kv => kv.Value.Length);

            var coverage = CopyNumberEstimator.ComputeCoverage(index, unitigs, k);
            var result = CopyNumberEstimator.Estimate(coverage, lengths);

            Assert.Equal(1, result.CopyNumber[2]);
        }

        /// <summary>
        /// k-mer 長より短い unitig はカバレッジを測れないが、
        /// コピー数0にして経路から締め出してはいけない(配列自体は存在する)。
        /// </summary>
        [Fact]
        public void Estimate_UnitigShorterThanKmer_GetsCopyNumberOneRatherThanZero()
        {
            const int k = 21;
            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = 1 };

            var normal = RandomSequence(300, seed: 8);
            var tooShort = RandomSequence(10, seed: 9);

            using var index = new TrustedKmerIndex(this._tempDir);
            var bytes = normal.Select(Util.GetSimpleNucleotideID).ToArray();
            for (var i = 0; i + k <= bytes.Length; i++)
            {
                for (var rep = 0; rep < 20; rep++)
                {
                    index.Add(bytes.AsSpan(i, k), workerIndex: 0);
                }
            }
            _ = index.Cutoff(bounds: 2);

            Dictionary<int, string> unitigs = new() { [1] = normal, [2] = tooShort };
            var lengths = unitigs.ToDictionary(kv => kv.Key, kv => kv.Value.Length);

            var coverage = CopyNumberEstimator.ComputeCoverage(index, unitigs, k);
            var result = CopyNumberEstimator.Estimate(coverage, lengths);

            Assert.Equal(0, coverage[2]);
            Assert.Equal(1, result.CopyNumber[2]);
        }
    }
}
