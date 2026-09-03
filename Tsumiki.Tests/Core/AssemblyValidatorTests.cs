using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// アセンブリが観測された k-mer とその出現回数に対して辻褄が合っているかを
    /// 確かめる自己検査の検証。リファレンス配列なしで「取りこぼし」と
    /// 「出しすぎ」を検出できることを固定する。
    ///
    /// 「出しすぎ」の検出は特に重要で、総延長が実際のゲノムサイズより大きく
    /// なる原因はほぼこれ(実際、修正前は同じ配列を順鎖と逆鎖の両方で出力して
    /// いて総長がちょうど2.009倍に膨れていた)。
    /// </summary>
    public class AssemblyValidatorTests : IDisposable
    {
        private readonly string _tempDir;

        public AssemblyValidatorTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_validator_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private const int K = 21;

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        private TrustedKmerIndex BuildIndex(int depth, params string[] sequences)
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = K, ThreadCount = 1 };
            var index = new TrustedKmerIndex(this._tempDir);
            foreach (var seq in sequences)
            {
                var bytes = seq.Select(Util.GetSimpleNucleotideID).ToArray();
                for (var i = 0; i + K <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < depth; rep++)
                    {
                        index.Add(bytes.AsSpan(i, K), workerIndex: 0);
                    }
                }
            }
            _ = index.Cutoff(bounds: 2);
            return index;
        }

        private string WriteFasta(string name, params string[] sequences)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new FastaWriter(path);
            var id = 1;
            foreach (var seq in sequences)
            {
                writer.Write($"NODE{id++}", seq);
            }
            return path;
        }

        [Fact]
        public void Validate_AssemblyThatExactlyReproducesTheInput_ReportsNoMissingAndNoExcess()
        {
            var truth = RandomSequence(600, seed: 101);
            using var index = this.BuildIndex(depth: 20, truth);

            var path = this.WriteFasta("perfect.fasta", truth);
            var result = AssemblyValidator.Validate(path, index, K, singleCopyBaseline: 20);

            Assert.Equal(0, result.MissingKmers);
            Assert.Equal(0, result.ExcessInstances);
        }

        [Fact]
        public void Validate_AssemblyMissingHalfTheSequence_ReportsTheMissingKmers()
        {
            var truth = RandomSequence(600, seed: 102);
            using var index = this.BuildIndex(depth: 20, truth);

            // 後半を落としたアセンブリ。
            var path = this.WriteFasta("truncated.fasta", truth[..300]);
            var result = AssemblyValidator.Validate(path, index, K, singleCopyBaseline: 20);

            Assert.True(result.MissingKmers > 0, "truncated assembly should report missing k-mers");
            // 600bp の k-mer は 580 個、そのうち前半 300bp に含まれるのは 280 個。
            Assert.Equal(580 - 280, result.MissingKmers);
            Assert.InRange(result.MissingPercent, 45, 55);
        }

        /// <summary>
        /// 単一コピーの配列を2回出力してしまった場合、カバレッジは1コピー分しか
        /// 無いので「出しすぎ」として検出されなければならない。これが検出できないと、
        /// 総延長が水増しされていることに気付けない。
        /// </summary>
        [Fact]
        public void Validate_SingleCopySequenceEmittedTwice_ReportsItAsExcess()
        {
            var truth = RandomSequence(600, seed: 103);
            using var index = this.BuildIndex(depth: 20, truth);

            var path = this.WriteFasta("duplicated.fasta", truth, truth);
            var result = AssemblyValidator.Validate(path, index, K, singleCopyBaseline: 20);

            Assert.Equal(0, result.MissingKmers);
            // 各 k-mer が期待の2倍出ているので、延べ数の半分が余分。
            Assert.Equal(580, result.OverRepresentedKmers);
            Assert.Equal(580, result.ExcessInstances);
            Assert.InRange(result.ExcessPercent, 45, 55);
        }

        /// <summary>
        /// 逆相補で出力されていても同じ配列とみなされること(正規化の確認)。
        /// これが効いていないと、逆鎖側の contig がすべて「取りこぼし」に見えてしまう。
        /// </summary>
        [Fact]
        public void Validate_ReverseComplementedAssembly_IsTreatedAsTheSameSequence()
        {
            var truth = RandomSequence(600, seed: 104);
            using var index = this.BuildIndex(depth: 20, truth);

            var path = this.WriteFasta("revcomp.fasta", Util.ReverseComprement(truth));
            var result = AssemblyValidator.Validate(path, index, K, singleCopyBaseline: 20);

            Assert.Equal(0, result.MissingKmers);
            Assert.Equal(0, result.ExcessInstances);
        }

        /// <summary>
        /// 2コピー分のカバレッジがある反復配列を2回出力するのは正しい。
        /// これを「出しすぎ」と誤判定してはいけない。
        /// </summary>
        [Fact]
        public void Validate_TwoCopyRepeatEmittedTwice_IsNotCountedAsExcess()
        {
            var single = RandomSequence(600, seed: 105);
            var repeat = RandomSequence(200, seed: 106);

            ConfigurationManager.Arguments = new Parameters { Kmer = K, ThreadCount = 1 };
            using var index = new TrustedKmerIndex(this._tempDir);

            void Add(string seq, int depth)
            {
                var bytes = seq.Select(Util.GetSimpleNucleotideID).ToArray();
                for (var i = 0; i + K <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < depth; rep++)
                    {
                        index.Add(bytes.AsSpan(i, K), workerIndex: 0);
                    }
                }
            }

            Add(single, 20);
            Add(repeat, 40); // 2コピー相当のカバレッジ
            _ = index.Cutoff(bounds: 2);

            var path = this.WriteFasta("repeat_twice.fasta", single, repeat, repeat);
            var result = AssemblyValidator.Validate(path, index, K, singleCopyBaseline: 20);

            Assert.Equal(0, result.MissingKmers);
            Assert.Equal(0, result.ExcessInstances);
        }
    }
}
