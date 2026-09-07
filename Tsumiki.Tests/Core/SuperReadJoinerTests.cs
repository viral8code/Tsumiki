using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアエンドの2本を、間の未読区間ごと1本の合成リード(SuperRead)へ
    /// 統合する処理の検証。read1 の末尾 k-mer から RC(read2) の先頭 k-mer まで、
    /// 信頼できる k-mer 集合の中で経路がちょうど1本に定まったときだけ統合する。
    /// </summary>
    public class SuperReadJoinerTests : IDisposable
    {
        private readonly string _tempDir;

        public SuperReadJoinerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_superread_tests_" + Guid.NewGuid().ToString("N"));
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
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = kmerLength, A_スレッド数 = 1 };
            var index = new TrustedKmerIndex(this._tempDir);
            foreach (var seq in sequences)
            {
                var bytes = seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + kmerLength <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < 3; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, kmerLength), p_ワーカー番号: 0);
                    }
                }
            }
            _ = index.V_カットオフ(p_カットオフ: 2);
            return index;
        }

        [Fact]
        public void Get_合成配列_UniquePathBetweenTheMates_RestoresTheTrueFragment()
        {
            const int k = 21;
            var truth = RandomSequence(200, seed: 20260907);

            using var index = this.BuildIndex(k, truth);

            var read1 = truth[..80];
            var read2 = Util.V_逆相補(truth[120..200]); // RC(read2) == truth[120..200]

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Equal(truth, result);
        }

        /// <summary>
        /// 橋渡しする経路が複数ある場合、どれが正しいか決められない。
        /// 誤った配列で繋ぐより、統合を諦めて元のペアのまま残すほうが安全。
        /// </summary>
        [Fact]
        public void Get_合成配列_MultiplePathsFitTheBridge_ReturnsNullRatherThanGuessing()
        {
            const int k = 21;
            var prefix = RandomSequence(80, seed: 11);
            var suffix = RandomSequence(80, seed: 12);
            var middleA = RandomSequence(40, seed: 13);
            var middleB = RandomSequence(40, seed: 14);

            using var index = this.BuildIndex(k, prefix + middleA + suffix, prefix + middleB + suffix);

            var read1 = prefix;
            var read2 = Util.V_逆相補(suffix);

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Null(result);
        }

        [Fact]
        public void Get_合成配列_NoPathConnectsTheMates_ReturnsNull()
        {
            const int k = 21;
            var left = RandomSequence(80, seed: 21);
            var right = RandomSequence(80, seed: 22);

            // 左右それぞれの k-mer は入れるが、両者を繋ぐ配列は入れない。
            using var index = this.BuildIndex(k, left, right);

            var read1 = left;
            var read2 = Util.V_逆相補(right);

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Null(result);
        }

        [Fact]
        public void Get_合成配列_ReadShorterThanK_ReturnsNull()
        {
            const int k = 21;
            var truth = RandomSequence(100, seed: 30);
            using var index = this.BuildIndex(k, truth);

            var result = SuperReadJoiner.Get_合成配列(truth[..10], Util.V_逆相補(truth[50..]), index, k);

            Assert.Null(result);
        }

        /// <summary>
        /// read1 と RC(read2) がそのまま隣接する(橋渡しの長さが 0 の)ケースも、
        /// 特別扱いなく正しく1本に統合できること。
        /// </summary>
        [Fact]
        public void Get_合成配列_MatesMeetWithNoGap_JoinsWithoutInsertingExtraSequence()
        {
            const int k = 21;
            var truth = RandomSequence(160, seed: 40);
            using var index = this.BuildIndex(k, truth);

            var read1 = truth[..80];
            var read2 = Util.V_逆相補(truth[80..]);

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Equal(truth, result);
        }
    }
}
