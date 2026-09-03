using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer のカウントが、シャード数(スレッド数)に依らず正確であることを確認する。
    ///
    /// k-mer をワーカー単位ではなくハッシュ値でシャードへ振り分けるようにした際、
    /// 実データでカウントがちょうど2倍になる不具合が出た(ヒストグラムが
    /// 偶数のカウントしか持たない、という形で表面化した)。
    /// スレッド数を変えて同じ答えになることを固定しておく。
    /// </summary>
    public class CountingShardTests : IDisposable
    {
        private readonly string _tempDir;

        public CountingShardTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_shard_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        [Theory]
        [InlineData(1)]
        [InlineData(2)]
        [InlineData(8)]
        public void GetCoverage_CountsEachAdditionExactlyOnce_RegardlessOfShardCount(int threadCount)
        {
            const int k = 21;
            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = threadCount };

            using var index = new TrustedKmerIndex(this._tempDir);

            var seq = "ACGGTCATTGACCTAGGATCA"; // 21塩基
            var kmer = seq.Select(Util.GetSimpleNucleotideID).ToArray();

            // ちょうど7回登録する(奇数にして「2倍になっていないか」を確実に見る)。
            for (var i = 0; i < 7; i++)
            {
                index.Add(kmer.AsSpan(), workerIndex: i % threadCount);
            }

            _ = index.Cutoff(bounds: 2);

            Assert.Equal(7UL, index.GetCoverage(kmer));
        }

        /// <summary>
        /// 多数の異なる k-mer を、それぞれ異なる回数だけ登録しても
        /// 正確に数えられること(シャード分割とマージの整合性)。
        /// </summary>
        [Theory]
        [InlineData(1)]
        [InlineData(8)]
        public void GetCoverage_IsExactAcrossManyDistinctKmers(int threadCount)
        {
            const int k = 21;
            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = threadCount };

            using var index = new TrustedKmerIndex(this._tempDir);

            var rng = new Random(1234);
            var sequence = string.Concat(Enumerable.Range(0, 500).Select(_ => "ACGT"[rng.Next(4)]));
            var bytes = sequence.Select(Util.GetSimpleNucleotideID).ToArray();

            // 位置 i の k-mer を (i % 5) + 2 回登録する。
            var expected = new Dictionary<int, ulong>();
            for (var i = 0; i + k <= bytes.Length; i++)
            {
                var times = (ulong)((i % 5) + 2);
                expected[i] = times;
                for (ulong t = 0; t < times; t++)
                {
                    index.Add(bytes.AsSpan(i, k), workerIndex: (int)(t % (ulong)threadCount));
                }
            }

            _ = index.Cutoff(bounds: 2);

            foreach (var (position, times) in expected)
            {
                Assert.Equal(times, index.GetCoverage(bytes.AsSpan(position, k)));
            }
        }
    }
}
