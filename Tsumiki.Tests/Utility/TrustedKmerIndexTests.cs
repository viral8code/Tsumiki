using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    public class TrustedKmerIndexTests : IDisposable
    {
        private readonly string _tempDir;

        public TrustedKmerIndexTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_trusted_kmer_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private static byte[] ToBytes(string kmer)
        {
            return [.. kmer.Select(c => c switch
            {
                'A' => Consts.NucleotideID.A,
                'C' => Consts.NucleotideID.C,
                'G' => Consts.NucleotideID.G,
                'T' => Consts.NucleotideID.T,
                _ => throw new InvalidOperationException(),
            })];
        }

        [Fact]
        public void Contains_ReturnsTrueForInsertedKmer_InEitherOrientation_AndFalseForAbsentKmer()
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = 8, ThreadCount = 1 };

            using var index = new TrustedKmerIndex(this._tempDir);
            var inserted = ToBytes("ACGTACGT");
            var revComp = ToBytes(Util.ReverseComprement("ACGTACGT"));
            var neverInserted = ToBytes("TTTTTTTT");

            // カットオフ(2)を超えるよう複数回登録する。
            for (var i = 0; i < 5; i++)
            {
                index.Add(inserted.AsSpan(), workerIndex: 0);
            }

            _ = index.Cutoff(bounds: 2);

            Assert.True(index.Contains(inserted));
            Assert.True(index.Contains(revComp)); // 正規化されるため逆鎖側からの問い合わせでもヒットする
            Assert.False(index.Contains(neverInserted));
        }

        [Fact]
        public void Cutoff_ExcludesKmersBelowThreshold()
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = 8, ThreadCount = 1 };

            using var index = new TrustedKmerIndex(this._tempDir);
            var belowThreshold = ToBytes("GGGGCCCC");

            index.Add(belowThreshold.AsSpan(), workerIndex: 0); // 1回だけ = カットオフ2未満

            _ = index.Cutoff(bounds: 2);

            Assert.False(index.Contains(belowThreshold));
        }

        /// <summary>
        /// k=33(k&lt;=32の高速経路が使えない)でも、従来通り厳密な
        /// HashSet&lt;KmerKey&gt;経路で正しく動作することを確認する回帰テスト。
        /// </summary>
        [Fact]
        public void Contains_WorksForKmerLongerThan32_UsingKmerKeyFallbackPath()
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = 33, ThreadCount = 1 };

            using var index = new TrustedKmerIndex(this._tempDir);
            var seq = "ACGTACGTACGTACGTACGTACGTACGTACGTA"; // 34塩基(kmer長33を1つ取れる)
            var inserted = ToBytes(seq[..33]);
            var revComp = ToBytes(Util.ReverseComprement(seq[..33]));
            var neverInserted = ToBytes(new string('T', 33));

            for (var i = 0; i < 5; i++)
            {
                index.Add(inserted.AsSpan(), workerIndex: 0);
            }

            _ = index.Cutoff(bounds: 2);

            Assert.True(index.Contains(inserted));
            Assert.True(index.Contains(revComp));
            Assert.False(index.Contains(neverInserted));
        }
    }
}
