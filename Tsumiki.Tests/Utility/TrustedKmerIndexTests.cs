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

        /// <summary>
        /// 長さ len の直鎖配列(分岐なし)の全k-merをカットオフ以上登録する。
        /// GraphSimplifierのテストとも共通で使える小さなヘルパー。
        /// </summary>
        private TrustedKmerIndex BuildLinearIndex(string seq, int kmerLength)
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = kmerLength, ThreadCount = 1 };
            var index = new TrustedKmerIndex(this._tempDir);
            var bytes = ToBytes(seq);
            for (var i = 0; i + kmerLength <= bytes.Length; i++)
            {
                for (var rep = 0; rep < 3; rep++)
                {
                    index.Add(bytes.AsSpan(i, kmerLength), workerIndex: 0);
                }
            }
            _ = index.Cutoff(bounds: 2);
            return index;
        }

        [Fact]
        public void CountOutEdges_IsOneInsideALinearSequence_AndZeroAtItsEnd()
        {
            const string seq = "GCTAAAGACAATTACATAACATAC"; // 24bp、非周期的(内部にk=8の重複なし)
            const int k = 8;
            using var index = this.BuildLinearIndex(seq, k);
            var bytes = ToBytes(seq);

            // 途中のk-mer: ちょうど1通りだけ後続がある。
            Assert.Equal(1, index.CountOutEdges(bytes.AsSpan(0, k)));

            // 配列の末尾k-mer: これ以上後続がない(out-degree 0)。
            Assert.Equal(0, index.CountOutEdges(bytes.AsSpan(bytes.Length - k, k)));
        }

        [Fact]
        public void CountInEdges_IsZeroAtSequenceStart()
        {
            const string seq = "GCTAAAGACAATTACATAACATAC";
            const int k = 8;
            using var index = this.BuildLinearIndex(seq, k);
            var bytes = ToBytes(seq);

            Assert.Equal(0, index.CountInEdges(bytes.AsSpan(0, k)));
        }

        [Fact]
        public void RemoveTrusted_RemovesBothOrientations()
        {
            const string seq = "GCTAAAGACAATTACATAACATAC";
            const int k = 8;
            using var index = this.BuildLinearIndex(seq, k);
            var kmer = ToBytes(seq[..k]);
            var revComp = ToBytes(Util.ReverseComprement(seq[..k]));

            Assert.True(index.Contains(kmer));

            index.RemoveTrusted(kmer);

            Assert.False(index.Contains(kmer));
            Assert.False(index.Contains(revComp));
        }

        [Fact]
        public void EnumerateTrustedKmers_YieldsExpectedCount()
        {
            const string seq = "GCTAAAGACAATTACATAACATAC"; // 24bp、非周期的(内部にk=8の重複なし)
            const int k = 8;
            using var index = this.BuildLinearIndex(seq, k);

            // 24bp・k=8 の非周期的な直鎖配列は 24-8+1=17 個のユニークk-mer位置を持ち、
            // 内部に重複(自己一致・逆相補との一致)がないよう検証済みの配列なので、
            // 正規化後もちょうど17件になるはず。
            var count = index.EnumerateTrustedKmers().Count();
            Assert.Equal(17, count);
        }

        [Fact]
        public void FindFirstKmers_FindsTheSingleStartOfALinearSequence()
        {
            // "ACGT"の繰り返しだと逆相補と自己一致してしまい分岐点が
            // 複雑になるため、非周期的な配列を使う。
            const string seq = "GCTAAAGACAATTACATAACATAC"; // 非周期的
            const int k = 8;
            using var index = this.BuildLinearIndex(seq, k);
            var bytes = ToBytes(seq);

            var firstKmers = index.FindFirstKmers();

            // 開始k-merとして、配列の先頭(またはその正規化された逆鎖)が
            // 含まれているはず。
            var startKmer = bytes.AsSpan(0, k).ToArray();
            var startRevComp = ToBytes(Util.ReverseComprement(seq[..k]));
            Assert.Contains(firstKmers, fk => fk.SequenceEqual(startKmer) || fk.SequenceEqual(startRevComp));
        }
    }
}
