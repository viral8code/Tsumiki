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
                'A' => Consts.塩基ID.A,
                'C' => Consts.塩基ID.C,
                'G' => Consts.塩基ID.G,
                'T' => Consts.塩基ID.T,
                _ => throw new InvalidOperationException(),
            })];
        }

        [Fact]
        public void Contains_ReturnsTrueForInsertedKmer_InEitherOrientation_AndFalseForAbsentKmer()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            using var index = new TrustedKmerIndex(this._tempDir);
            var inserted = ToBytes("ACGTACGT");
            var revComp = ToBytes(Util.V_逆相補("ACGTACGT"));
            var neverInserted = ToBytes("TTTTTTTT");

            // カットオフ(2)を超えるよう複数回登録する。
            for (var i = 0; i < 5; i++)
            {
                index.V_登録(inserted.AsSpan(), p_ワーカー番号: 0);
            }

            _ = index.V_カットオフ(p_カットオフ: 2);

            Assert.True(index.Get_含まれるか(inserted));
            Assert.True(index.Get_含まれるか(revComp)); // 正規化されるため逆鎖側からの問い合わせでもヒットする
            Assert.False(index.Get_含まれるか(neverInserted));
        }

        [Fact]
        public void Cutoff_ExcludesKmersBelowThreshold()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            using var index = new TrustedKmerIndex(this._tempDir);
            var belowThreshold = ToBytes("GGGGCCCC");

            index.V_登録(belowThreshold.AsSpan(), p_ワーカー番号: 0); // 1回だけ = カットオフ2未満

            _ = index.V_カットオフ(p_カットオフ: 2);

            Assert.False(index.Get_含まれるか(belowThreshold));
        }

        [Fact]
        public void GetCoverage_SumsForwardAndReverseStrandCounts()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            using var index = new TrustedKmerIndex(this._tempDir);
            var forward = ToBytes("ACGTACGT");
            var revComp = ToBytes(Util.V_逆相補("ACGTACGT"));

            // 順鎖を3回、逆鎖を2回登録する。カウント段階では別キー扱いだが、
            // カットオフ後の正規化されたエントリでは合算されているはず。
            for (var i = 0; i < 3; i++)
            {
                index.V_登録(forward.AsSpan(), p_ワーカー番号: 0);
            }
            for (var i = 0; i < 2; i++)
            {
                index.V_登録(revComp.AsSpan(), p_ワーカー番号: 0);
            }

            _ = index.V_カットオフ(p_カットオフ: 2);

            Assert.Equal(5UL, index.Get_カバレッジ(forward));
            Assert.Equal(5UL, index.Get_カバレッジ(revComp)); // 正規化されるため同じ値
        }

        [Fact]
        public void GetCoverage_ReturnsZeroForAbsentKmer()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 8, A_スレッド数 = 1 };

            using var index = new TrustedKmerIndex(this._tempDir);
            index.V_登録(ToBytes("AAAAAAAA").AsSpan(), p_ワーカー番号: 0);
            index.V_登録(ToBytes("AAAAAAAA").AsSpan(), p_ワーカー番号: 0);

            _ = index.V_カットオフ(p_カットオフ: 2);

            Assert.Equal(0UL, index.Get_カバレッジ(ToBytes("TTTTGGGG")));
        }

        /// <summary>
        /// k が 32 を超え 64 以下のとき使われる UInt128 経路(_trustedKmersMid)と、
        /// 64 を超えたときの KmerKey フォールバック経路のそれぞれで、
        /// 正規化(順鎖・逆鎖のどちらから問い合わせても同じ結果)と
        /// カバレッジ合算が正しく行われることを確認する。
        ///
        /// 150bp リードでは k=31 のままだと 31bp 以上の反復配列がすべて潰れ
        /// contig N50 が伸びないため、k=63 前後で正しく動くことは品質上重要。
        /// </summary>
        [Theory]
        [InlineData(33)] // UInt128 経路の下限
        [InlineData(63)] // 150bp リードでの実用値
        [InlineData(64)] // UInt128 経路の上限(ちょうど128bitを使い切る)
        [InlineData(65)] // KmerKey フォールバック経路
        public void Contains_And_GetCoverage_WorkForKmerLongerThan32(int k)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k, A_スレッド数 = 1 };

            using var index = new TrustedKmerIndex(this._tempDir);
            // 逆相補と自己一致しないよう、非周期的な塩基列を決定的に生成する。
            var seq = string.Concat(Enumerable.Range(0, k).Select(i => "ACGGTCATTGAC"[(i * 7) % 12]));
            var inserted = ToBytes(seq);
            var revComp = ToBytes(Util.V_逆相補(seq));
            var neverInserted = ToBytes(new string('T', k));

            for (var i = 0; i < 3; i++)
            {
                index.V_登録(inserted.AsSpan(), p_ワーカー番号: 0);
            }
            for (var i = 0; i < 2; i++)
            {
                index.V_登録(revComp.AsSpan(), p_ワーカー番号: 0);
            }

            _ = index.V_カットオフ(p_カットオフ: 2);

            Assert.True(index.Get_含まれるか(inserted));
            Assert.True(index.Get_含まれるか(revComp));
            Assert.False(index.Get_含まれるか(neverInserted));

            // 順鎖3回 + 逆鎖2回 が同一の正規化キーへ合算されているはず。
            Assert.Equal(5UL, index.Get_カバレッジ(inserted));
            Assert.Equal(5UL, index.Get_カバレッジ(revComp));
        }

        /// <summary>
        /// k=63 の直鎖配列で、EnumerateTrustedKmers が UInt128 経路でも
        /// 正しく塩基列へ復元でき(UnpackMid)、隣接判定(CountOutEdges)が
        /// 成立することを確認する。パック/アンパックの往復が壊れていると
        /// unitig 構築が丸ごと機能しなくなるため、経路ごとに固定しておく。
        /// </summary>
        [Fact]
        public void EnumerateAndDegrees_RoundTripThroughUInt128Path()
        {
            const int k = 63;
            // 70塩基の非周期的な配列(k=63のk-merが8個取れる)。
            var seq = string.Concat(Enumerable.Range(0, 70).Select(i => "ACGGTCATTGACCTA"[(i * 11) % 15]));

            using var index = this.BuildLinearIndex(seq, k);

            var kmers = index.Get_信頼kmer一覧().ToList();
            Assert.NotEmpty(kmers);
            Assert.All(kmers, km => Assert.Equal(k, km.Length));
            // 復元した k-mer は必ず集合に含まれていなければならない。
            Assert.All(kmers, km => Assert.True(index.Get_含まれるか(km)));

            var bytes = ToBytes(seq);
            Assert.Equal(1, index.Get_出次数(bytes.AsSpan(0, k)));
            Assert.Equal(0, index.Get_出次数(bytes.AsSpan(bytes.Length - k, k)));
            Assert.Equal(0, index.Get_入次数(bytes.AsSpan(0, k)));
        }

        /// <summary>
        /// 長さ len の直鎖配列(分岐なし)の全k-merをカットオフ以上登録する。
        /// GraphSimplifierのテストとも共通で使える小さなヘルパー。
        /// </summary>
        private TrustedKmerIndex BuildLinearIndex(string seq, int kmerLength)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = kmerLength, A_スレッド数 = 1 };
            var index = new TrustedKmerIndex(this._tempDir);
            var bytes = ToBytes(seq);
            for (var i = 0; i + kmerLength <= bytes.Length; i++)
            {
                for (var rep = 0; rep < 3; rep++)
                {
                    index.V_登録(bytes.AsSpan(i, kmerLength), p_ワーカー番号: 0);
                }
            }
            _ = index.V_カットオフ(p_カットオフ: 2);
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
            Assert.Equal(1, index.Get_出次数(bytes.AsSpan(0, k)));

            // 配列の末尾k-mer: これ以上後続がない(out-degree 0)。
            Assert.Equal(0, index.Get_出次数(bytes.AsSpan(bytes.Length - k, k)));
        }

        [Fact]
        public void CountInEdges_IsZeroAtSequenceStart()
        {
            const string seq = "GCTAAAGACAATTACATAACATAC";
            const int k = 8;
            using var index = this.BuildLinearIndex(seq, k);
            var bytes = ToBytes(seq);

            Assert.Equal(0, index.Get_入次数(bytes.AsSpan(0, k)));
        }

        [Fact]
        public void RemoveTrusted_RemovesBothOrientations()
        {
            const string seq = "GCTAAAGACAATTACATAACATAC";
            const int k = 8;
            using var index = this.BuildLinearIndex(seq, k);
            var kmer = ToBytes(seq[..k]);
            var revComp = ToBytes(Util.V_逆相補(seq[..k]));

            Assert.True(index.Get_含まれるか(kmer));

            index.V_除去(kmer);

            Assert.False(index.Get_含まれるか(kmer));
            Assert.False(index.Get_含まれるか(revComp));
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
            var count = index.Get_信頼kmer一覧().Count();
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

            var firstKmers = index.Get_開始kmer一覧();

            // 開始k-merとして、配列の先頭(またはその正規化された逆鎖)が
            // 含まれているはず。
            var startKmer = bytes.AsSpan(0, k).ToArray();
            var startRevComp = ToBytes(Util.V_逆相補(seq[..k]));
            Assert.Contains(firstKmers, fk => fk.SequenceEqual(startKmer) || fk.SequenceEqual(startRevComp));
        }
    }
}
