using Tsumiki.Common;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// パック値を転がしながら進める walk が、従来の実装と同じ結果を返すことを固定する<br/>
    /// 転がし更新は unitig 構築の時間のほとんどを占めていた O(k) の詰め直しを
    /// 省くためのもので、結果は 1 塩基たりとも変わってはいけない<br/>
    /// 2 つの実装が並存する以上、等価性の確認は必須になる
    /// </summary>
    public class UnitigWalkTests : IDisposable
    {
        private readonly string _tempDir;

        public UnitigWalkTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_walk_tests_" + Guid.NewGuid().ToString("N"));
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
            var l_作業 = Path.Combine(this._tempDir, Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_作業);
            var index = new TrustedKmerIndex(l_作業);
            foreach (var seq in sequences)
            {
                var bytes = seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + kmerLength <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < 5; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, kmerLength), p_ワーカー番号: 0);
                    }
                }
            }
            _ = index.V_カットオフ(p_カットオフ: 2);
            return index;
        }

        /// <summary>
        /// 従来実装と転がし実装が、すべての開始点で同じ配列を返すこと
        /// </summary>
        private void V_両実装が一致する(int p_k長, params string[] p_配列)
        {
            using var l_索引 = this.BuildIndex(p_k長, p_配列);
            var l_開始kmer = l_索引.Get_開始kmer一覧();
            Assert.NotEmpty(l_開始kmer);

            var l_従来 = new UnitigMaker(l_索引);
            var l_転がし = new UnitigWalk(l_索引, p_k長);
            HashSet<UInt128> l_訪問済み = [];

            foreach (var l_開始 in l_開始kmer)
            {
                var l_期待 = l_従来.Get_ユニティグ(l_開始).A_配列;
                var l_実際 = string.Concat(
                    l_転がし.Get_塩基列(l_開始, l_訪問済み).Select(Util.Get_塩基文字));
                Assert.Equal(l_期待, l_実際);
            }
        }

        [Theory]
        [InlineData(21)]
        [InlineData(31)]
        [InlineData(33)]
        [InlineData(63)]
        public void Walk_MatchesTheOriginalImplementation_OnALinearSequence(int kmerLength)
        {
            this.V_両実装が一致する(kmerLength, RandomSequence(3_000, seed: 801));
        }

        /// <summary>
        /// k=32 と k=33 は内部表現(ulong と UInt128)の境界<br/>
        /// 転がしのマスクとシフトがここで壊れやすい
        /// </summary>
        [Theory]
        [InlineData(31)]
        [InlineData(32)]
        [InlineData(33)]
        [InlineData(64)]
        public void Walk_MatchesTheOriginalImplementation_AtTheRepresentationBoundary(int kmerLength)
        {
            this.V_両実装が一致する(kmerLength, RandomSequence(2_000, seed: 802));
        }

        /// <summary>
        /// 分岐がある場合、両実装が同じ位置で walk を止めること
        /// </summary>
        [Fact]
        public void Walk_MatchesTheOriginalImplementation_WithBranches()
        {
            var l_共通 = RandomSequence(1_000, seed: 811);
            var l_枝1 = RandomSequence(800, seed: 812);
            var l_枝2 = RandomSequence(800, seed: 813);

            this.V_両実装が一致する(31, l_共通 + l_枝1, l_共通 + l_枝2);
        }

        /// <summary>
        /// 反復配列を含む場合<br/>
        /// 合流点の入次数判定が両実装で一致すること
        /// </summary>
        [Fact]
        public void Walk_MatchesTheOriginalImplementation_WithARepeat()
        {
            var a = RandomSequence(800, seed: 821);
            var r = RandomSequence(300, seed: 822);
            var b = RandomSequence(800, seed: 823);
            var c = RandomSequence(800, seed: 824);

            this.V_両実装が一致する(31, a + r + b + r + c);
        }

        /// <summary>
        /// 環状配列<br/>
        /// 循環検出の打ち切り位置が両実装で一致すること<br/>
        /// 完全な環には開始 k-mer が存在しない(どの k-mer も入次数 1 で、
        /// その予測元の出次数も 1)ため、任意の k-mer から walk して比べる
        /// </summary>
        [Fact]
        public void Walk_MatchesTheOriginalImplementation_OnACircularSequence()
        {
            const int k = 31;
            var l_環 = RandomSequence(1_500, seed: 831);
            var l_配列 = l_環 + l_環[..(k - 1)];

            using var l_索引 = this.BuildIndex(k, l_配列);
            Assert.Empty(l_索引.Get_開始kmer一覧());

            var l_開始 = l_配列[..k].Select(Util.Get_塩基ID).ToArray();
            var l_期待 = new UnitigMaker(l_索引).Get_ユニティグ(l_開始).A_配列;
            var l_実際 = string.Concat(
                new UnitigWalk(l_索引, k).Get_塩基列(l_開始, []).Select(Util.Get_塩基文字));

            Assert.Equal(l_期待, l_実際);
        }

        /// <summary>
        /// 逆相補側から始めても一致すること(正規形の判定が転がしでも正しいこと)
        /// </summary>
        [Fact]
        public void Walk_MatchesTheOriginalImplementation_FromTheReverseStrand()
        {
            var l_配列 = RandomSequence(2_000, seed: 841);
            this.V_両実装が一致する(31, l_配列, Util.V_逆相補(l_配列));
        }

        [Fact]
        public void Walk_IsNotUsedBeyondItsRepresentationLimit()
        {
            Assert.True(UnitigWalk.Get_扱えるか(64));
            Assert.False(UnitigWalk.Get_扱えるか(65));
        }
    }
}
