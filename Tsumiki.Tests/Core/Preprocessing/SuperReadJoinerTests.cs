using Tsumiki.Common;
using Tsumiki.Core.Preprocessing;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアエンドの 2 本を、間の未読区間ごと 1 本の合成リード(SuperRead)へ
    /// 統合する処理の検証<br/>
    /// read1 の末尾 k-mer から RC(read2) の先頭 k-mer まで、
    /// 信頼できる k-mer 集合の中で経路がちょうど 1 本に定まったときだけ統合する
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
        /// 橋渡しする経路が複数ある場合、どれが正しいか決められない<br/>
        /// 誤った配列で繋ぐより、統合を諦めて元のペアのまま残すほうが安全
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

            // 左右それぞれの k-mer は入れるが、両者を繋ぐ配列は入れない
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
        /// 特別扱いなく正しく 1 本に統合できること
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
    
        /// <summary>
        /// 断片がリード長の 2 倍を下回るライブラリでは read1 と RC(read2) が重なる<br/>
        /// 橋渡しに必要な長さが負になり経路探索では解けないが、重なりそのものが
        /// 断片を決めるので統合できる
        /// </summary>
        [Fact]
        public void Get_合成配列_OverlappingMates_RestoresTheTrueFragment()
        {
            const int k = 21;
            var truth = RandomSequence(225, seed: 20260909);

            using var index = this.BuildIndex(k, truth);

            var read1 = truth[..150];
            var read2 = Util.V_逆相補(truth[75..225]); // RC(read2) == truth[75..225]

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Equal(truth, result);
        }

        /// <summary>
        /// 重なりに許容範囲内の不一致が残っていると、繋いだ配列の継ぎ目には
        /// どのリードにも無い k-mer が生まれる<br/>
        /// そこで弾く
        /// </summary>
        [Fact]
        public void Get_合成配列_OverlapWithAnUnseenSeam_IsRejected()
        {
            const int k = 21;
            var truth = RandomSequence(225, seed: 20260910);

            using var index = this.BuildIndex(k, truth);

            var l_文字 = truth[..150].ToCharArray();
            l_文字[145] = l_文字[145] == 'A' ? 'C' : 'A';
            var read1 = new string(l_文字);
            var read2 = Util.V_逆相補(truth[75..225]);

            Assert.Null(SuperReadJoiner.Get_合成配列(read1, read2, index, k));
        }

        /// <summary>
        /// 断片が read1 に収まっている(アダプタ読み抜け)場合は、繋いでも
        /// 長さが伸びないので重なりでは統合しない
        /// </summary>
        [Fact]
        public void Get_合成配列_FragmentShorterThanTheRead_IsNotJoinedByOverlap()
        {
            const int k = 21;
            var truth = RandomSequence(120, seed: 20260911);

            using var index = this.BuildIndex(k, truth);

            var read1 = truth;
            var read2 = Util.V_逆相補(truth);

            Assert.Null(SuperReadJoiner.Get_合成配列(read1, read2, index, k));
        }

        /// <summary>
        /// 重なりが最小長に満たないときは、偶然の一致と区別できないので繋がない
        /// </summary>
        [Fact]
        public void Get_合成配列_OverlapShorterThanTheMinimum_IsNotJoinedByOverlap()
        {
            const int k = 121;
            var truth = RandomSequence(280, seed: 20260912);

            // 重なりは 20bp (Consts.ペア結合の最小重なり長 = 40 未満)
            // k を read1 より長くして、経路探索の側も走らないようにする
            var read1 = truth[..150];
            var read2 = Util.V_逆相補(truth[130..280]);

            using var index = this.BuildIndex(k, truth);

            Assert.Null(SuperReadJoiner.Get_合成配列(read1, read2, index, k));
        }

        /// <summary>
        /// 反復配列の中では、周期のぶんだけずれた位置も同じくらい良く合う<br/>
        /// どれか一つに決められないので、重なりでは繋がない
        /// </summary>
        [Fact]
        public void Get_合成配列_OverlapThatFitsAtSeveralOffsets_IsNotJoined()
        {
            const int k = 21;
            var l_単位 = RandomSequence(30, seed: 20260918);
            var l_左 = RandomSequence(40, seed: 20260919);
            var l_右 = RandomSequence(30, seed: 20260920);
            var truth = l_左 + string.Concat(Enumerable.Repeat(l_単位, 5)) + l_右;

            using var index = this.BuildIndex(k, truth);

            var read1 = truth[..150];
            var read2 = Util.V_逆相補(truth[70..220]);

            var (_, l_重なりで結合したか, l_曖昧で捨てた数) =
                SuperReadJoiner.Get_合成配列_内訳つき(read1, read2, index, k, null);

            Assert.False(l_重なりで結合したか);
            Assert.Equal(1, l_曖昧で捨てた数);
        }

        /// <summary>
        /// 重なりが最小長に満たないときも、偶然の一致と区別できないので繋がない
        /// </summary>
        [Fact]
        public void Get_合成配列_OverlapBelowTheRaisedMinimum_IsNotJoined()
        {
            const int k = 121;
            var truth = RandomSequence(260, seed: 20260921);

            // 重なりは 40bp
            // Consts.ペア結合の最小重なり長 (60) に届かない
            var read1 = truth[..150];
            var read2 = Util.V_逆相補(truth[110..260]);

            using var index = this.BuildIndex(k, truth);

            Assert.Null(SuperReadJoiner.Get_合成配列(read1, read2, index, k));
        }
}
}
