using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアエンドの 2 本を、間の未読区間ごと 1 本の合成リード (SuperRead) へ
    /// 統合する処理の検証
    /// </summary>
    /// <remarks>
    /// read1 の末尾 k-mer から RC(read2) の先頭 k-mer まで、
    /// 信頼できる k-mer 集合の中で経路がちょうど 1 本に定まったときだけ統合する
    /// </remarks>
    public class SuperReadJoinerTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public SuperReadJoinerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_superread_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_ランダム配列(int p_長さ, int p_シード)
        {
            var l_乱数 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }

        /// <summary>
        /// 与えた配列から信頼できる k-mer 集合を組み立てる
        /// </summary>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex V_構築_インデックス(int p_k長, params string[] p_配列群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 1 };
            var l_インデックス = new TrustedKmerIndex(this._tempDir);
            foreach (var l_配列 in p_配列群)
            {
                var l_塩基列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
                {
                    for (var l_回 = 0; l_回 < 3; l_回++)
                    {
                        l_インデックス.V_登録(l_塩基列.AsSpan(i, p_k長), p_ワーカー番号: 0);
                    }
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2);
            return l_インデックス;
        }

        /// <summary>
        /// ペア間の経路が一意に定まるとき、真のフラグメントを復元することを検証する
        /// </summary>
        [Fact]
        public void Get_合成配列_ペア間の経路が一意なら真のフラグメントを復元する()
        {
            const int k = 21;
            var truth = V_生成_ランダム配列(200, p_シード: 20260907);

            using var index = this.V_構築_インデックス(k, truth);

            var read1 = truth[..80];
            var read2 = Util.V_逆相補(truth[120..200]); // RC(read2) == truth[120..200]

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Equal(truth, result);
        }

        /// <summary>
        /// 橋渡しする経路が複数ある場合、どれが正しいか決められない
        /// </summary>
        /// <remarks>
        /// 誤った配列で繋ぐより、統合を諦めて元のペアのまま残すほうが安全
        /// </remarks>
        [Fact]
        public void Get_合成配列_橋渡しの経路が複数あるときは推測せずnullを返す()
        {
            const int k = 21;
            var prefix = V_生成_ランダム配列(80, p_シード: 11);
            var suffix = V_生成_ランダム配列(80, p_シード: 12);
            var middleA = V_生成_ランダム配列(40, p_シード: 13);
            var middleB = V_生成_ランダム配列(40, p_シード: 14);

            using var index = this.V_構築_インデックス(k, prefix + middleA + suffix, prefix + middleB + suffix);

            var read1 = prefix;
            var read2 = Util.V_逆相補(suffix);

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Null(result);
        }

        /// <summary>
        /// ペアを繋ぐ経路が無ければ null を返すことを検証する
        /// </summary>
        [Fact]
        public void Get_合成配列_ペアを繋ぐ経路が無ければnullを返す()
        {
            const int k = 21;
            var left = V_生成_ランダム配列(80, p_シード: 21);
            var right = V_生成_ランダム配列(80, p_シード: 22);

            // 左右それぞれの k-mer は入れるが、両者を繋ぐ配列は入れない
            using var index = this.V_構築_インデックス(k, left, right);

            var read1 = left;
            var read2 = Util.V_逆相補(right);

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Null(result);
        }

        /// <summary>
        /// リードが k より短いときは null を返すことを検証する
        /// </summary>
        [Fact]
        public void Get_合成配列_リードがk未満のときはnullを返す()
        {
            const int k = 21;
            var truth = V_生成_ランダム配列(100, p_シード: 30);
            using var index = this.V_構築_インデックス(k, truth);

            var result = SuperReadJoiner.Get_合成配列(truth[..10], Util.V_逆相補(truth[50..]), index, k);

            Assert.Null(result);
        }

        /// <summary>
        /// read1 と RC(read2) がそのまま隣接する (橋渡しの長さが 0 の) ケースも、
        /// 特別扱いなく正しく 1 本に統合できること
        /// </summary>
        [Fact]
        public void Get_合成配列_ペアが隙間なく接するときも正しく統合する()
        {
            const int k = 21;
            var truth = V_生成_ランダム配列(160, p_シード: 40);
            using var index = this.V_構築_インデックス(k, truth);

            var read1 = truth[..80];
            var read2 = Util.V_逆相補(truth[80..]);

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Equal(truth, result);
        }
    
        /// <summary>
        /// 断片がリード長の 2 倍を下回るライブラリでは read1 と RC(read2) が重なる
        /// </summary>
        /// <remarks>
        /// 橋渡しに必要な長さが負になり経路探索では解けないが、重なりそのものが
        /// 断片を決めるので統合できる
        /// </remarks>
        [Fact]
        public void Get_合成配列_ペアが重なるときは重なりから真のフラグメントを復元する()
        {
            const int k = 21;
            var truth = V_生成_ランダム配列(225, p_シード: 20260909);

            using var index = this.V_構築_インデックス(k, truth);

            var read1 = truth[..150];
            var read2 = Util.V_逆相補(truth[75..225]); // RC(read2) == truth[75..225]

            var result = SuperReadJoiner.Get_合成配列(read1, read2, index, k);

            Assert.Equal(truth, result);
        }

        /// <summary>
        /// 重なりに許容範囲内の不一致が残っていると、繋いだ配列の継ぎ目には
        /// どのリードにも無い k-mer が生まれる
        /// </summary>
        /// <remarks>
        /// そこで弾く
        /// </remarks>
        [Fact]
        public void Get_合成配列_継ぎ目に未観測のkmerができる重なりは棄却する()
        {
            const int k = 21;
            var truth = V_生成_ランダム配列(225, p_シード: 20260910);

            using var index = this.V_構築_インデックス(k, truth);

            var l_文字 = truth[..150].ToCharArray();
            l_文字[145] = l_文字[145] == 'A' ? 'C' : 'A';
            var read1 = new string(l_文字);
            var read2 = Util.V_逆相補(truth[75..225]);

            Assert.Null(SuperReadJoiner.Get_合成配列(read1, read2, index, k));
        }

        /// <summary>
        /// 断片が read1 に収まっている (アダプタ読み抜け) 場合は、繋いでも
        /// 長さが伸びないので重なりでは統合しない
        /// </summary>
        [Fact]
        public void Get_合成配列_断片がリードより短いときは重なりで統合しない()
        {
            const int k = 21;
            var truth = V_生成_ランダム配列(120, p_シード: 20260911);

            using var index = this.V_構築_インデックス(k, truth);

            var read1 = truth;
            var read2 = Util.V_逆相補(truth);

            Assert.Null(SuperReadJoiner.Get_合成配列(read1, read2, index, k));
        }

        /// <summary>
        /// 重なりが最小長に満たないときは、偶然の一致と区別できないので繋がない
        /// </summary>
        [Fact]
        public void Get_合成配列_重なりが最小長未満のときは統合しない()
        {
            const int k = 121;
            var truth = V_生成_ランダム配列(280, p_シード: 20260912);

            // 重なりは 20 bp (Consts.ペア結合の最小重なり長 = 40 未満)
            // k を read1 より長くして、経路探索の側も走らないようにする
            var read1 = truth[..150];
            var read2 = Util.V_逆相補(truth[130..280]);

            using var index = this.V_構築_インデックス(k, truth);

            Assert.Null(SuperReadJoiner.Get_合成配列(read1, read2, index, k));
        }

        /// <summary>
        /// 反復配列の中では、周期のぶんだけずれた位置も同じくらい良く合う
        /// </summary>
        /// <remarks>
        /// どれか一つに決められないので、重なりでは繋がない
        /// </remarks>
        [Fact]
        public void Get_合成配列_複数のオフセットで同程度に合う重なりは統合しない()
        {
            const int k = 21;
            var l_単位 = V_生成_ランダム配列(30, p_シード: 20260918);
            var l_左 = V_生成_ランダム配列(40, p_シード: 20260919);
            var l_右 = V_生成_ランダム配列(30, p_シード: 20260920);
            var truth = l_左 + string.Concat(Enumerable.Repeat(l_単位, 5)) + l_右;

            using var index = this.V_構築_インデックス(k, truth);

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
        public void Get_合成配列_引き上げた最小重なり長に届かないときは統合しない()
        {
            const int k = 121;
            var truth = V_生成_ランダム配列(260, p_シード: 20260921);

            // 重なりは 40 bp
            // Consts.ペア結合の最小重なり長 (60) に届かない
            var read1 = truth[..150];
            var read2 = Util.V_逆相補(truth[110..260]);

            using var index = this.V_構築_インデックス(k, truth);

            Assert.Null(SuperReadJoiner.Get_合成配列(read1, read2, index, k));
        }
}
}
