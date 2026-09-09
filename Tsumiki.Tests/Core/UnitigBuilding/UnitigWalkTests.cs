using Tsumiki.Common;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// パック値を転がしながら進める walk が、従来の実装と同じ結果を返すことを固定する
    /// </summary>
    /// <remarks>
    /// 転がし更新は unitig 構築の時間のほとんどを占めていた O(k) の詰め直しを
    /// 省くためのもので、結果は 1 塩基たりとも変わってはいけない<br/>
    /// 2 つの実装が並存する以上、等価性の確認は必須になる
    /// </remarks>
    public class UnitigWalkTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public UnitigWalkTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_walk_tests_" + Guid.NewGuid().ToString("N"));
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
        /// <param name="p_length">作る長さ</param>
        /// <param name="p_seed">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_length, int p_seed)
        {
            var l_rng = new Random(p_seed);
            return string.Concat(Enumerable.Range(0, p_length).Select(_ => "ACGT"[l_rng.Next(4)]));
        }

        /// <summary>
        /// 与えた配列から信頼できる k-mer 集合を組み立てる
        /// </summary>
        /// <param name="p_kmerLength">k 長</param>
        /// <param name="p_sequences">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex V_構築_索引(int p_kmerLength, params string[] p_sequences)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_kmerLength, A_スレッド数 = 1 };
            var l_作業 = Path.Combine(this._tempDir, Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_作業);
            var l_index = new TrustedKmerIndex(l_作業);
            foreach (var l_seq in p_sequences)
            {
                var l_bytes = l_seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + p_kmerLength <= l_bytes.Length; i++)
                {
                    for (var rep = 0; rep < 5; rep++)
                    {
                        l_index.V_登録(l_bytes.AsSpan(i, p_kmerLength), p_ワーカー番号: 0);
                    }
                }
            }
            _ = l_index.V_カットオフ(p_カットオフ: 2);
            return l_index;
        }

        /// <summary>
        /// 従来実装と転がし実装が、すべての開始点で同じ配列を返すこと
        /// </summary>
        private void V_両実装が一致する(int p_k長, params string[] p_配列)
        {
            using var l_索引 = this.V_構築_索引(p_k長, p_配列);
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

        /// <summary>
        /// 直鎖配列で、転がし実装と従来実装が一致する
        /// </summary>
        [Theory]
        [InlineData(21)]
        [InlineData(31)]
        [InlineData(33)]
        [InlineData(63)]
        public void 直鎖配列で転がし実装と従来実装が一致する(int p_kmerLength)
        {
            this.V_両実装が一致する(p_kmerLength, V_生成_乱数配列(3_000, p_seed: 801));
        }

        /// <summary>
        /// k=32 と k=33 は内部表現 (ulong と UInt128) の境界
        /// </summary>
        /// <remarks>
        /// 転がしのマスクとシフトがここで壊れやすい
        /// </remarks>
        [Theory]
        [InlineData(31)]
        [InlineData(32)]
        [InlineData(33)]
        [InlineData(64)]
        public void 内部表現の境界で転がし実装と従来実装が一致する(int p_kmerLength)
        {
            this.V_両実装が一致する(p_kmerLength, V_生成_乱数配列(2_000, p_seed: 802));
        }

        /// <summary>
        /// 分岐がある場合、両実装が同じ位置で walk を止めること
        /// </summary>
        [Fact]
        public void 分岐がある場合に転がし実装と従来実装が一致する()
        {
            var l_共通 = V_生成_乱数配列(1_000, p_seed: 811);
            var l_枝1 = V_生成_乱数配列(800, p_seed: 812);
            var l_枝2 = V_生成_乱数配列(800, p_seed: 813);

            this.V_両実装が一致する(31, l_共通 + l_枝1, l_共通 + l_枝2);
        }

        /// <summary>
        /// 反復配列を含む場合
        /// </summary>
        /// <remarks>
        /// 合流点の入次数判定が両実装で一致すること
        /// </remarks>
        [Fact]
        public void 反復配列を含む場合に転がし実装と従来実装が一致する()
        {
            var l_a = V_生成_乱数配列(800, p_seed: 821);
            var l_r = V_生成_乱数配列(300, p_seed: 822);
            var l_b = V_生成_乱数配列(800, p_seed: 823);
            var l_c = V_生成_乱数配列(800, p_seed: 824);

            this.V_両実装が一致する(31, l_a + l_r + l_b + l_r + l_c);
        }

        /// <summary>
        /// 環状配列
        /// </summary>
        /// <remarks>
        /// 循環検出の打ち切り位置が両実装で一致すること<br/>
        /// 完全な環には開始 k-mer が存在しない (どの k-mer も入次数 1 で、
        /// その予測元の出次数も 1) ため、任意の k-mer から walk して比べる
        /// </remarks>
        [Fact]
        public void 環状配列で転がし実装と従来実装が一致する()
        {
            const int k = 31;
            var l_環 = V_生成_乱数配列(1_500, p_seed: 831);
            var l_配列 = l_環 + l_環[..(k - 1)];

            using var l_索引 = this.V_構築_索引(k, l_配列);
            Assert.Empty(l_索引.Get_開始kmer一覧());

            var l_開始 = l_配列[..k].Select(Util.Get_塩基ID).ToArray();
            var l_期待 = new UnitigMaker(l_索引).Get_ユニティグ(l_開始).A_配列;
            var l_実際 = string.Concat(
                new UnitigWalk(l_索引, k).Get_塩基列(l_開始, []).Select(Util.Get_塩基文字));

            Assert.Equal(l_期待, l_実際);
        }

        /// <summary>
        /// 逆相補側から始めても一致すること (正規形の判定が転がしでも正しいこと)
        /// </summary>
        [Fact]
        public void 逆相補側から始めても転がし実装と従来実装が一致する()
        {
            var l_配列 = V_生成_乱数配列(2_000, p_seed: 841);
            this.V_両実装が一致する(31, l_配列, Util.V_逆相補(l_配列));
        }

        /// <summary>
        /// 内部表現の上限を超える k 長では、扱えないと判定される
        /// </summary>
        [Fact]
        public void 内部表現の上限を超えるk長では扱えないと判定される()
        {
            Assert.True(UnitigWalk.Get_扱えるか(64));
            Assert.False(UnitigWalk.Get_扱えるか(65));
        }
    }
}
