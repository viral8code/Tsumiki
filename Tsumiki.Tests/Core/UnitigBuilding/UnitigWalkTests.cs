using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// パック値を転がしながら進める walk が、従来の実装と同じ結果を返すことを固定する
    /// </summary>
    /// <remarks>
    /// 転がし更新は unitig 構築の時間のほとんどを占めていた O (k) の詰め直しを省くためのもので、結果は 1 塩基たりとも変わってはいけない<br/>
    /// 2 つの実装が並存する以上、等価性の確認は必須になる
    /// </remarks>
    public class UnitigWalkTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public UnitigWalkTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_walk_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 直鎖配列で、転がし実装と従来実装が一致する
        /// </summary>
        /// <param name="p_kmer長"></param>
        [Theory]
        [InlineData(21)]
        [InlineData(31)]
        [InlineData(33)]
        [InlineData(63)]
        public void V_直鎖配列で転がし実装と従来実装が一致する(int p_kmer長)
        {
            this.V_両実装が一致する(p_kmer長, V_生成_乱数配列(3_000, p_乱数種: 801));
        }

        /// <summary>
        /// k=32 と k=33 は内部表現 (ulong と UInt128) の境界
        /// </summary>
        /// <remarks>
        /// 転がしのマスクとシフトがここで壊れやすい
        /// </remarks>
        /// <param name="p_kmer長"></param>
        [Theory]
        [InlineData(31)]
        [InlineData(32)]
        [InlineData(33)]
        [InlineData(64)]
        public void V_内部表現の境界で転がし実装と従来実装が一致する(int p_kmer長)
        {
            this.V_両実装が一致する(p_kmer長, V_生成_乱数配列(2_000, p_乱数種: 802));
        }

        /// <summary>
        /// 分岐がある場合、両実装が同じ位置で walk を止めること
        /// </summary>
        [Fact]
        public void V_分岐がある場合に転がし実装と従来実装が一致する()
        {
            var l_共通 = V_生成_乱数配列(1_000, p_乱数種: 811);
            var l_枝1 = V_生成_乱数配列(800, p_乱数種: 812);
            var l_枝2 = V_生成_乱数配列(800, p_乱数種: 813);

            this.V_両実装が一致する(31, l_共通 + l_枝1, l_共通 + l_枝2);
        }

        /// <summary>
        /// 反復配列を含む場合
        /// </summary>
        /// <remarks>
        /// 合流点の入次数判定が両実装で一致すること
        /// </remarks>
        [Fact]
        public void V_反復配列を含む場合に転がし実装と従来実装が一致する()
        {
            var l_先頭配列 = V_生成_乱数配列(800, p_乱数種: 821);
            var l_反復配列 = V_生成_乱数配列(300, p_乱数種: 822);
            var l_中間配列 = V_生成_乱数配列(800, p_乱数種: 823);
            var l_末尾配列 = V_生成_乱数配列(800, p_乱数種: 824);

            this.V_両実装が一致する(31, l_先頭配列 + l_反復配列 + l_中間配列 + l_反復配列 + l_末尾配列);
        }

        /// <summary>
        /// 環状配列
        /// </summary>
        /// <remarks>
        /// 循環検出の打ち切り位置が両実装で一致すること<br/>
        /// 完全な環には開始 k-mer が存在しない (どの k-mer も入次数 1 で、その予測元の出次数も 1) ため、任意の k-mer から walk して比べる
        /// </remarks>
        [Fact]
        public void V_環状配列で転がし実装と従来実装が一致する()
        {
            const int l_k長 = 31;
            var l_環 = V_生成_乱数配列(1_500, p_乱数種: 831);
            var l_配列 = l_環 + l_環[..(l_k長 - 1)];

            using var l_索引 = this.V_構築_索引(l_k長, l_配列);
            Assert.Empty(l_索引.Get_開始kmer一覧());

            var l_開始 = l_配列[..l_k長].Select(Util.Get_塩基ID).ToArray();
            var l_期待 = new UnitigMaker(l_索引).Get_ユニティグ(l_開始).A_配列;
            var l_実際 = string.Concat(new UnitigWalk(l_索引, l_k長).Get_塩基列(l_開始, []).Select(Util.Get_塩基文字));

            Assert.Equal(l_期待, l_実際);
        }

        /// <summary>
        /// 逆相補側から始めても一致すること (正規形の判定が転がしでも正しいこと)
        /// </summary>
        [Fact]
        public void V_逆相補側から始めても転がし実装と従来実装が一致する()
        {
            var l_配列 = V_生成_乱数配列(2_000, p_乱数種: 841);
            this.V_両実装が一致する(31, l_配列, Util.V_逆相補(l_配列));
        }

        /// <summary>
        /// 内部表現の上限を超える k 長では、扱えないと判定される
        /// </summary>
        [Fact]
        public void V_内部表現の上限を超えるk長では扱えないと判定される()
        {
            Assert.True(UnitigWalk.Get_扱えるか(64));
            Assert.False(UnitigWalk.Get_扱えるか(65));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_乱数種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_長さ, int p_乱数種)
        {
            var l_乱数生成器 = new Random(p_乱数種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数生成器.Next(4)]));
        }

        /// <summary>
        /// 与えた配列から信頼できる k-mer 集合を組み立てる
        /// </summary>
        /// <param name="p_kmer長">k 長</param>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        private TrustedKmerIndex V_構築_索引(int p_kmer長, params string[] p_配列群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_kmer長, A_スレッド数 = 1 };
            var l_作業 = Path.Combine(this._作業ディレクトリ, Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_作業);
            var l_インデックス = new TrustedKmerIndex(l_作業);
            foreach (var l_配列 in p_配列群)
            {
                var l_バイト列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + p_kmer長 <= l_バイト列.Length; i++)
                {
                    for (var l_反復回数 = 0; l_反復回数 < 5; l_反復回数++)
                    {
                        l_インデックス.V_登録(l_バイト列.AsSpan(i, p_kmer長));
                    }
                }
            }
            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);
            return l_インデックス;
        }

        /// <summary>
        /// 従来実装と転がし実装が、すべての開始点で同じ配列を返すこと
        /// </summary>
        /// <param name="p_k長"></param>
        /// <param name="p_配列"></param>
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
                var l_実際 = string.Concat(l_転がし.Get_塩基列(l_開始, l_訪問済み).Select(Util.Get_塩基文字));
                Assert.Equal(l_期待, l_実際);
            }
        }

        #endregion

    }
}
