using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// カットオフを掛ける前に出現回数ヒストグラムだけを取り出せること
    /// </summary>
    /// <remarks>
    /// -kc を自動決定するには、カットオフを決める前にスペクトルを見る必要がある<br/>
    /// 事前走査と本体の走査は同じ統合ファイルを使い回す<br/>
    /// 統合をやり直すとマージソートのディスク I/O が丸ごと二重になるため
    /// </remarks>
    public class TrustedKmerIndexHistogramTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 21;

        #endregion

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
        public TrustedKmerIndexHistogramTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_indexhist_tests_" + Guid.NewGuid().ToString("N"));
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
        /// カットオフ前のヒストグラムは、登録した回数と一致する
        /// </summary>
        [Fact]
        public void V_カットオフ前のヒストグラムは登録した回数と一致する()
        {
            using var l_インデックス = this.V_構築_索引(out var l_期待値);

            var l_histogram = l_インデックス.Get_出現回数ヒストグラム();

            Assert.Equal(l_期待値, l_histogram);
        }

        /// <summary>
        /// 事前走査でヒストグラムを取っても、その後のカットオフが正しく動くこと
        /// </summary>
        /// <remarks>
        /// 統合ファイルを使い回す実装なので、事前走査がファイルを消費・削除してしまうと本体が壊れる
        /// </remarks>
        [Fact]
        public void V_ヒストグラム取得は統合ファイルを消費せずカットオフも動く()
        {
            using var l_インデックス = this.V_構築_索引(out var l_期待値);

            var l_事前 = l_インデックス.Get_出現回数ヒストグラム();
            _ = l_インデックス.V_カットオフ(p_カットオフ: 3UL);

            // カットオフ本体が集計したヒストグラムも同じでなければならない
            Assert.Equal(l_事前, l_インデックス.A_出現回数ヒストグラム);
            Assert.Equal(l_期待値, l_インデックス.A_出現回数ヒストグラム);
        }

        /// <summary>
        /// 事前走査を挟まなかった場合も A_出現回数ヒストグラム は埋まること (-kc 明示指定時はこちらの経路しか通らない)
        /// </summary>
        [Fact]
        public void V_事前走査を挟まなくてもヒストグラムは記録される()
        {
            using var l_インデックス = this.V_構築_索引(out var l_期待値);

            _ = l_インデックス.V_カットオフ(p_カットオフ: 2UL);

            Assert.Equal(l_期待値, l_インデックス.A_出現回数ヒストグラム);
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
        /// 位置 i の k-mer を (i % 4) + 1 回登録し、その分布がそのままヒストグラムに現れることを確かめる
        /// </summary>
        /// <param name="p_期待ヒストグラム"></param>
        /// <returns></returns>
        private TrustedKmerIndex V_構築_索引(out Dictionary<ulong, long> p_期待ヒストグラム)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 4 };
            var l_インデックス = new TrustedKmerIndex(this._作業ディレクトリ);

            var l_バイト列 = V_生成_乱数配列(400, p_乱数種: 4_242).Select(Util.Get_塩基ID).ToArray();
            p_期待ヒストグラム = [];
            for (var i = 0; i + k長 <= l_バイト列.Length; i++)
            {
                var l_times = (ulong)((i % 4) + 1);
                for (var t = 0UL; t < l_times; t++)
                {
                    l_インデックス.V_登録(l_バイト列.AsSpan(i, k長));
                }
                p_期待ヒストグラム[l_times] = p_期待ヒストグラム.GetValueOrDefault(l_times, 0L) + 1L;
            }
            return l_インデックス;
        }

        #endregion

    }
}
