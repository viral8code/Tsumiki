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
    /// 統合をやり直すと
    /// マージソートのディスク I/O が丸ごと二重になるため
    /// </remarks>
    public class TrustedKmerIndexHistogramTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public TrustedKmerIndexHistogramTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_indexhist_tests_" + Guid.NewGuid().ToString("N"));
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
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 21;

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
        /// 位置 i の k-mer を (i % 4) + 1 回登録し、その分布がそのまま
        /// ヒストグラムに現れることを確かめる
        /// </summary>
        private TrustedKmerIndex V_構築_索引(out Dictionary<ulong, long> p_期待ヒストグラム)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 4 };
            var l_index = new TrustedKmerIndex(this._tempDir);

            var l_bytes = V_生成_乱数配列(400, p_seed: 4242).Select(Util.Get_塩基ID).ToArray();
            p_期待ヒストグラム = [];
            for (var i = 0; i + k長 <= l_bytes.Length; i++)
            {
                var l_times = (ulong)((i % 4) + 1);
                for (var t = 0UL; t < l_times; t++)
                {
                    l_index.V_登録(l_bytes.AsSpan(i, k長), p_ワーカー番号: (int)(t % 4));
                }
                p_期待ヒストグラム[l_times] = p_期待ヒストグラム.GetValueOrDefault(l_times, 0L) + 1;
            }
            return l_index;
        }

        /// <summary>
        /// カットオフ前のヒストグラムは、登録した回数と一致する
        /// </summary>
        [Fact]
        public void カットオフ前のヒストグラムは登録した回数と一致する()
        {
            using var l_index = this.V_構築_索引(out var l_expected);

            var l_histogram = l_index.Get_出現回数ヒストグラム();

            Assert.Equal(l_expected, l_histogram);
        }

        /// <summary>
        /// 事前走査でヒストグラムを取っても、その後のカットオフが正しく動くこと
        /// </summary>
        /// <remarks>
        /// 統合ファイルを使い回す実装なので、事前走査がファイルを消費・削除して
        /// しまうと本体が壊れる
        /// </remarks>
        [Fact]
        public void ヒストグラム取得は統合ファイルを消費せずカットオフも動く()
        {
            using var l_index = this.V_構築_索引(out var l_expected);

            var l_事前 = l_index.Get_出現回数ヒストグラム();
            _ = l_index.V_カットオフ(p_カットオフ: 3);

            // カットオフ本体が集計したヒストグラムも同じでなければならない
            Assert.Equal(l_事前, l_index.A_出現回数ヒストグラム);
            Assert.Equal(l_expected, l_index.A_出現回数ヒストグラム);
        }

        /// <summary>
        /// 事前走査を挟まなかった場合も A_出現回数ヒストグラム は埋まること
        /// (-kc 明示指定時はこちらの経路しか通らない)
        /// </summary>
        [Fact]
        public void 事前走査を挟まなくてもヒストグラムは記録される()
        {
            using var l_index = this.V_構築_索引(out var l_expected);

            _ = l_index.V_カットオフ(p_カットオフ: 2);

            Assert.Equal(l_expected, l_index.A_出現回数ヒストグラム);
        }
    }
}
