using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// カットオフを掛ける前に出現回数ヒストグラムだけを取り出せること。
    /// -kc を自動決定するには、カットオフを決める前にスペクトルを見る必要がある。
    ///
    /// 事前走査と本体の走査は同じ統合ファイルを使い回す。統合をやり直すと
    /// マージソートのディスク I/O が丸ごと二重になるため。
    /// </summary>
    public class TrustedKmerIndexHistogramTests : IDisposable
    {
        private readonly string _tempDir;

        public TrustedKmerIndexHistogramTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_indexhist_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private const int K = 21;

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        /// <summary>
        /// 位置 i の k-mer を (i % 4) + 1 回登録し、その分布がそのまま
        /// ヒストグラムに現れることを確かめる。
        /// </summary>
        private TrustedKmerIndex BuildIndex(out Dictionary<ulong, long> p_期待ヒストグラム)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 4 };
            var index = new TrustedKmerIndex(this._tempDir);

            var bytes = RandomSequence(400, seed: 4242).Select(Util.Get_塩基ID).ToArray();
            p_期待ヒストグラム = [];
            for (var i = 0; i + K <= bytes.Length; i++)
            {
                var times = (ulong)((i % 4) + 1);
                for (var t = 0UL; t < times; t++)
                {
                    index.V_登録(bytes.AsSpan(i, K), p_ワーカー番号: (int)(t % 4));
                }
                p_期待ヒストグラム[times] = p_期待ヒストグラム.GetValueOrDefault(times, 0L) + 1;
            }
            return index;
        }

        [Fact]
        public void GetHistogram_BeforeCutoff_MatchesTheRegisteredMultiplicities()
        {
            using var index = this.BuildIndex(out var expected);

            var histogram = index.Get_出現回数ヒストグラム();

            Assert.Equal(expected, histogram);
        }

        /// <summary>
        /// 事前走査でヒストグラムを取っても、その後のカットオフが正しく動くこと。
        /// 統合ファイルを使い回す実装なので、事前走査がファイルを消費・削除して
        /// しまうと本体が壊れる。
        /// </summary>
        [Fact]
        public void GetHistogram_DoesNotConsumeTheMergedFile_CutoffStillWorks()
        {
            using var index = this.BuildIndex(out var expected);

            var 事前 = index.Get_出現回数ヒストグラム();
            _ = index.V_カットオフ(p_カットオフ: 3);

            // カットオフ本体が集計したヒストグラムも同じでなければならない。
            Assert.Equal(事前, index.A_出現回数ヒストグラム);
            Assert.Equal(expected, index.A_出現回数ヒストグラム);
        }

        /// <summary>
        /// 事前走査を挟まなかった場合も A_出現回数ヒストグラム は埋まること
        /// (-kc 明示指定時はこちらの経路しか通らない)。
        /// </summary>
        [Fact]
        public void Cutoff_WithoutAPriorScan_StillRecordsTheHistogram()
        {
            using var index = this.BuildIndex(out var expected);

            _ = index.V_カットオフ(p_カットオフ: 2);

            Assert.Equal(expected, index.A_出現回数ヒストグラム);
        }
    }
}
