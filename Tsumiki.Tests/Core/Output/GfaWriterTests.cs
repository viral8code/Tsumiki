using Tsumiki.Commons;
using Tsumiki.Cores.Output;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// unitig グラフの GFA1 出力の検証
    /// </summary>
    /// <remarks>
    /// 決められない分岐は今まで打ち切り点になるだけで理由が出力に残らな
    /// かった<br/>
    /// GFA として書き出せば、Bandage 等のビューアでグラフの形が
    /// 直接見えるようになる
    /// </remarks>
    public class GfaWriterTests : IDisposable
    {
        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int 曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int K = 8;

        // BeamSearchExtenderTests と同じ構成: A が B/C へ分岐する

        /// <summary>
        /// 分岐元
        /// </summary>
        private const string ユニティグA = "TGGCAAGTCACTCTCGACCGA";

        /// <summary>
        /// 分岐先の片方
        /// </summary>
        private const string ユニティグB = "CGACCGAACGGCGCCGGATC";

        /// <summary>
        /// 分岐先のもう片方
        /// </summary>
        private const string ユニティグC = "CGACCGACTGTAATTCTACC";

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public GfaWriterTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_gfa_tests_" + Guid.NewGuid().ToString("N"));
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
        /// 分岐を持つ検証用のユニティググラフを組み立てる
        /// </summary>
        /// <returns>ユニティグ一覧とグラフ</returns>
        private static (List<string> A_ユニティグ一覧, UnitigGraph A_グラフ) V_構築()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };
            List<string> l_ユニティグ一覧 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> l_kmer辞書 = [];

            var l_id = 1;
            foreach (var l_seq in new[] { ユニティグA, ユニティグB, ユニティグC })
            {
                l_ユニティグ一覧.Add(l_seq);
                l_ユニティグ一覧.Add(Util.V_逆相補(l_seq));
                for (var i = K; i <= l_seq.Length; i++)
                {
                    var l_開始位置 = i - K;
                    var l_key = new KmerKey(l_seq.AsSpan(l_開始位置, K));
                    V_登録(l_kmer辞書, l_key, l_id, l_開始位置);
                    V_登録(l_kmer辞書, l_key.Get_逆相補(), -l_id, l_seq.Length - i);
                }
                l_id++;
            }
            return (l_ユニティグ一覧, UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, K, 曖昧kmer番号));
        }

        /// <summary>
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_key">登録する k-mer</param>
        /// <param name="p_id">ユニティグ ID</param>
        /// <param name="p_位置">ユニティグ内の開始位置</param>
        private static void V_登録(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_key, int p_id, int p_位置)
        {
            if (p_辞書.TryGetValue(p_key, out var l_既存))
            {
                if (l_既存.Item1 is 曖昧kmer番号 || l_既存.Item1 == p_id)
                {
                    return;
                }
                p_辞書[p_key] = (曖昧kmer番号, 0);
                return;
            }
            p_辞書[p_key] = (p_id, p_位置);
        }

        /// <summary>
        /// ユニティグごとに順方向配列と長さを持つ S 行を 1 つ書くことを確かめる
        /// </summary>
        [Fact]
        public void V_出力_ユニティグごとに順方向配列と長さを持つS行を1つ書く()
        {
            var (l_ユニティグ一覧, l_グラフ) = V_構築();
            var l_パス = Path.Combine(this._tempDir, "graph.gfa");

            GfaWriter.V_出力(l_パス, l_ユニティグ一覧, l_グラフ, K);

            var l_行 = File.ReadAllLines(l_パス);
            Assert.Equal("H\tVN:Z:1.0", l_行[0]);

            var l_S行 = l_行.Where(l => l.StartsWith("S\t")).ToList();
            Assert.Equal(3, l_S行.Count);
            Assert.Contains($"S\t1\t{ユニティグA}\tLN:i:{ユニティグA.Length}", l_S行);
            Assert.Contains($"S\t2\t{ユニティグB}\tLN:i:{ユニティグB.Length}", l_S行);
            Assert.Contains($"S\t3\t{ユニティグC}\tLN:i:{ユニティグC.Length}", l_S行);
        }

        /// <summary>
        /// 物理的な隣接ごとに L 行を 1 本だけ書くことを確かめる
        /// </summary>
        [Fact]
        public void V_出力_物理的な隣接ごとにL行を1本だけ書く()
        {
            var (l_ユニティグ一覧, l_グラフ) = V_構築();
            var l_パス = Path.Combine(this._tempDir, "graph_links.gfa");

            GfaWriter.V_出力(l_パス, l_ユニティグ一覧, l_グラフ, K);

            var l_L行 = File.ReadAllLines(l_パス).Where(l => l.StartsWith("L\t")).ToList();

            // A は B・C の両方へ分岐する (2 つの物理的な隣接)
            // 各隣接は v→w と w^1→v^1 の双子として内部的には 2 回現れるが、
            // GFA には 1 本ずつしか出ないこと
            Assert.Equal(2, l_L行.Count);
            Assert.Contains(l_L行, l => l == $"L\t1\t+\t2\t+\t{K - 1}M");
            Assert.Contains(l_L行, l => l == $"L\t1\t+\t3\t+\t{K - 1}M");
        }

        /// <summary>
        /// コピー数を渡せば CN タグを含めることを確かめる
        /// </summary>
        [Fact]
        public void V_出力_コピー数を渡せばCNタグを含める()
        {
            var (l_ユニティグ一覧, l_グラフ) = V_構築();
            var l_パス = Path.Combine(this._tempDir, "graph_cn.gfa");
            Dictionary<int, int> l_コピー数 = new() { [1] = 1, [2] = 2, [3] = 1 };

            GfaWriter.V_出力(l_パス, l_ユニティグ一覧, l_グラフ, K, l_コピー数);

            var l_S行 = File.ReadAllLines(l_パス).Where(l => l.StartsWith("S\t")).ToList();
            Assert.Contains(l_S行, l => l.StartsWith("S\t2\t") && l.EndsWith("CN:i:2"));
            Assert.Contains(l_S行, l => l.StartsWith("S\t1\t") && l.EndsWith("CN:i:1"));
        }

        /// <summary>
        /// コピー数を渡さなければ CN タグを省くことを確かめる
        /// </summary>
        [Fact]
        public void V_出力_コピー数を渡さなければCNタグを省く()
        {
            var (l_ユニティグ一覧, l_グラフ) = V_構築();
            var l_パス = Path.Combine(this._tempDir, "graph_nocn.gfa");

            GfaWriter.V_出力(l_パス, l_ユニティグ一覧, l_グラフ, K);

            var l_S行 = File.ReadAllLines(l_パス).Where(l => l.StartsWith("S\t")).ToList();
            Assert.DoesNotContain(l_S行, l => l.Contains("CN:i:"));
        }
    }
}
