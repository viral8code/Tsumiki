using Tsumiki.Common;
using Tsumiki.Core.Output;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;

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
        private const int AmbiguousKmer = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int K = 8;

        // BeamSearchExtenderTests と同じ構成: A が B/C へ分岐する

        /// <summary>
        /// 分岐元
        /// </summary>
        private const string UnitigA = "TGGCAAGTCACTCTCGACCGA";

        /// <summary>
        /// 分岐先の片方
        /// </summary>
        private const string UnitigB = "CGACCGAACGGCGCCGGATC";

        /// <summary>
        /// 分岐先のもう片方
        /// </summary>
        private const string UnitigC = "CGACCGACTGTAATTCTACC";

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

        private static (List<string> UnitigList, UnitigGraph Graph) Build()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };
            List<string> unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> kmerDict = [];

            var id = 1;
            foreach (var seq in new[] { UnitigA, UnitigB, UnitigC })
            {
                unitigList.Add(seq);
                unitigList.Add(Util.V_逆相補(seq));
                for (var i = K; i <= seq.Length; i++)
                {
                    var startPos = i - K;
                    var key = new KmerKey(seq.AsSpan(startPos, K));
                    Register(kmerDict, key, id, startPos);
                    Register(kmerDict, key.Get_逆相補(), -id, seq.Length - i);
                }
                id++;
            }
            return (unitigList, UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer));
        }

        /// <summary>
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="dict">登録先の辞書</param>
        /// <param name="key">登録する k-mer</param>
        /// <param name="id">ユニティグ ID</param>
        /// <param name="position">ユニティグ内の開始位置</param>
        private static void Register(Dictionary<KmerKey, (int, int)> dict, KmerKey key, int id, int position)
        {
            if (dict.TryGetValue(key, out var existing))
            {
                if (existing.Item1 is AmbiguousKmer || existing.Item1 == id)
                {
                    return;
                }
                dict[key] = (AmbiguousKmer, 0);
                return;
            }
            dict[key] = (id, position);
        }

        [Fact]
        public void V_出力_WritesOneSegmentPerUnitigWithItsForwardSequenceAndLength()
        {
            var (unitigList, graph) = Build();
            var path = Path.Combine(this._tempDir, "graph.gfa");

            GfaWriter.V_出力(path, unitigList, graph, K);

            var lines = File.ReadAllLines(path);
            Assert.Equal("H\tVN:Z:1.0", lines[0]);

            var sLines = lines.Where(l => l.StartsWith("S\t")).ToList();
            Assert.Equal(3, sLines.Count);
            Assert.Contains($"S\t1\t{UnitigA}\tLN:i:{UnitigA.Length}", sLines);
            Assert.Contains($"S\t2\t{UnitigB}\tLN:i:{UnitigB.Length}", sLines);
            Assert.Contains($"S\t3\t{UnitigC}\tLN:i:{UnitigC.Length}", sLines);
        }

        [Fact]
        public void V_出力_WritesOneLinkPerPhysicalAdjacency_NotOnePerDirectedTwinPair()
        {
            var (unitigList, graph) = Build();
            var path = Path.Combine(this._tempDir, "graph_links.gfa");

            GfaWriter.V_出力(path, unitigList, graph, K);

            var lLines = File.ReadAllLines(path).Where(l => l.StartsWith("L\t")).ToList();

            // A は B・C の両方へ分岐する (2 つの物理的な隣接)
            // 各隣接は v→w と w^1→v^1 の双子として内部的には 2 回現れるが、
            // GFA には 1 本ずつしか出ないこと
            Assert.Equal(2, lLines.Count);
            Assert.Contains(lLines, l => l == $"L\t1\t+\t2\t+\t{K - 1}M");
            Assert.Contains(lLines, l => l == $"L\t1\t+\t3\t+\t{K - 1}M");
        }

        [Fact]
        public void V_出力_IncludesCopyNumberTag_WhenProvided()
        {
            var (unitigList, graph) = Build();
            var path = Path.Combine(this._tempDir, "graph_cn.gfa");
            Dictionary<int, int> copyNumber = new() { [1] = 1, [2] = 2, [3] = 1 };

            GfaWriter.V_出力(path, unitigList, graph, K, copyNumber);

            var sLines = File.ReadAllLines(path).Where(l => l.StartsWith("S\t")).ToList();
            Assert.Contains(sLines, l => l.StartsWith("S\t2\t") && l.EndsWith("CN:i:2"));
            Assert.Contains(sLines, l => l.StartsWith("S\t1\t") && l.EndsWith("CN:i:1"));
        }

        [Fact]
        public void V_出力_OmitsCopyNumberTag_WhenNotProvided()
        {
            var (unitigList, graph) = Build();
            var path = Path.Combine(this._tempDir, "graph_nocn.gfa");

            GfaWriter.V_出力(path, unitigList, graph, K);

            var sLines = File.ReadAllLines(path).Where(l => l.StartsWith("S\t")).ToList();
            Assert.DoesNotContain(sLines, l => l.Contains("CN:i:"));
        }
    }
}
