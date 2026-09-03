using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 先読み(ビームサーチ)による分岐解決の検証。
    ///
    /// 相互一意性の判定は「その1歩だけ」を見るため、分岐の直後だけを見ると
    /// 五分五分に見えるが、2〜3本先まで進めると片方だけがペアエンドの証拠と
    /// 整合する、という状況を取りこぼす。ここでは
    ///   A →(B or C)、B → D、C → E
    /// という形で、A の直後には証拠が無く D の位置に初めて証拠が現れる構成を作り、
    /// 先読みによって A → B が選ばれることを確認する。
    /// </summary>
    public class BeamSearchExtenderTests
    {
        private const int AmbiguousKmer = int.MinValue;
        private const int K = 8;

        // k=8 で5本すべてを通じて重複する正規化 k-mer が無いことを確認済みの構成。
        private const string UnitigA = "TGGCAAGTCACTCTCGACCGA";
        private const string UnitigB = "CGACCGAACGGCGCCGGATC";
        private const string UnitigC = "CGACCGACTGTAATTCTACC";
        private const string UnitigD = "CCGGATCAAAGCCACGGCTAG";
        private const string UnitigE = "TTCTACCAAAGGCTAGTATGA";

        private static (List<string> UnitigList, UnitigGraph Graph) Build()
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = K, ThreadCount = 1 };
            List<string> unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> kmerDict = [];

            var id = 1;
            foreach (var seq in new[] { UnitigA, UnitigB, UnitigC, UnitigD, UnitigE })
            {
                unitigList.Add(seq);
                unitigList.Add(Util.ReverseComprement(seq));
                for (var i = K; i <= seq.Length; i++)
                {
                    var startPos = i - K;
                    var key = new KmerKey(seq.AsSpan(startPos, K));
                    Register(kmerDict, key, id, startPos);
                    Register(kmerDict, key.ReverseComprement(), -id, seq.Length - i);
                }
                id++;
            }
            return (unitigList, UnitigGraph.Build(unitigList, kmerDict, K, AmbiguousKmer));
        }

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

        private static int[] NoMerges(UnitigGraph graph)
        {
            var merge = new int[graph.VertexCount];
            Array.Fill(merge, -1);
            return merge;
        }

        [Fact]
        public void Extend_ResolvesABranch_WhenTheEvidenceOnlyAppearsOneStepLater()
        {
            var (unitigList, graph) = Build();
            var a = ContigMaker.VertexIndex(1);
            var b = ContigMaker.VertexIndex(2);
            var d = ContigMaker.VertexIndex(4);

            // 前提: A は B と C の両方へ伸びられる(1歩だけでは決められない)。
            Assert.Equal(2, graph.OutEdges[a].Count);

            // 証拠は A の直後(B/C)ではなく、その次の D に現れる。
            Dictionary<(int, int), ulong> pairLink = new() { [(a, d)] = 30 };
            Dictionary<int, int> copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var merge = NoMerges(graph);
            var committed = BeamSearchExtender.Extend(
                graph, unitigList, merge, pairLink, copyNumber,
                insertSize: 400, dominanceThreshold: 0.8m, minimumEvidence: 5);

            Assert.True(committed > 0, "lookahead should have resolved at least one junction");
            Assert.Equal(b, merge[a]);
            // 逆鎖側も対称に設定されていること。
            Assert.Equal(a ^ 1, merge[b ^ 1]);
        }

        /// <summary>
        /// どちらの枝にも同程度の証拠がある場合は、僅差で選ばずに繋がない。
        /// ビームサーチの利点は「広く探して有力な仮説が一致する部分にだけ
        /// コミットする」ことにあり、五分五分の分岐で1本を選ぶことではない。
        /// </summary>
        [Fact]
        public void Extend_DoesNothing_WhenBothBranchesAreEquallySupported()
        {
            var (unitigList, graph) = Build();
            var a = ContigMaker.VertexIndex(1);
            var d = ContigMaker.VertexIndex(4);
            var e = ContigMaker.VertexIndex(5);

            Dictionary<(int, int), ulong> pairLink = new() { [(a, d)] = 20, [(a, e)] = 19 };
            Dictionary<int, int> copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var merge = NoMerges(graph);
            _ = BeamSearchExtender.Extend(
                graph, unitigList, merge, pairLink, copyNumber,
                insertSize: 400, dominanceThreshold: 0.8m, minimumEvidence: 5);

            Assert.Equal(-1, merge[a]);
        }

        /// <summary>
        /// ペアエンドの証拠がまったく無ければ、根拠が無いので繋がない。
        /// </summary>
        [Fact]
        public void Extend_DoesNothing_WhenThereIsNoPairEvidenceAtAll()
        {
            var (unitigList, graph) = Build();
            var a = ContigMaker.VertexIndex(1);

            Dictionary<(int, int), ulong> pairLink = [];
            Dictionary<int, int> copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var merge = NoMerges(graph);
            _ = BeamSearchExtender.Extend(
                graph, unitigList, merge, pairLink, copyNumber,
                insertSize: 400, dominanceThreshold: 0.8m, minimumEvidence: 5);

            Assert.Equal(-1, merge[a]);
        }

        /// <summary>
        /// 証拠はあるが少なすぎる場合、偶然の一致で繋いでしまわないよう見送る。
        /// </summary>
        [Fact]
        public void Extend_DoesNothing_WhenEvidenceIsBelowTheMinimum()
        {
            var (unitigList, graph) = Build();
            var a = ContigMaker.VertexIndex(1);
            var d = ContigMaker.VertexIndex(4);

            Dictionary<(int, int), ulong> pairLink = new() { [(a, d)] = 2 };
            Dictionary<int, int> copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var merge = NoMerges(graph);
            _ = BeamSearchExtender.Extend(
                graph, unitigList, merge, pairLink, copyNumber,
                insertSize: 400, dominanceThreshold: 0.8m, minimumEvidence: 10);

            Assert.Equal(-1, merge[a]);
        }

        /// <summary>
        /// 既に別の結合が入っている行き先へは、それを壊してまで繋がない
        /// (相互一意性を保つ)。
        /// </summary>
        [Fact]
        public void Extend_DoesNotStealATargetThatAlreadyHasAnIncomingMerge()
        {
            var (unitigList, graph) = Build();
            var a = ContigMaker.VertexIndex(1);
            var b = ContigMaker.VertexIndex(2);
            var d = ContigMaker.VertexIndex(4);

            Dictionary<(int, int), ulong> pairLink = new() { [(a, d)] = 30 };
            Dictionary<int, int> copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var merge = NoMerges(graph);
            // B には既に(別の経路からの)結合が入っていることにする。
            merge[b ^ 1] = ContigMaker.VertexIndex(5) ^ 1;

            _ = BeamSearchExtender.Extend(
                graph, unitigList, merge, pairLink, copyNumber,
                insertSize: 400, dominanceThreshold: 0.8m, minimumEvidence: 5);

            Assert.Equal(-1, merge[a]);
        }
    }
}
