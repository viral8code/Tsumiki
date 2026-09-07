using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 提案D: BeamSearchExtender のスコアを生カウントではなく期待本数との比で
    /// 測ることの検証。
    ///
    /// 分岐元 A から、短い unitig B と長い unitig C へ分岐する構成を作る。
    /// C は B よりずっと長いため、同じフラグメント長分布のもとでは
    /// 「両端が収まる開始位置」の窓が B よりずっと広い(=同じ観測本数でも
    /// 期待本数は C のほうが大きい)。生カウントでは C がわずかに優勢に
    /// 見えても優勢閾値を超えないケースで、期待本数との比を取ると
    /// 実際には B のほうが強い証拠であるとして正しく選ばれることを確認する。
    /// </summary>
    public class BeamSearchExtenderCalibrationTests
    {
        private const int AmbiguousKmer = int.MinValue;
        private const int K = 21;
        private const int AnchorLength = K - 1;

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
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

        /// <summary>
        /// A(250bp) が分岐元、B(35bp, 短い)と C(2020bp, 長い)がその行き先。
        /// A の末尾20塩基(=k-1)を B・C 両方の先頭が共有することで分岐にする。
        /// </summary>
        private static (List<string> UnitigList, UnitigGraph Graph, int A, int B, int C) Build()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };

            var anchor = RandomSequence(AnchorLength, seed: 20260908);
            var unitigA = RandomSequence(230, seed: 1) + anchor; // 250bp
            var unitigB = anchor + RandomSequence(15, seed: 2); // 35bp(短い)
            var unitigC = anchor + RandomSequence(2000, seed: 3); // 2020bp(長い)

            List<string> unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> kmerDict = [];

            var id = 1;
            foreach (var seq in new[] { unitigA, unitigB, unitigC })
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

            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);
            return (unitigList, graph, 1, 2, 3);
        }

        private static int[] NoMerges(UnitigGraph graph)
        {
            var merge = new int[graph.A_出辺.Count];
            Array.Fill(merge, -1);
            return merge;
        }

        /// <summary>中央 150 のフラグメント長標本(密度較正に使う)。</summary>
        private static List<int> Get_同一ユニティグ標本(int p_件数 = 300) => [.. Enumerable.Repeat(150, p_件数)];

        [Fact]
        public void Extend_WithoutCalibration_RawCountsDoNotDominateEnoughToDecide()
        {
            var (unitigList, graph, aId, bId, cId) = Build();
            var a = ContigMaker.Get_頂点番号(aId);
            var b = ContigMaker.Get_頂点番号(bId);
            var c = ContigMaker.Get_頂点番号(cId);

            // 生カウントでは C(4) が B(3) よりわずかに多いが、
            // 優勢比 4/7=0.571 は閾値0.8を超えない。
            Dictionary<(int, int), ulong> pairLink = new() { [(a, b)] = 3, [(a, c)] = 4 };
            Dictionary<int, int> copyNumber = new() { [aId] = 1, [bId] = 1, [cId] = 1 };

            var merge = NoMerges(graph);
            _ = BeamSearchExtender.V_延長_先読み(
                graph, unitigList, merge, pairLink, copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8m, p_最小証拠数: 3, p_較正器: null);

            Assert.Equal(-1, merge[a]);
        }

        /// <summary>
        /// 同じ観測本数(3 vs 4)でも、C は B よりずっと長いぶん期待本数も
        /// 大きい。期待本数との比を取ると B(短い)のほうが実際には強い証拠で
        /// あるとわかり、優勢閾値を超えて A→B が選ばれるはず。
        /// </summary>
        [Fact]
        public void Extend_WithCalibration_PrefersTheShortFlank_ThatRawCountsCouldNotDecide()
        {
            var (unitigList, graph, aId, bId, cId) = Build();
            var a = ContigMaker.Get_頂点番号(aId);
            var b = ContigMaker.Get_頂点番号(bId);
            var c = ContigMaker.Get_頂点番号(cId);

            Dictionary<(int, int), ulong> pairLink = new() { [(a, b)] = 3, [(a, c)] = 4 };
            Dictionary<int, int> copyNumber = new() { [aId] = 1, [bId] = 1, [cId] = 1 };

            var 較正器 = 証拠較正器.Get_較正器(
                Get_同一ユニティグ標本(), p_リード長: 30,
                [(long)unitigList[a].Length, (long)unitigList[b].Length, (long)unitigList[c].Length]);
            Assert.True(較正器.A_使えるか);

            var merge = NoMerges(graph);
            var committed = BeamSearchExtender.V_延長_先読み(
                graph, unitigList, merge, pairLink, copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8m, p_最小証拠数: 3, p_較正器: 較正器);

            Assert.True(committed > 0, "calibrated lookahead should have resolved the junction toward the short flank");
            Assert.Equal(b, merge[a]);
            Assert.Equal(a ^ 1, merge[b ^ 1]);
        }
    }
}
