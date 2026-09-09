using Tsumiki.Common;
using Tsumiki.Core.Evidence;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 提案 D: BeamSearchExtender のスコアを生カウントではなく期待本数との比で
    /// 測ることの検証
    /// </summary>
    /// <remarks>
    /// 分岐元 A から、短い unitig B と長い unitig C へ分岐する構成を作る<br/>
    /// C は B よりずっと長いため、同じフラグメント長分布のもとでは
    /// 「両端が収まる開始位置」の窓が B よりずっと広い (=同じ観測本数でも
    /// 期待本数は C のほうが大きい)<br/>
    /// 生カウントでは C がわずかに優勢に
    /// 見えても優勢閾値を超えないケースで、期待本数との比を取ると
    /// 実際には B のほうが強い証拠であるとして正しく選ばれることを確認する
    /// </remarks>
    public class BeamSearchExtenderCalibrationTests
    {
        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int 曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 21;

        /// <summary>
        /// アンカーとして使う長さ
        /// </summary>
        private const int アンカー長 = k長 - 1;

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
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_dict">登録先の辞書</param>
        /// <param name="p_key">登録する k-mer</param>
        /// <param name="p_id">ユニティグ ID</param>
        /// <param name="p_position">ユニティグ内の開始位置</param>
        private static void V_登録_kmer(Dictionary<KmerKey, (int, int)> p_dict, KmerKey p_key, int p_id, int p_position)
        {
            if (p_dict.TryGetValue(p_key, out var l_existing))
            {
                if (l_existing.Item1 is 曖昧kmer番号 || l_existing.Item1 == p_id)
                {
                    return;
                }
                p_dict[p_key] = (曖昧kmer番号, 0);
                return;
            }
            p_dict[p_key] = (p_id, p_position);
        }

        /// <summary>
        /// A(250 bp) が分岐元、B(35 bp, 短い) と C(2020 bp, 長い) がその行き先
        /// </summary>
        /// <remarks>
        /// A の末尾 20 塩基 (=k-1) を B・C 両方の先頭が共有することで分岐にする
        /// </remarks>
        private static (List<string> A_ユニティグ一覧, UnitigGraph A_グラフ, int A_分岐元ID, int A_短い方ID, int A_長い方ID) V_構築()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };

            var l_anchor = V_生成_乱数配列(アンカー長, p_seed: 20260908);
            var l_unitigA = V_生成_乱数配列(230, p_seed: 1) + l_anchor; // 250bp
            var l_unitigB = l_anchor + V_生成_乱数配列(15, p_seed: 2); // 35bp(短い)
            var l_unitigC = l_anchor + V_生成_乱数配列(2000, p_seed: 3); // 2020bp(長い)

            List<string> l_unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> l_kmerDict = [];

            var l_id = 1;
            foreach (var l_seq in new[] { l_unitigA, l_unitigB, l_unitigC })
            {
                l_unitigList.Add(l_seq);
                l_unitigList.Add(Util.V_逆相補(l_seq));
                for (var i = k長; i <= l_seq.Length; i++)
                {
                    var l_startPos = i - k長;
                    var l_key = new KmerKey(l_seq.AsSpan(l_startPos, k長));
                    V_登録_kmer(l_kmerDict, l_key, l_id, l_startPos);
                    V_登録_kmer(l_kmerDict, l_key.Get_逆相補(), -l_id, l_seq.Length - i);
                }
                l_id++;
            }

            var l_graph = UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k長, 曖昧kmer番号);
            return (l_unitigList, l_graph, 1, 2, 3);
        }

        /// <summary>
        /// どこも結合していない状態の結合表を作る
        /// </summary>
        /// <param name="p_graph">対象のユニティググラフ</param>
        /// <returns>結合表</returns>
        private static int[] V_構築_未結合表(UnitigGraph p_graph)
        {
            var l_merge = new int[p_graph.A_出辺.Count];
            Array.Fill(l_merge, -1);
            return l_merge;
        }

        /// <summary>
        /// 中央 150 のフラグメント長標本 (密度較正に使う)
        /// </summary>
        private static List<int> Get_同一ユニティグ標本(int p_件数 = 300) => [.. Enumerable.Repeat(150, p_件数)];

        /// <summary>
        /// 較正なしでは、生カウントの差だけでは優勢と判定されない
        /// </summary>
        [Fact]
        public void 較正なしでは生カウントの差だけでは優勢と判定されない()
        {
            var (l_unitigList, l_graph, l_aId, l_bId, l_cId) = V_構築();
            var l_a = ContigMaker.Get_頂点番号(l_aId);
            var l_b = ContigMaker.Get_頂点番号(l_bId);
            var l_c = ContigMaker.Get_頂点番号(l_cId);

            // 生カウントでは C(4) が B(3) よりわずかに多いが、
            // 優勢比 4/7=0.571 は閾値 0.8 を超えない
            Dictionary<(int, int), ulong> l_pairLink = new() { [(l_a, l_b)] = 3, [(l_a, l_c)] = 4 };
            Dictionary<int, int> l_copyNumber = new() { [l_aId] = 1, [l_bId] = 1, [l_cId] = 1 };

            var l_merge = V_構築_未結合表(l_graph);
            _ = BeamSearchExtender.V_延長_先読み(
                l_graph, l_unitigList, l_merge, l_pairLink, l_copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 3, p_較正器: null);

            Assert.Equal(-1, l_merge[l_a]);
        }

        /// <summary>
        /// 同じ観測本数 (3 vs 4) でも、C は B よりずっと長いぶん期待本数も
        /// 大きい
        /// </summary>
        /// <remarks>
        /// 期待本数との比を取ると B(短い) のほうが実際には強い証拠で
        /// あるとわかり、優勢閾値を超えて A→B が選ばれるはず
        /// </remarks>
        [Fact]
        public void 較正すると生カウントでは決められなかった短い側の分岐が選ばれる()
        {
            var (l_unitigList, l_graph, l_aId, l_bId, l_cId) = V_構築();
            var l_a = ContigMaker.Get_頂点番号(l_aId);
            var l_b = ContigMaker.Get_頂点番号(l_bId);
            var l_c = ContigMaker.Get_頂点番号(l_cId);

            Dictionary<(int, int), ulong> l_pairLink = new() { [(l_a, l_b)] = 3, [(l_a, l_c)] = 4 };
            Dictionary<int, int> l_copyNumber = new() { [l_aId] = 1, [l_bId] = 1, [l_cId] = 1 };

            var l_較正器 = 証拠較正器.Get_較正器(
                Get_同一ユニティグ標本(), p_リード長: 30,
                [(long)l_unitigList[l_a].Length, (long)l_unitigList[l_b].Length, (long)l_unitigList[l_c].Length]);
            Assert.True(l_較正器.A_使えるか);

            var l_merge = V_構築_未結合表(l_graph);
            var l_committed = BeamSearchExtender.V_延長_先読み(
                l_graph, l_unitigList, l_merge, l_pairLink, l_copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 3, p_較正器: l_較正器);

            Assert.True(l_committed > 0, "calibrated lookahead should have resolved the junction toward the short flank");
            Assert.Equal(l_b, l_merge[l_a]);
            Assert.Equal(l_a ^ 1, l_merge[l_b ^ 1]);
        }
    }
}
