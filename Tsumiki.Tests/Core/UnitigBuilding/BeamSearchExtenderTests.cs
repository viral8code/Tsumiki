using Tsumiki.Common;
using Tsumiki.Core.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 先読み (ビームサーチ) による分岐解決の検証
    /// </summary>
    /// <remarks>
    /// 相互一意性の判定は「その 1 歩だけ」を見るため、分岐の直後だけを見ると
    /// 五分五分に見えるが、2〜3 本先まで進めると片方だけがペアエンドの証拠と
    /// 整合する、という状況を取りこぼす<br/>
    /// ここでは
    /// A →(B or C)、B → D、C → E
    /// という形で、A の直後には証拠が無く D の位置に初めて証拠が現れる構成を作り、
    /// 先読みによって A → B が選ばれることを確認する
    /// </remarks>
    public class BeamSearchExtenderTests
    {
        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int 曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 8;

        // k=8 で 5 本すべてを通じて重複する正規化 k-mer が無いことを確認済みの構成

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
        /// B の先へ続く配列
        /// </summary>
        private const string ユニティグD = "CCGGATCAAAGCCACGGCTAG";

        /// <summary>
        /// C の先へ続く配列
        /// </summary>
        private const string ユニティグE = "TTCTACCAAAGGCTAGTATGA";

        /// <summary>
        /// A →(B or C)、B → D、C → E という形の分岐構造を持つグラフを作る
        /// </summary>
        /// <returns>ユニティグ一覧とグラフ</returns>
        private static (List<string> A_ユニティグ一覧, UnitigGraph A_グラフ) V_構築()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            List<string> l_unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> l_kmerDict = [];

            var l_id = 1;
            foreach (var l_seq in new[] { ユニティグA, ユニティグB, ユニティグC, ユニティグD, ユニティグE })
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
            return (l_unitigList, UnitigGraph.Get_グラフ(l_unitigList, l_kmerDict, k長, 曖昧kmer番号));
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
        /// 証拠が 1 歩先にしか現れない分岐を、先読みによって解決できる
        /// </summary>
        [Fact]
        public void 証拠が1歩先にしか現れない分岐を先読みで解決する()
        {
            var (l_unitigList, l_graph) = V_構築();
            var l_a = ContigMaker.Get_頂点番号(1);
            var l_b = ContigMaker.Get_頂点番号(2);
            var l_d = ContigMaker.Get_頂点番号(4);

            // 前提: A は B と C の両方へ伸びられる (1 歩だけでは決められない)
            Assert.Equal(2, l_graph.A_出辺[l_a].Count);

            // 証拠は A の直後 (B/C) ではなく、その次の D に現れる
            Dictionary<(int, int), ulong> l_pairLink = new() { [(l_a, l_d)] = 30 };
            Dictionary<int, int> l_copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var l_merge = V_構築_未結合表(l_graph);
            var l_committed = BeamSearchExtender.V_延長_先読み(
                l_graph, l_unitigList, l_merge, l_pairLink, l_copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.True(l_committed > 0, "lookahead should have resolved at least one junction");
            Assert.Equal(l_b, l_merge[l_a]);
            // 逆鎖側も対称に設定されていること
            Assert.Equal(l_a ^ 1, l_merge[l_b ^ 1]);
        }

        /// <summary>
        /// どちらの枝にも同程度の証拠がある場合は、僅差で選ばずに繋がない
        /// </summary>
        /// <remarks>
        /// ビームサーチの利点は「広く探して有力な仮説が一致する部分にだけ
        /// コミットする」ことにあり、五分五分の分岐で 1 本を選ぶことではない
        /// </remarks>
        [Fact]
        public void 両方の枝に同程度の証拠がある場合は繋がない()
        {
            var (l_unitigList, l_graph) = V_構築();
            var l_a = ContigMaker.Get_頂点番号(1);
            var l_d = ContigMaker.Get_頂点番号(4);
            var l_e = ContigMaker.Get_頂点番号(5);

            Dictionary<(int, int), ulong> l_pairLink = new() { [(l_a, l_d)] = 20, [(l_a, l_e)] = 19 };
            Dictionary<int, int> l_copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var l_merge = V_構築_未結合表(l_graph);
            _ = BeamSearchExtender.V_延長_先読み(
                l_graph, l_unitigList, l_merge, l_pairLink, l_copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.Equal(-1, l_merge[l_a]);
        }

        /// <summary>
        /// ペアエンドの証拠がまったく無ければ、根拠が無いので繋がない
        /// </summary>
        [Fact]
        public void ペア証拠が全く無い場合は繋がない()
        {
            var (l_unitigList, l_graph) = V_構築();
            var l_a = ContigMaker.Get_頂点番号(1);

            Dictionary<(int, int), ulong> l_pairLink = [];
            Dictionary<int, int> l_copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var l_merge = V_構築_未結合表(l_graph);
            _ = BeamSearchExtender.V_延長_先読み(
                l_graph, l_unitigList, l_merge, l_pairLink, l_copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.Equal(-1, l_merge[l_a]);
        }

        /// <summary>
        /// 証拠はあるが少なすぎる場合、偶然の一致で繋いでしまわないよう見送る
        /// </summary>
        [Fact]
        public void 証拠数が最小値未満の場合は繋がない()
        {
            var (l_unitigList, l_graph) = V_構築();
            var l_a = ContigMaker.Get_頂点番号(1);
            var l_d = ContigMaker.Get_頂点番号(4);

            Dictionary<(int, int), ulong> l_pairLink = new() { [(l_a, l_d)] = 2 };
            Dictionary<int, int> l_copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var l_merge = V_構築_未結合表(l_graph);
            _ = BeamSearchExtender.V_延長_先読み(
                l_graph, l_unitigList, l_merge, l_pairLink, l_copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 10);

            Assert.Equal(-1, l_merge[l_a]);
        }

        /// <summary>
        /// いま反復配列 (多コピー) の上にいて、単一コピーの足場が 1 つも取れない
        /// 場合は、どのコピーにいるのか分からないので進む方向を選べない
        /// </summary>
        /// <remarks>
        /// 反復の内部から読まれたリードはどのコピー由来か区別できない<br/>
        /// それが
        /// 反復が解けない理由そのものなので、そこを起点にしたペアの証拠は
        /// どの行き先にも付いてしまう<br/>
        /// 標本数が少ないと偶然の偏りが閾値を超えて
        /// 誤った側が選ばれる<br/>
        /// これは実際に起きた: 反復入りの合成ゲノム (A-R-B-R-C、R は150 bp の
        /// 2 コピー反復) で、R 自身を足場にしたために A-R-C という中間の B を
        /// 飛ばした contig が出力されていた (真値照合で発覚)
        /// </remarks>
        [Fact]
        public void 単一コピーの足場が無い反復上では繋がない()
        {
            var (l_unitigList, l_graph) = V_構築();
            var l_a = ContigMaker.Get_頂点番号(1);
            var l_d = ContigMaker.Get_頂点番号(4);

            // A 自身が 2 コピーの反復
            // 足場に使える単一コピーの unitig が無い
            Dictionary<int, int> l_copyNumber = new() { [1] = 2, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };
            // 片側にだけ強い (しかし信用してはいけない) 証拠を置く
            Dictionary<(int, int), ulong> l_pairLink = new() { [(l_a, l_d)] = 30 };

            var l_merge = V_構築_未結合表(l_graph);
            var l_committed = BeamSearchExtender.V_延長_先読み(
                l_graph, l_unitigList, l_merge, l_pairLink, l_copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.Equal(0, l_committed);
            Assert.Equal(-1, l_merge[l_a]);
        }

        /// <summary>
        /// 既に別の結合が入っている行き先へは、それを壊してまで繋がない
        /// (相互一意性を保つ)
        /// </summary>
        [Fact]
        public void 既に結合済みの行き先を奪って繋がない()
        {
            var (l_unitigList, l_graph) = V_構築();
            var l_a = ContigMaker.Get_頂点番号(1);
            var l_b = ContigMaker.Get_頂点番号(2);
            var l_d = ContigMaker.Get_頂点番号(4);

            Dictionary<(int, int), ulong> l_pairLink = new() { [(l_a, l_d)] = 30 };
            Dictionary<int, int> l_copyNumber = new() { [1] = 1, [2] = 1, [3] = 1, [4] = 1, [5] = 1 };

            var l_merge = V_構築_未結合表(l_graph);
            // B には既に (別の経路からの) 結合が入っていることにする
            l_merge[l_b ^ 1] = ContigMaker.Get_頂点番号(5) ^ 1;

            _ = BeamSearchExtender.V_延長_先読み(
                l_graph, l_unitigList, l_merge, l_pairLink, l_copyNumber,
                p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.Equal(-1, l_merge[l_a]);
        }
    }
}
