using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 先読み (ビームサーチ) による分岐解決の検証
    /// </summary>
    public class BeamSearchExtenderTests
    {
        #region 定数

        /// <summary>
        /// 項目 lookahead should have re ast one junction
        /// </summary>
        private const string C_項目_lookahead_should_have_re_ast_one_junction = "lookahead should have resolved at least one junction";

        /// <summary>
        /// 項目 区間の上限により B が候補に含まれ 証拠に基づいて選ばれるはず
        /// </summary>
        private const string C_項目_区間の上限により_B_が候補に含まれ_証拠に基づいて選ばれるはず = "区間の上限により B が候補に含まれ、証拠に基づいて選ばれるはず";

        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int C_曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k長 = 8;

        /// <summary>
        /// 分岐元
        /// </summary>
        private const string C_unitigA = "TGGCAAGTCACTCTCGACCGA";

        /// <summary>
        /// 分岐先の片方
        /// </summary>
        private const string C_unitigB = "CGACCGAACGGCGCCGGATC";

        /// <summary>
        /// 分岐先のもう片方
        /// </summary>
        private const string C_unitigC = "CGACCGACTGTAATTCTACC";

        /// <summary>
        /// B の先へ続く配列
        /// </summary>
        private const string C_unitigD = "CCGGATCAAAGCCACGGCTAG";

        /// <summary>
        /// C の先へ続く配列
        /// </summary>
        private const string C_unitigE = "TTCTACCAAAGGCTAGTATGA";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 証拠が 1 歩先にしか現れない分岐を、先読みによって解決できる
        /// </summary>
        [Fact]
        public void V_証拠が1歩先にしか現れない分岐を先読みで解決する()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_終端配列 = ContigMaker.Get_頂点番号(4);
            Assert.Equal(2, l_グラフ.A_出辺[l_先頭配列].Count);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 30UL
            };
            Dictionary<int, int> l_copyNumber = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            var l_merge = V_構築_未結合表(l_グラフ);
            var l_committed = BeamSearchExtender.V_延長_先読み(l_グラフ, l_unitig一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.True(l_committed > 0, C_項目_lookahead_should_have_re_ast_one_junction);
            Assert.Equal(l_中間配列, l_merge[l_先頭配列]);
            Assert.Equal(l_先頭配列 ^ 1, l_merge[l_中間配列 ^ 1]);
        }

        /// <summary>
        /// どちらの枝にも同程度の証拠がある場合は、僅差で選ばずに繋がない
        /// </summary>
        [Fact]
        public void V_両方の枝に同程度の証拠がある場合は繋がない()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_終端配列 = ContigMaker.Get_頂点番号(4);
            var l_e = ContigMaker.Get_頂点番号(5);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 20UL,
                [(l_先頭配列, l_e)] = 19UL
            };
            Dictionary<int, int> l_copyNumber = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            var l_merge = V_構築_未結合表(l_グラフ);
            _ = BeamSearchExtender.V_延長_先読み(l_グラフ, l_unitig一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(-1, l_merge[l_先頭配列]);
        }

        /// <summary>
        /// ペアエンドの証拠がまったく無ければ、根拠が無いので繋がない
        /// </summary>
        [Fact]
        public void V_ペア証拠が全く無い場合は繋がない()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            Dictionary<(int, int), ulong> l_ペア連結 = [];
            Dictionary<int, int> l_copyNumber = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            var l_merge = V_構築_未結合表(l_グラフ);
            _ = BeamSearchExtender.V_延長_先読み(l_グラフ, l_unitig一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(-1, l_merge[l_先頭配列]);
        }

        /// <summary>
        /// 証拠はあるが少なすぎる場合、偶然の一致で繋いでしまわないよう見送る
        /// </summary>
        [Fact]
        public void V_証拠数が最小値未満の場合は繋がない()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_終端配列 = ContigMaker.Get_頂点番号(4);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 2UL
            };
            Dictionary<int, int> l_copyNumber = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            var l_merge = V_構築_未結合表(l_グラフ);
            _ = BeamSearchExtender.V_延長_先読み(l_グラフ, l_unitig一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 10UL);
            Assert.Equal(-1, l_merge[l_先頭配列]);
        }

        /// <summary>
        /// いま反復配列 (多コピー) の上にいて、単一コピーの足場が 1 つも取れない場合は、どのコピーにいるのか分からないので進む方向を選べない
        /// </summary>
        [Fact]
        public void V_単一コピーの足場が無い反復上では繋がない()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_終端配列 = ContigMaker.Get_頂点番号(4);
            Dictionary<int, int> l_copyNumber = new()
            {
                [1] = 2,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 30UL
            };
            var l_merge = V_構築_未結合表(l_グラフ);
            var l_committed = BeamSearchExtender.V_延長_先読み(l_グラフ, l_unitig一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(0, l_committed);
            Assert.Equal(-1, l_merge[l_先頭配列]);
        }

        /// <summary>
        /// 点推定ではコピー数 0 (=候補から除外) の行き先でも、コピー数区間の上限が正なら候補に含め、 実際に強い証拠があれば選べることを確かめる
        /// </summary>
        [Fact]
        public void V_点推定がゼロでも区間の上限が正なら候補に含める()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_終端配列 = ContigMaker.Get_頂点番号(4);
            Dictionary<int, int> l_copyNumber = new()
            {
                [1] = 1,
                [2] = 0,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            Dictionary<int, コピー数区間> l_区間 = new()
            {
                [2] = new コピー数区間(0, 1)
            };
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 30UL
            };
            var l_merge = V_構築_未結合表(l_グラフ);
            var l_committed = BeamSearchExtender.V_延長_先読み(l_グラフ, l_unitig一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5UL, p_コピー数区間: l_区間);
            Assert.True(l_committed > 0, C_項目_区間の上限により_B_が候補に含まれ_証拠に基づいて選ばれるはず);
            Assert.Equal(l_中間配列, l_merge[l_先頭配列]);
        }

        /// <summary>
        /// コピー数区間を渡さない場合、点推定がそのまま予算になり、点推定ゼロの行き先は候補に含めない
        /// </summary>
        [Fact]
        public void V_コピー数区間を渡さない場合は点推定がそのまま予算になる()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_終端配列 = ContigMaker.Get_頂点番号(4);
            Dictionary<int, int> l_copyNumber = new()
            {
                [1] = 1,
                [2] = 0,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 30UL
            };
            var l_merge = V_構築_未結合表(l_グラフ);
            var l_committed = BeamSearchExtender.V_延長_先読み(l_グラフ, l_unitig一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(0, l_committed);
        }

        /// <summary>
        /// 既に別の結合が入っている行き先へは、それを壊してまで繋がない (相互一意性を保つ)
        /// </summary>
        [Fact]
        public void V_既に結合済みの行き先を奪って繋がない()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_終端配列 = ContigMaker.Get_頂点番号(4);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 30UL
            };
            Dictionary<int, int> l_copyNumber = new()
            {
                [1] = 1,
                [2] = 1,
                [3] = 1,
                [4] = 1,
                [5] = 1
            };
            var l_merge = V_構築_未結合表(l_グラフ);
            l_merge[l_中間配列 ^ 1] = ContigMaker.Get_頂点番号(5) ^ 1;
            _ = BeamSearchExtender.V_延長_先読み(l_グラフ, l_unitig一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(-1, l_merge[l_先頭配列]);
        }

        /// <summary>
        /// 出次数 1 の始点から入次数 1 の終点への辺だけを、コピー数に依らず結合できる辺とみなす
        /// </summary>
        [Fact]
        public void Is構造上一意な辺_分岐を持たない辺だけを一意とみなす()
        {
            var (_, l_グラフ) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_終端配列 = ContigMaker.Get_頂点番号(4);
            Assert.True(l_グラフ.Is構造上一意な辺(l_中間配列, l_終端配列));
            Assert.False(l_グラフ.Is構造上一意な辺(l_先頭配列, l_中間配列));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// A → (B or C) 、B → D、C → E という形の分岐構造を持つグラフを作る
        /// </summary>
        /// <returns>unitig 一覧とグラフ</returns>
        private static (List<string> A_unitig一覧, UnitigGraph A_グラフ) V_構築()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
                A_スレッド数 = 1
            };
            List<string> l_unitig一覧 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_unitigID, int A_位置)> l_kmer辞書 = [];
            var l_ID = 1;
            foreach (var l_配列 in new[]
            {
                C_unitigA,
                C_unitigB,
                C_unitigC,
                C_unitigD,
                C_unitigE
            }

            )
            {
                l_unitig一覧.Add(l_配列);
                l_unitig一覧.Add(Util.V_逆相補(l_配列));
                for (var i = C_k長; i <= l_配列.Length; i++)
                {
                    var l_開始位置 = i - C_k長;
                    var l_キー = new KmerKey(l_配列.AsSpan(l_開始位置, C_k長));
                    V_登録_kmer(l_kmer辞書, l_キー, l_ID, l_開始位置);
                    V_登録_kmer(l_kmer辞書, l_キー.Get_逆相補(), -l_ID, l_配列.Length - i);
                }

                l_ID++;
            }

            return (l_unitig一覧, UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号));
        }

        /// <summary>
        /// k-mer を、それが載る unitig と開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_キー">登録する k-mer</param>
        /// <param name="p_ID">unitig ID</param>
        /// <param name="p_位置">unitig 内の開始位置</param>
        private static void V_登録_kmer(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存値))
            {
                if (l_既存値.Item1 is C_曖昧kmer番号 || l_既存値.Item1 == p_ID)
                {
                    return;
                }

                p_辞書[p_キー] = (C_曖昧kmer番号, 0);
                return;
            }

            p_辞書[p_キー] = (p_ID, p_位置);
        }

        /// <summary>
        /// どこも結合していない状態の結合表を作る
        /// </summary>
        /// <param name="p_グラフ">対象の unitig グラフ</param>
        /// <returns>結合表</returns>
        private static int[] V_構築_未結合表(UnitigGraph p_グラフ)
        {
            var l_merge = new int[p_グラフ.A_出辺.Count];
            Array.Fill(l_merge, -1);
            return l_merge;
        }

        #endregion
    }
}
