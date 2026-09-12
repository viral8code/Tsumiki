using Tsumiki.Commons;
using Tsumiki.Cores.Evidence;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 提案 D: BeamSearchExtender のスコアを生カウントではなく期待本数との比で測ることの検証
    /// </summary>
    /// <remarks>
    /// 分岐元 A から、短い unitig B と長い unitig C へ分岐する構成を作る<br/>
    /// C は B よりずっと長いため、同じフラグメント長分布のもとでは「両端が収まる開始位置」の窓が B よりずっと広い (=同じ観測本数でも期待本数は C のほうが大きい) <br/>
    /// 生カウントでは C がわずかに優勢に見えても優勢閾値を超えないケースで、期待本数との比を取ると実際には B のほうが強い証拠であるとして正しく選ばれることを確認する
    /// </remarks>
    public class BeamSearchExtenderCalibrationTests
    {
        #region 定数

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

        #endregion

        #region 公開メソッド

        /// <summary>
        /// フラグメントが届かない先読み経路を支持として加算しない
        /// </summary>
        [Fact]
        public void V_先読み距離_到達不能な証拠を除外()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長 };
            var l_較正器 = 証拠較正器.Get_較正器(Enumerable.Repeat(400, 100).ToArray(), 100, [2000L]);
            List<string> l_配列 = [new string('A', 500), new string('C', 500)];
            Dictionary<(int, int), ulong> l_証拠 = new() { [(0, 1)] = 10UL };
            var l_近い = BeamSearchExtender.Get_スコア([(0, 0)], 1, l_配列, l_証拠, l_較正器, 0);
            var l_遠い = BeamSearchExtender.Get_スコア([(0, 0)], 1, l_配列, l_証拠, l_較正器, 1000);
            Assert.True(l_近い.A_正規化 > 0D);
            Assert.Equal(0D, l_遠い.A_正規化);
        }

        /// <summary>
        /// 較正なしでは、生カウントの差だけでは優勢と判定されない
        /// </summary>
        [Fact]
        public void V_較正なしでは生カウントの差だけでは優勢と判定されない()
        {
            var (l_ユニティグ一覧, l_グラフ, l_aId, l_bId, l_cId) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(l_aId);
            var l_中間配列 = ContigMaker.Get_頂点番号(l_bId);
            var l_末尾配列 = ContigMaker.Get_頂点番号(l_cId);

            // 生カウントでは C (4) が B (3) よりわずかに多いが、
            // 優勢比 4/7=0.571 は閾値 0.8 を超えない
            Dictionary<(int, int), ulong> l_ペア連結 = new() { [(l_先頭配列, l_中間配列)] = 3UL, [(l_先頭配列, l_末尾配列)] = 4UL };
            Dictionary<int, int> l_copyNumber = new() { [l_aId] = 1, [l_bId] = 1, [l_cId] = 1 };

            var l_merge = V_構築_未結合表(l_グラフ);
            _ = BeamSearchExtender.V_延長_先読み(l_グラフ, l_ユニティグ一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 3UL, p_較正器: null);

            Assert.Equal(-1, l_merge[l_先頭配列]);
        }

        /// <summary>
        /// 同じ観測本数 (3 vs 4) でも、C は B よりずっと長いぶん期待本数も大きい
        /// </summary>
        /// <remarks>
        /// 期待本数との比を取ると B (短い) のほうが実際には強い証拠であるとわかり、優勢閾値を超えて A→B が選ばれるはず
        /// </remarks>
        [Fact]
        public void V_較正すると生カウントでは決められなかった短い側の分岐が選ばれる()
        {
            var (l_ユニティグ一覧, l_グラフ, l_aId, l_bId, l_cId) = V_構築();
            var l_先頭配列 = ContigMaker.Get_頂点番号(l_aId);
            var l_中間配列 = ContigMaker.Get_頂点番号(l_bId);
            var l_末尾配列 = ContigMaker.Get_頂点番号(l_cId);

            Dictionary<(int, int), ulong> l_ペア連結 = new() { [(l_先頭配列, l_中間配列)] = 3UL, [(l_先頭配列, l_末尾配列)] = 4UL };
            Dictionary<int, int> l_copyNumber = new() { [l_aId] = 1, [l_bId] = 1, [l_cId] = 1 };

            var l_較正器 = 証拠較正器.Get_較正器(Get_同一ユニティグ標本(), p_リード長: 30, [(long)l_ユニティグ一覧[l_先頭配列].Length, (long)l_ユニティグ一覧[l_中間配列].Length, (long)l_ユニティグ一覧[l_末尾配列].Length]);
            Assert.True(l_較正器.A_使えるか);

            var l_merge = V_構築_未結合表(l_グラフ);
            var l_committed = BeamSearchExtender.V_延長_先読み(l_グラフ, l_ユニティグ一覧, l_merge, l_ペア連結, l_copyNumber, p_インサートサイズ: 400, p_優勢閾値: 0.8M, p_最小証拠数: 3UL, p_較正器: l_較正器);

            Assert.True(l_committed > 0, "calibrated lookahead should have resolved the junction toward the short flank");
            Assert.Equal(l_中間配列, l_merge[l_先頭配列]);
            Assert.Equal(l_先頭配列 ^ 1, l_merge[l_中間配列 ^ 1]);
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
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_キー">登録する k-mer</param>
        /// <param name="p_ID">ユニティグ ID</param>
        /// <param name="p_位置">ユニティグ内の開始位置</param>
        private static void V_登録_kmer(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存値))
            {
                if (l_既存値.Item1 is 曖昧kmer番号 || l_既存値.Item1 == p_ID)
                {
                    return;
                }
                p_辞書[p_キー] = (曖昧kmer番号, 0);
                return;
            }
            p_辞書[p_キー] = (p_ID, p_位置);
        }

        /// <summary>
        /// A (250 bp) が分岐元、B (35 bp, 短い) と C (2020 bp, 長い) がその行き先
        /// </summary>
        /// <remarks>
        /// A の末尾 20 塩基 (=k-1) を B ・ C 両方の先頭が共有することで分岐にする
        /// </remarks>
        /// <returns></returns>
        private static (List<string> A_ユニティグ一覧, UnitigGraph A_グラフ, int A_分岐元ID, int A_短い方ID, int A_長い方ID) V_構築()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };

            var l_anchor = V_生成_乱数配列(アンカー長, p_乱数種: 20260908);
            var l_unitigA = V_生成_乱数配列(230, p_乱数種: 1) + l_anchor; // 250 bp
            var l_unitigB = l_anchor + V_生成_乱数配列(15, p_乱数種: 2); // 35 bp (短い)
            var l_unitigC = l_anchor + V_生成_乱数配列(2000, p_乱数種: 3); // 2020 bp (長い)

            List<string> l_ユニティグ一覧 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_ユニティグID, int A_位置)> l_kmer辞書 = [];

            var l_ID = 1;
            foreach (var l_配列 in new[] { l_unitigA, l_unitigB, l_unitigC })
            {
                l_ユニティグ一覧.Add(l_配列);
                l_ユニティグ一覧.Add(Util.V_逆相補(l_配列));
                for (var i = k長; i <= l_配列.Length; i++)
                {
                    var l_開始位置 = i - k長;
                    var l_キー = new KmerKey(l_配列.AsSpan(l_開始位置, k長));
                    V_登録_kmer(l_kmer辞書, l_キー, l_ID, l_開始位置);
                    V_登録_kmer(l_kmer辞書, l_キー.Get_逆相補(), -l_ID, l_配列.Length - i);
                }
                l_ID++;
            }

            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);
            return (l_ユニティグ一覧, l_グラフ, 1, 2, 3);
        }

        /// <summary>
        /// どこも結合していない状態の結合表を作る
        /// </summary>
        /// <param name="p_グラフ">対象のユニティググラフ</param>
        /// <returns>結合表</returns>
        private static int[] V_構築_未結合表(UnitigGraph p_グラフ)
        {
            var l_merge = new int[p_グラフ.A_出辺.Count];
            Array.Fill(l_merge, -1);
            return l_merge;
        }

        /// <summary>
        /// 中央 150 のフラグメント長標本 (密度較正に使う)
        /// </summary>
        /// <param name="p_件数"></param>
        /// <returns></returns>
        private static List<int> Get_同一ユニティグ標本(int p_件数 = 300) => [.. Enumerable.Repeat(150, p_件数)];

        #endregion

    }
}
