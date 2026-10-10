using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core.Evaluation
{
    /// <summary>
    /// 継ぎ目の候補の拾い方と、リードの証拠の数え方を固定する
    /// </summary>
    public class JunctionRiskEvaluatorTests
    {
        #region 定数

        /// <summary>
        /// 塩基配列 NNNNN
        /// </summary>
        private const string C_塩基配列_NNNNN = "NNNNN";

        /// <summary>
        /// 項目 s
        /// </summary>
        private const string C_項目_s = "s";

        /// <summary>
        /// 項目 S1 circular
        /// </summary>
        private const string C_項目_S1_circular = "S1_circular";

        /// <summary>
        /// 項目 S1 1
        /// </summary>
        private const string C_項目_S1_1 = "S1_1";

        /// <summary>
        /// 項目 S1 2
        /// </summary>
        private const string C_項目_S1_2 = "S1_2";

        /// <summary>
        /// 項目 S2 1
        /// </summary>
        private const string C_項目_S2_1 = "S2_1";

        /// <summary>
        /// 項目 S2 2
        /// </summary>
        private const string C_項目_S2_2 = "S2_2";

        /// <summary>
        /// 項目 S2
        /// </summary>
        private const string C_項目_S2 = "S2";

        /// <summary>
        /// 手で作る断片の長さ
        /// </summary>
        private const int C_断片長 = 800;

        /// <summary>
        /// 手で作る深さ
        /// </summary>
        private const int C_深さ = 30;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// テストに必要な共有状態を初期化する
        /// </summary>
        public JunctionRiskEvaluatorTests()
        {
            ConfigurationManager.A_実行時引数 = new Parameters();
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列の中で 2 回現れる区間と N に印を付け、固有の区間には付けないことを確かめる
        /// </summary>
        [Fact]
        public void Get_反復の印_2回現れる区間とNに印を付ける()
        {
            var l_反復 = Get_乱数配列(100, 1);
            var l_配列 = Get_乱数配列(300, 2) + l_反復 + Get_乱数配列(300, 3) + l_反復 + C_塩基配列_NNNNN + Get_乱数配列(300, 4);
            var l_印 = JunctionRiskEvaluator.Get_反復の印([l_配列])[0];
            Assert.All(Enumerable.Range(300, 100), p => Assert.True(l_印[p]));
            Assert.All(Enumerable.Range(700, 100), p => Assert.True(l_印[p]));
            Assert.All(Enumerable.Range(800, 5), p => Assert.True(l_印[p]));
            Assert.All(Enumerable.Range(0, 250), p => Assert.False(l_印[p]));
            Assert.All(Enumerable.Range(450, 200), p => Assert.False(l_印[p]));
        }

        /// <summary>
        /// 反復の区間ごとに、区間を覆う断片の数と、区間の手前で相方が組めていない錨の数を数えることを確かめる (区間の端は、両コピーの外側の塩基が偶然一致すると数塩基広がる)
        /// </summary>
        [Fact]
        public void Get_候補群_跨ぐ組と外れ錨を区間ごとに数える()
        {
            var l_反復 = Get_乱数配列(200, 11);
            var l_配列 = Get_乱数配列(1_000, 12) + l_反復 + Get_乱数配列(2_000, 13) + l_反復 + Get_乱数配列(1_000, 14);
            var l_断片 = Enumerable.Range(0, 10).Select(_ => (700, 700 + C_断片長)).Concat(Enumerable.Range(0, 200).Select(i => (1_300 + (i * 5), 1_300 + (i * 5) + C_断片長))).Order().ToArray();
            var l_外れ錨 = Enumerable.Range(0, 7).Select(i => (3_000 + (i * 20), 1)).ToArray();
            var l_候補群 = JunctionRiskEvaluator.Get_候補群(Get_証拠(l_配列, l_断片, l_外れ錨), [C_項目_s]);
            Assert.Equal(2, l_候補群.Count);
            var l_跨がれた = l_候補群.Single(x => x.A_開始 < 1_500);
            var l_外れた = l_候補群.Single(x => x.A_開始 > 1_500);
            Assert.InRange(l_跨がれた.A_開始, 995, 1_000);
            Assert.InRange(l_跨がれた.A_終了, 1_200, 1_205);
            Assert.Equal(10, l_跨がれた.A_跨ぐ組);
            Assert.Equal(0, l_跨がれた.A_外れ錨);
            Assert.InRange(l_外れた.A_開始, 3_195, 3_200);
            Assert.InRange(l_外れた.A_終了, 3_400, 3_405);
            Assert.Equal(0, l_外れた.A_跨ぐ組);
            Assert.Equal(7, l_外れた.A_外れ錨);
            Assert.All(l_候補群, x => Assert.True(x.A_期待の組 > 0D));
            Assert.All(l_候補群, x => Assert.Equal(1D, x.A_コピー数, 3));
        }

        /// <summary>
        /// 跨ぐ組が無く外れ錨の多い継ぎ目に、跨ぐ組の揃った継ぎ目より高い誤りの確率を付けることを確かめる
        /// </summary>
        [Fact]
        public void Get_評価_跨ぐ組が無く外れ錨が多いほど確率が高い()
        {
            var l_支えあり = new 継ぎ目候補(C_項目_s, 1_000, 1_200, false, 1D, 0, 20, 20D, 0, 0D, 1D, 1D);
            var l_支えなし = l_支えあり with
            {
                A_跨ぐ組 = 0,
                A_外れ錨 = 40
            };
            var l_低い = JunctionRiskModel.Get_評価(l_支えあり);
            var l_高い = JunctionRiskModel.Get_評価(l_支えなし);
            Assert.Equal(継ぎ目の組.組で跨げる, l_低い.A_組);
            Assert.True(l_低い.A_誤りの確率 < 0.01D);
            Assert.True(l_高い.A_誤りの確率 > 0.5D);
        }

        /// <summary>
        /// N のギャップ、組で跨げる継ぎ目、組で跨げない継ぎ目を分けることを確かめる
        /// </summary>
        [Fact]
        public void Get_組_ギャップと期待の組で分ける()
        {
            var l_候補 = new 継ぎ目候補(C_項目_s, 1_000, 1_200, false, 1D, 0, 0, 5D, 0, 0D, 1D, 1D);
            Assert.Equal(継ぎ目の組.組で跨げる, JunctionRiskModel.Get_組(l_候補));
            Assert.Equal(継ぎ目の組.組で跨げない, JunctionRiskModel.Get_組(l_候補 with { A_期待の組 = 0.5D }));
            Assert.Equal(継ぎ目の組.ギャップ, JunctionRiskModel.Get_組(l_候補 with { A_Isギャップ = true }));
        }

        /// <summary>
        /// 反復の区間で切るときは区間を両側の片に残し、切った環状配列からは環状の目印を外すことを確かめる
        /// </summary>
        [Fact]
        public void Get_切った配列群_反復の区間を両側に残して切る()
        {
            var l_左 = Get_乱数配列(300, 31);
            var l_反復 = Get_乱数配列(100, 32);
            var l_右 = Get_乱数配列(300, 33);
            var l_片群 = JunctionRiskEvaluator.Get_切った配列群(C_項目_S1_circular, l_左 + l_反復 + l_右, [(300, 400)]);
            Assert.Equal([(C_項目_S1_1, l_左 + l_反復), (C_項目_S1_2, l_反復 + l_右)], l_片群);
        }

        /// <summary>
        /// 区間に N のギャップがあれば N の連続の手前と後ろで切り、切る所が無ければそのまま返すことを確かめる
        /// </summary>
        [Fact]
        public void Get_切った配列群_ギャップはNの手前と後ろで切る()
        {
            var l_左 = Get_乱数配列(300, 41);
            var l_右 = Get_乱数配列(300, 42);
            var l_配列 = l_左 + new string('N', 100) + l_右;
            Assert.Equal([(C_項目_S2_1, l_左), (C_項目_S2_2, l_右)], JunctionRiskEvaluator.Get_切った配列群(C_項目_S2, l_配列, [(270, 430)]));
            Assert.Equal([(C_項目_S2, l_配列)], JunctionRiskEvaluator.Get_切った配列群(C_項目_S2, l_配列, []));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 一様な深さ (反復の中は一意に当たらない) と、与えた断片・外れ錨を持つ証拠を作る
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_断片">始まりの順の断片</param>
        /// <param name="p_外れ錨">位置の順の外れ錨</param>
        /// <returns></returns>
        private static 継ぎ目の証拠 Get_証拠(string p_配列, (int, int)[] p_断片, (int, int)[] p_外れ錨)
        {
            var l_印 = JunctionRiskEvaluator.Get_反復の印([p_配列])[0];
            return new 継ぎ目の証拠
            {
                A_配列群 = [p_配列],
                A_深さ = [Enumerable.Repeat(C_深さ, p_配列.Length).ToArray()],
                A_一意な深さ = [l_印.Select(x => x ? 0 : C_深さ).ToArray()],
                A_左の切れ端 = [new int[p_配列.Length + 1]],
                A_右の切れ端 = [new int[p_配列.Length + 1]],
                A_読み = [[]],
                A_断片 = [p_断片],
                A_外れ錨 = [p_外れ錨],
            };
        }

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ"></param>
        /// <param name="p_種"></param>
        /// <returns></returns>
        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            const string l_塩基 = Consts.塩基文字;
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => l_塩基[l_乱数.Next(4)]));
        }

        #endregion
    }
}
