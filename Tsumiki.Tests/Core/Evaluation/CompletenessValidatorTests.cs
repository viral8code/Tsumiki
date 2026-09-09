using Tsumiki.Core.Evaluation;
using Tsumiki.Core;
using Tsumiki.Model.Evaluation;
using Tsumiki.Model.Polishing;
using Tsumiki.Model.Reporting;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 完全長の判定を固定する。
    ///
    /// 要点は「材料が無いことを合格にしない」ことと「不合格と判定不能を
    /// 混同しない」こと。どちらを崩しても、根拠の無い完全長が通ってしまう。
    /// </summary>
    public class CompletenessValidatorTests
    {
        /// <summary>リードに裏付けの無い位置が一つも無い検査結果。</summary>
        private static 支持検査結果 Get_良好な支持()
        {
            return new 支持検査結果(A_r長: 31, A_調べた位置数: 100000, A_支持のない位置数: 0, A_区間: []);
        }

        private static 整合性検査結果 Get_良好な自己検査()
        {
            // 取りこぼし 1%、出しすぎ 0%。
            return new 整合性検査結果(1000, 1000, 1000, 10, 0, 0);
        }

        private static ポリッシュ統計 Get_良好な深度()
        {
            return new ポリッシュ統計(1, 1000, 1000, 0, 0, 0, 1000, 80);
        }

        private static IReadOnlyList<環状閉鎖検証結果> Get_裏付けのある閉鎖()
        {
            return [new 環状閉鎖検証結果("scaffold1_circular", 1000, 12, 5)];
        }

        private static 検査判定 Get_判定(完全性判定結果 p_判定, string p_キー)
        {
            return p_判定.A_検査項目.Single(x => x.A_キー == p_キー).A_判定;
        }

        [Fact]
        public void Get_判定結果_全ての検査を通れば完全長になる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: Get_良好な自己検査(),
                p_閉鎖検証: Get_裏付けのある閉鎖(),
                p_ポリッシュ: Get_良好な深度(),
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            Assert.True(l_判定.A_完全長か);
            Assert.Equal(品質保証レベル.完全長, l_判定.A_品質保証レベル);
            Assert.Empty(l_判定.A_未達理由);
        }

        [Fact]
        public void Get_判定結果_閉じ目を検証していなければ完全長にはしない()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: Get_良好な自己検査(),
                p_閉鎖検証: null,
                p_ポリッシュ: Get_良好な深度(),
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            Assert.False(l_判定.A_完全長か);
            Assert.Equal(品質保証レベル.接合点が支持済み, l_判定.A_品質保証レベル);
            Assert.Equal(検査判定.判定不能, Get_判定(l_判定, "circular_closure"));
            Assert.Contains(未達理由.環状閉鎖を検証していない, l_判定.A_未達理由);
        }

        [Fact]
        public void Get_判定結果_閉じ目に裏付けが無い場合は不合格として区別する()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: Get_良好な自己検査(),
                p_閉鎖検証: [new 環状閉鎖検証結果("scaffold1_circular", 1000, 1, 5)],
                p_ポリッシュ: Get_良好な深度(),
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            Assert.Equal(検査判定.不合格, Get_判定(l_判定, "circular_closure"));
            Assert.Contains(未達理由.閉じ目がリードで裏付けられない, l_判定.A_未達理由);
        }

        [Fact]
        public void Get_判定結果_環状の配列が1本も無ければ不合格になる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: Get_良好な自己検査(),
                p_閉鎖検証: [],
                p_ポリッシュ: Get_良好な深度(),
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            Assert.Equal(検査判定.不合格, Get_判定(l_判定, "circular_closure"));
            Assert.Contains(未達理由.環状に閉じていない, l_判定.A_未達理由);
        }

        [Fact]
        public void Get_判定結果_未解決のギャップが残っていればペア整合で止まる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 3,
                p_整合性: Get_良好な自己検査(),
                p_閉鎖検証: Get_裏付けのある閉鎖(),
                p_ポリッシュ: Get_良好な深度(),
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.マッピング整合, l_判定.A_品質保証レベル);
            Assert.Contains(未達理由.未解決のギャップが残る, l_判定.A_未達理由);
        }

        [Fact]
        public void Get_判定結果_決めきれない分岐が残っていれば接合点の支持で止まる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: Get_良好な自己検査(),
                p_閉鎖検証: Get_裏付けのある閉鎖(),
                p_ポリッシュ: Get_良好な深度(),
                p_曖昧箇所: [new 曖昧箇所(63, 曖昧箇所の種別.僅差, "unitig7+", 1.2, 1.1, 9, 0.95)],
                p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.ペア整合, l_判定.A_品質保証レベル);
            Assert.Equal(検査判定.不合格, Get_判定(l_判定, "junction_support"));
            Assert.Equal(検査判定.不合格, Get_判定(l_判定, "no_alternative_path"));
        }

        [Fact]
        public void Get_判定結果_深度を測っていなければグラフ整合で止まる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: Get_良好な自己検査(),
                p_閉鎖検証: Get_裏付けのある閉鎖(),
                p_ポリッシュ: null,
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.グラフ整合, l_判定.A_品質保証レベル);
            Assert.Equal(検査判定.判定不能, Get_判定(l_判定, "coverage_continuity"));
            Assert.Contains(未達理由.深度を測っていない, l_判定.A_未達理由);
        }

        [Fact]
        public void Get_判定結果_自己検査ができていなければ出力のみに留まる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: null,
                p_閉鎖検証: Get_裏付けのある閉鎖(),
                p_ポリッシュ: Get_良好な深度(),
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.出力のみ, l_判定.A_品質保証レベル);
            Assert.Equal(検査判定.判定不能, Get_判定(l_判定, "graph_coverage"));
            Assert.Contains(未達理由.自己検査を行えなかった, l_判定.A_未達理由);
        }

        [Fact]
        public void Get_判定結果_取りこぼしが多ければグラフ被覆で落ちる()
        {
            // 取りこぼし 20%。
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: new 整合性検査結果(1000, 1000, 1000, 200, 0, 0),
                p_閉鎖検証: Get_裏付けのある閉鎖(),
                p_ポリッシュ: Get_良好な深度(),
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.出力のみ, l_判定.A_品質保証レベル);
            Assert.Contains(未達理由.取りこぼしが多い, l_判定.A_未達理由);
        }

        [Fact]
        public void Get_理由コード_全ての理由に固有のコードが付く()
        {
            var l_コード = Enum.GetValues<未達理由>()
                .Select(CompletenessValidator.Get_理由コード)
                .ToList();

            Assert.Equal(l_コード.Count, l_コード.Distinct().Count());
            Assert.DoesNotContain(l_コード, string.IsNullOrWhiteSpace);
        }
    }
}
