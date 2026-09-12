using Tsumiki.Cores.Evaluation;
using Tsumiki.Core;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Polishing;
using Tsumiki.Models.Reporting;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 完全長の判定を固定する
    /// </summary>
    /// <remarks>
    /// 要点は「材料が無いことを合格にしない」ことと「不合格と判定不能を混同しない」こと<br/>
    /// どちらを崩しても、根拠の無い完全長が通ってしまう
    /// </remarks>
    public class CompletenessValidatorTests
    {
        #region 公開メソッド

        /// <summary>
        /// 全ての検査を通れば完全長になることを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_全ての検査を通れば完全長になる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: Get_良好な自己検査(), p_閉鎖検証: Get_裏付けのある閉鎖(), p_ポリッシュ: Get_良好な深度(), p_曖昧箇所: [], p_支持検査: Get_良好な支持());

            Assert.True(l_判定.A_完全長か);
            Assert.Equal(品質保証レベル.完全長, l_判定.A_品質保証レベル);
            Assert.Empty(l_判定.A_未達理由);
        }

        /// <summary>
        /// リード支持が不合格または未検査なら完全長と判定しないことを確かめる
        /// </summary>
        /// <param name="p_検査済みか">支持のない位置を検出した検査結果を渡すか</param>
        [Theory]
        [InlineData(true)]
        [InlineData(false)]
        public void Get_判定結果_リード支持が合格しなければグラフ整合で止まる(bool p_検査済みか)
        {
            支持検査結果? l_支持 = p_検査済みか ? new 支持検査結果(31, 100000L, 1L, [new 支持のない区間("contig", 50, 50)]) : null;
            var l_判定 = CompletenessValidator.Get_判定結果(0, Get_良好な自己検査(), Get_裏付けのある閉鎖(), Get_良好な深度(), [], l_支持);

            Assert.False(l_判定.A_完全長か);
            Assert.Equal(品質保証レベル.グラフ整合, l_判定.A_品質保証レベル);
            Assert.Equal(p_検査済みか ? 検査判定.不合格 : 検査判定.判定不能, Get_判定(l_判定, "read_support"));
            Assert.Contains(p_検査済みか ? 未達理由.リードに裏付けの無い箇所がある : 未達理由.リードの支持を調べていない, l_判定.A_未達理由);
        }

        /// <summary>
        /// 閉じ目を検証していなければ完全長にはしないことを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_閉じ目を検証していなければ完全長にはしない()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: Get_良好な自己検査(), p_閉鎖検証: null, p_ポリッシュ: Get_良好な深度(), p_曖昧箇所: [], p_支持検査: Get_良好な支持());

            Assert.False(l_判定.A_完全長か);
            Assert.Equal(品質保証レベル.接合点が支持済み, l_判定.A_品質保証レベル);
            Assert.Equal(検査判定.判定不能, Get_判定(l_判定, "circular_closure"));
            Assert.Contains(未達理由.環状閉鎖を検証していない, l_判定.A_未達理由);
        }

        /// <summary>
        /// 閉じ目に裏付けが無い場合は判定不能ではなく不合格として区別することを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_閉じ目に裏付けが無い場合は不合格として区別する()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: Get_良好な自己検査(), p_閉鎖検証: [new 環状閉鎖検証結果("scaffold1_circular", 1000, 1, 5)], p_ポリッシュ: Get_良好な深度(), p_曖昧箇所: [], p_支持検査: Get_良好な支持());

            Assert.Equal(検査判定.不合格, Get_判定(l_判定, "circular_closure"));
            Assert.Contains(未達理由.閉じ目がリードで裏付けられない, l_判定.A_未達理由);
        }

        /// <summary>
        /// 環状の配列が 1 本も無ければ不合格になることを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_環状の配列が1本も無ければ不合格になる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: Get_良好な自己検査(), p_閉鎖検証: [], p_ポリッシュ: Get_良好な深度(), p_曖昧箇所: [], p_支持検査: Get_良好な支持());

            Assert.Equal(検査判定.不合格, Get_判定(l_判定, "circular_closure"));
            Assert.Contains(未達理由.環状に閉じていない, l_判定.A_未達理由);
        }

        /// <summary>
        /// 未解決のギャップが残っていればペア整合の段階で止まることを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_未解決のギャップが残っていればペア整合で止まる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 3, p_整合性: Get_良好な自己検査(), p_閉鎖検証: Get_裏付けのある閉鎖(), p_ポリッシュ: Get_良好な深度(), p_曖昧箇所: [], p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.マッピング整合, l_判定.A_品質保証レベル);
            Assert.Contains(未達理由.未解決のギャップが残る, l_判定.A_未達理由);
        }

        /// <summary>
        /// 決めきれない分岐が残っていれば接合点の支持の段階で止まることを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_決めきれない分岐が残っていれば接合点の支持で止まる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: Get_良好な自己検査(), p_閉鎖検証: Get_裏付けのある閉鎖(), p_ポリッシュ: Get_良好な深度(), p_曖昧箇所: [new 曖昧箇所(63, 曖昧箇所の種別.僅差, "unitig7+", 1.2D, 1.1D, 9L, 0.95D)], p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.ペア整合, l_判定.A_品質保証レベル);
            Assert.Equal(検査判定.不合格, Get_判定(l_判定, "junction_support"));
            Assert.Equal(検査判定.不合格, Get_判定(l_判定, "no_alternative_path"));
        }

        /// <summary>
        /// 深度を測っていなければグラフ整合の段階で止まることを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_深度を測っていなければグラフ整合で止まる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: Get_良好な自己検査(), p_閉鎖検証: Get_裏付けのある閉鎖(), p_ポリッシュ: null, p_曖昧箇所: [], p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.グラフ整合, l_判定.A_品質保証レベル);
            Assert.Equal(検査判定.判定不能, Get_判定(l_判定, "coverage_continuity"));
            Assert.Contains(未達理由.深度を測っていない, l_判定.A_未達理由);
        }

        /// <summary>
        /// 自己検査ができていなければ出力のみに留まることを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_自己検査ができていなければ出力のみに留まる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: null, p_閉鎖検証: Get_裏付けのある閉鎖(), p_ポリッシュ: Get_良好な深度(), p_曖昧箇所: [], p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.出力のみ, l_判定.A_品質保証レベル);
            Assert.Equal(検査判定.判定不能, Get_判定(l_判定, "graph_coverage"));
            Assert.Contains(未達理由.自己検査を行えなかった, l_判定.A_未達理由);
        }

        /// <summary>
        /// 取りこぼしが多ければグラフ被覆の検査で落ちることを確かめる
        /// </summary>
        [Fact]
        public void Get_判定結果_取りこぼしが多ければグラフ被覆で落ちる()
        {
            // 取りこぼし 20%
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: new 整合性検査結果(1000L, 1000L, 1000L, 200L, 0L, 0L), p_閉鎖検証: Get_裏付けのある閉鎖(), p_ポリッシュ: Get_良好な深度(), p_曖昧箇所: [], p_支持検査: Get_良好な支持());

            Assert.Equal(品質保証レベル.出力のみ, l_判定.A_品質保証レベル);
            Assert.Contains(未達理由.取りこぼしが多い, l_判定.A_未達理由);
        }

        /// <summary>
        /// 全ての未達理由に固有のコードが付くことを確かめる
        /// </summary>
        [Fact]
        public void Get_理由コード_全ての理由に固有のコードが付く()
        {
            var l_コード = Enum.GetValues<未達理由>()
                .Select(CompletenessValidator.Get_理由コード)
                .ToList();

            Assert.Equal(l_コード.Count, l_コード.Distinct().Count());
            Assert.DoesNotContain(l_コード, string.IsNullOrWhiteSpace);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// リードに裏付けの無い位置が一つも無い検査結果
        /// </summary>
        /// <returns>支持検査の結果</returns>
        private static 支持検査結果 Get_良好な支持()
        {
            return new 支持検査結果(A_r長: 31, A_調べた位置数: 100000L, A_支持のない位置数: 0L, A_区間: []);
        }

        /// <summary>
        /// 取りこぼしも出しすぎも許容内に収まる自己検査の結果
        /// </summary>
        /// <returns>自己検査の結果</returns>
        private static 整合性検査結果 Get_良好な自己検査()
        {
            // 取りこぼし 1%、出しすぎ 0%
            return new 整合性検査結果(1000L, 1000L, 1000L, 10L, 0L, 0L);
        }

        /// <summary>
        /// 深度の落ち込みが許容内に収まるポリッシュの結果
        /// </summary>
        /// <returns>ポリッシュの結果</returns>
        private static ポリッシュ統計 Get_良好な深度()
        {
            return new ポリッシュ統計(1, 1000L, 1000L, 0L, 0L, 0L, 1000L, 80D);
        }

        /// <summary>
        /// 閉じ目をリードが跨いでいる環状閉鎖の検証結果
        /// </summary>
        /// <returns>環状閉鎖の検証結果</returns>
        private static IReadOnlyList<環状閉鎖検証結果> Get_裏付けのある閉鎖()
        {
            return [new 環状閉鎖検証結果("scaffold1_circular", 1000, 12, 5)];
        }

        /// <summary>
        /// 完全長の判定から、指定した検査項目の判定を取り出す
        /// </summary>
        /// <param name="p_判定">完全長の判定結果</param>
        /// <param name="p_キー">取り出す検査項目のキー</param>
        /// <returns>その項目の判定</returns>
        private static 検査判定 Get_判定(完全性判定結果 p_判定, string p_キー)
        {
            return p_判定.A_検査項目.Single(x => x.A_キー == p_キー).A_判定;
        }

        #endregion

    }
}
