using Tsumiki.Cores.Evaluation;
using Tsumiki.Models.Reporting;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 決めきれなかった箇所の記録を固定する
    /// </summary>
    public class AmbiguityRecorderTests
    {
        #region 定数

        /// <summary>
        /// 塩基配列 ACGTACGT
        /// </summary>
        private const string C_塩基配列_ACGTACGT = "ACGTACGT";

        /// <summary>
        /// 塩基配列 TTTTAAAA
        /// </summary>
        private const string C_塩基配列_TTTTAAAA = "TTTTAAAA";

        /// <summary>
        /// 塩基配列 ACGTACGA
        /// </summary>
        private const string C_塩基配列_ACGTACGA = "ACGTACGA";

        /// <summary>
        /// 項目 再実行前アンカー左
        /// </summary>
        private const string C_項目_再実行前アンカー左 = "再実行前アンカー左";

        /// <summary>
        /// 項目 再実行前アンカー右
        /// </summary>
        private const string C_項目_再実行前アンカー右 = "再実行前アンカー右";

        /// <summary>
        /// 項目 scaffold0 0 10
        /// </summary>
        private const string C_項目_scaffold0_0_10 = "scaffold0:0-10";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 同じアンカーの組からは常に同じ安定 ID が得られることを確かめる
        /// </summary>
        [Fact]
        public void Get_安定ID_同じアンカーなら同じIDになる()
        {
            var l_ID1 = AmbiguityRecorder.Get_安定ID(C_塩基配列_ACGTACGT, C_塩基配列_TTTTAAAA);
            var l_ID2 = AmbiguityRecorder.Get_安定ID(C_塩基配列_ACGTACGT, C_塩基配列_TTTTAAAA);
            Assert.Equal(l_ID1, l_ID2);
        }

        /// <summary>
        /// アンカーが違えば別の ID になることを確かめる
        /// </summary>
        [Fact]
        public void Get_安定ID_アンカーが違えば別のIDになる()
        {
            var l_ID1 = AmbiguityRecorder.Get_安定ID(C_塩基配列_ACGTACGT, C_塩基配列_TTTTAAAA);
            var l_ID2 = AmbiguityRecorder.Get_安定ID(C_塩基配列_ACGTACGA, C_塩基配列_TTTTAAAA);
            Assert.NotEqual(l_ID1, l_ID2);
        }

        /// <summary>
        /// 同じ k を再実行して現在の記録が上書きされても、履歴には再実行前の判定が残ることを確かめる
        /// </summary>
        [Fact]
        public void V_開始_同じkの再実行でも履歴は上書きされない()
        {
            const int l_k長 = 90211;
            var l_安定ID = AmbiguityRecorder.Get_安定ID(C_項目_再実行前アンカー左, C_項目_再実行前アンカー右);
            AmbiguityRecorder.V_開始(l_k長);
            AmbiguityRecorder.V_記録(曖昧箇所の種別.到達不能, C_項目_scaffold0_0_10, p_安定ID: l_安定ID);
            Assert.Contains(AmbiguityRecorder.Get_記録(l_k長), x => x.A_安定ID == l_安定ID);
            AmbiguityRecorder.V_開始(l_k長);
            Assert.DoesNotContain(AmbiguityRecorder.Get_記録(l_k長), x => x.A_安定ID == l_安定ID);
            Assert.Contains(AmbiguityRecorder.Get_履歴(), x => x.A_安定ID == l_安定ID);
        }

        #endregion
    }
}
