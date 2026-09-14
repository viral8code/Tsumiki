using Tsumiki.Cores.Evaluation;
using Tsumiki.Models.Reporting;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 決めきれなかった箇所の記録を固定する
    /// </summary>
    /// <remarks>
    /// 要点は「同じ k を再実行しても、再実行前の判定が履歴から失われない」ことと
    /// 「座標が変わっても同じ箇所を安定IDで追跡できる」こと
    /// </remarks>
    public class AmbiguityRecorderTests
    {
        #region 公開メソッド

        /// <summary>
        /// 同じアンカーの組からは常に同じ安定IDが得られることを確かめる
        /// </summary>
        [Fact]
        public void Get_安定ID_同じアンカーなら同じIDになる()
        {
            var l_ID1 = AmbiguityRecorder.Get_安定ID("ACGTACGT", "TTTTAAAA");
            var l_ID2 = AmbiguityRecorder.Get_安定ID("ACGTACGT", "TTTTAAAA");

            Assert.Equal(l_ID1, l_ID2);
        }

        /// <summary>
        /// アンカーが違えば別のIDになることを確かめる
        /// </summary>
        [Fact]
        public void Get_安定ID_アンカーが違えば別のIDになる()
        {
            var l_ID1 = AmbiguityRecorder.Get_安定ID("ACGTACGT", "TTTTAAAA");
            var l_ID2 = AmbiguityRecorder.Get_安定ID("ACGTACGA", "TTTTAAAA");

            Assert.NotEqual(l_ID1, l_ID2);
        }

        /// <summary>
        /// 同じ k を再実行して現在の記録が上書きされても、履歴には再実行前の判定が残ることを確かめる
        /// </summary>
        [Fact]
        public void V_開始_同じkの再実行でも履歴は上書きされない()
        {
            const int l_k長 = 90211;
            var l_安定ID = AmbiguityRecorder.Get_安定ID("再実行前アンカー左", "再実行前アンカー右");

            AmbiguityRecorder.V_開始(l_k長);
            AmbiguityRecorder.V_記録(曖昧箇所の種別.到達不能, "scaffold0:0-10", p_安定ID: l_安定ID);
            Assert.Contains(AmbiguityRecorder.Get_記録(l_k長), x => x.A_安定ID == l_安定ID);

            // 同じ k をもう一度実行すると、現在の記録ビューは上書きされる
            AmbiguityRecorder.V_開始(l_k長);
            Assert.DoesNotContain(AmbiguityRecorder.Get_記録(l_k長), x => x.A_安定ID == l_安定ID);

            // それでも累積履歴には再実行前の判定が残っている
            Assert.Contains(AmbiguityRecorder.Get_履歴(), x => x.A_安定ID == l_安定ID);
        }

        #endregion
    }
}
