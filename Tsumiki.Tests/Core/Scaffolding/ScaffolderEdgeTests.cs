using Tsumiki.Core.Scaffolding;
using Tsumiki.Core;
using Tsumiki.Model.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// スキャフォールド辺の採用条件<br/>
    /// 観測本数ではなく期待本数に対する比で
    /// 測ることを固定する<br/>
    /// 辺が長く距離が近いほど多く観測されるという幾何的な
    /// 偏りがあるため、固定の下限では期待が数本の場所と数百本の場所を
    /// 同じ物差しで測ってしまう
    /// </summary>
    public class ScaffolderEdgeTests
    {
        /// <summary>
        /// スキャフォールドの候補を組み立てる
        /// </summary>
        /// <param name="p_行き先">繋ぐ先の頂点</param>
        /// <param name="p_支持数">支持したペアの数</param>
        /// <param name="p_期待比">理想本数に対する比</param>
        /// <returns>スキャフォールドの候補</returns>
        private static スキャフォールド候補 候補(int p_行き先, ulong p_支持数, double p_期待比)
        {
            return new スキャフォールド候補(p_行き先, p_支持数, 300, p_期待比);
        }

        [Fact]
        public void Get_優勢な候補_TakesTheDominantEdge()
        {
            var l_結果 = Scaffolder.Get_優勢な候補(
                [候補(4, 100, 0.9), 候補(6, 12, 0.05)], p_優勢閾値: 0.8M, p_最小証拠数: 10);

            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        [Fact]
        public void Get_優勢な候補_RejectsWhenRivalsAreComparable()
        {
            var l_候補 = new List<スキャフォールド候補> { 候補(4, 20, 0.5), 候補(6, 18, 0.45) };

            Assert.Null(Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 10));
        }

        [Fact]
        public void Get_優勢な候補_RejectsWhenNoCandidateMeetsTheSupportFloor()
        {
            Assert.Null(Scaffolder.Get_優勢な候補(
                [候補(4, 9, 0.9)], p_優勢閾値: 0.8M, p_最小証拠数: 10));
        }

        /// <summary>
        /// 優劣は生の本数ではなく期待本数に対する比で決めること<br/>
        /// 本数は辺が長く距離が近いほど多くなるので、そのまま比べると
        /// 幾何的に有利なだけの辺が勝ってしまう
        /// </summary>
        [Fact]
        public void Get_優勢な候補_RanksByTheExpectedCountRatioNotTheRawCount()
        {
            var l_結果 = Scaffolder.Get_優勢な候補(
                [候補(4, 200, 0.05), 候補(6, 12, 0.95)], p_優勢閾値: 0.8M, p_最小証拠数: 10);

            Assert.NotNull(l_結果);
            Assert.Equal(6, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 本数が少なくても、その場所で期待される本数に見合っていれば採ること<br/>
        /// 短い contig 同士や広いギャップでは、正しい隣接でも本数は少なくなる
        /// </summary>
        [Fact]
        public void Get_優勢な候補_AcceptsAFewPairsWhenThatIsAllThatIsExpected()
        {
            var l_結果 = Scaffolder.Get_優勢な候補(
                [候補(4, 10, 0.95)], p_優勢閾値: 0.8M, p_最小証拠数: 10);

            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }
    }
}
