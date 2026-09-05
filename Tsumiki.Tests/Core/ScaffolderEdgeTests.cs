using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// スキャフォールド辺の採用条件。反復が作る偽の隣接を落とすための
    /// 2つの判定(優勢比の分母・支持ペアの距離の一貫性)を固定する。
    /// </summary>
    public class ScaffolderEdgeTests
    {
        private static (int, ulong, List<int>) 候補(int p_行き先, ulong p_支持数, params int[] p_標本)
        {
            return (p_行き先, p_支持数, [.. p_標本]);
        }

        [Fact]
        public void Get_優勢な候補_TakesTheDominantEdge()
        {
            var l_結果 = Scaffolder.Get_優勢な候補(
                [候補(4, 100, 300), 候補(6, 3, 300)], p_優勢閾値: 0.8m, p_最小証拠数: 10);

            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 支持数を満たす対抗辺が拮抗しているときは、どちらが正しいかの
        /// 根拠が無いので採用しない。
        /// </summary>
        [Fact]
        public void Get_優勢な候補_RejectsWhenRivalsAreComparable()
        {
            var l_候補 = new List<(int, ulong, List<int>)>
            {
                候補(4, 20, 300), 候補(6, 18, 300),
            };

            Assert.Null(Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8m, p_最小証拠数: 10));
        }

        /// <summary>
        /// 支持数に満たない対抗辺は分母から外す。ペアの支持は反復が捏造するため、
        /// 弱い対抗辺を曖昧さの証拠として扱うと連結が過度に失われる。
        /// </summary>
        [Fact]
        public void Get_優勢な候補_IgnoresRivalsBelowTheSupportFloor()
        {
            var l_候補 = new List<(int, ulong, List<int>)>
            {
                候補(4, 20, 300), 候補(6, 9, 300), 候補(8, 9, 300), 候補(10, 9, 300),
            };

            var l_結果 = Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8m, p_最小証拠数: 10);
            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        [Fact]
        public void Get_優勢な候補_RejectsWhenNoCandidateMeetsTheSupportFloor()
        {
            Assert.Null(Scaffolder.Get_優勢な候補(
                [候補(4, 9, 300)], p_優勢閾値: 0.8m, p_最小証拠数: 10));
        }

    }
}
