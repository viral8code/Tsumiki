using Tsumiki.Cores.Scaffolding;
using Tsumiki.Core;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// スキャフォールド辺の採用条件
    /// </summary>
    /// <remarks>
    /// 観測本数ではなく期待本数に対する比で測ることを固定する<br/>
    /// 辺が長く距離が近いほど多く観測されるという幾何的な偏りがあるため、固定の下限では期待が数本の場所と数百本の場所を同じ物差しで測ってしまう
    /// </remarks>
    public class ScaffolderEdgeTests
    {
        #region 公開メソッド

        /// <summary>
        /// 明確に優勢な候補を採用することを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_明確に優勢な候補を採用する()
        {
            var l_結果 = Scaffolder.Get_優勢な候補([V_構築_スキャフォールド候補(4, 100UL, 0.9D), V_構築_スキャフォールド候補(6, 12UL, 0.05D)], p_優勢閾値: 0.8M, p_最小証拠数: 10UL);

            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 競合する候補が拮抗しているときは採用しないことを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_競合する候補が拮抗しているときは採用しない()
        {
            var l_候補 = new List<スキャフォールド候補> { V_構築_スキャフォールド候補(4, 20UL, 0.5D), V_構築_スキャフォールド候補(6, 18UL, 0.45D) };

            Assert.Null(Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 10UL));
        }

        /// <summary>
        /// どの候補も最小証拠数に満たないときは採用しないことを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_最小証拠数に満たない候補は採用しない()
        {
            Assert.Null(Scaffolder.Get_優勢な候補([V_構築_スキャフォールド候補(4, 9UL, 0.9D)], p_優勢閾値: 0.8M, p_最小証拠数: 10UL));
        }

        /// <summary>
        /// 優劣は生の本数ではなく期待本数に対する比で決めること
        /// </summary>
        /// <remarks>
        /// 本数は辺が長く距離が近いほど多くなるので、そのまま比べると幾何的に有利なだけの辺が勝ってしまう
        /// </remarks>
        [Fact]
        public void Get_優勢な候補_生の本数ではなく期待本数に対する比で優劣を決める()
        {
            var l_結果 = Scaffolder.Get_優勢な候補([V_構築_スキャフォールド候補(4, 200UL, 0.05D), V_構築_スキャフォールド候補(6, 12UL, 0.95D)], p_優勢閾値: 0.8M, p_最小証拠数: 10UL);

            Assert.NotNull(l_結果);
            Assert.Equal(6, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 本数が少なくても、その場所で期待される本数に見合っていれば採ること
        /// </summary>
        /// <remarks>
        /// 短い contig 同士や広いギャップでは、正しい隣接でも本数は少なくなる
        /// </remarks>
        [Fact]
        public void Get_優勢な候補_期待本数に見合っていれば本数が少なくても採用する()
        {
            var l_結果 = Scaffolder.Get_優勢な候補([V_構築_スキャフォールド候補(4, 10UL, 0.95D)], p_優勢閾値: 0.8M, p_最小証拠数: 10UL);

            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// スキャフォールドの候補を組み立てる
        /// </summary>
        /// <param name="p_行き先">繋ぐ先の頂点</param>
        /// <param name="p_支持数">支持したペアの数</param>
        /// <param name="p_期待比">理想本数に対する比</param>
        /// <returns>スキャフォールドの候補</returns>
        private static スキャフォールド候補 V_構築_スキャフォールド候補(int p_行き先, ulong p_支持数, double p_期待比)
        {
            return new スキャフォールド候補(p_行き先, p_支持数, 300, p_期待比);
        }

        #endregion

    }
}
