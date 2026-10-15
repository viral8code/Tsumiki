using Tsumiki.Cores.Scaffolding;
using Tsumiki.Core;
using Tsumiki.Models.Scaffolding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// scaffold 辺の採用条件
    /// </summary>
    public class ScaffolderEdgeTests
    {
        #region 公開メソッド

        /// <summary>
        /// 明確に優勢な候補を採用することを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_明確に優勢な候補を採用する()
        {
            var l_結果 = Scaffolder.Get_優勢な候補([V_構築_scaffold候補(4, 100UL, 0.9D), V_構築_scaffold候補(6, 12UL, 0.05D)], p_優勢閾値: 0.8M, p_最小証拠数: 10UL);
            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 競合する候補が拮抗しているときは採用しないことを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_競合する候補が拮抗しているときは採用しない()
        {
            var l_候補 = new List<Scaffold候補>
            {
                V_構築_scaffold候補(4, 20UL, 0.5D),
                V_構築_scaffold候補(6, 18UL, 0.45D)
            };
            Assert.Null(Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 10UL));
        }

        /// <summary>
        /// どの候補も最小証拠数に満たないときは採用しないことを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_最小証拠数に満たない候補は採用しない()
        {
            Assert.Null(Scaffolder.Get_優勢な候補([V_構築_scaffold候補(4, 9UL, 0.9D)], p_優勢閾値: 0.8M, p_最小証拠数: 10UL));
        }

        /// <summary>
        /// 優劣は生の本数ではなく期待本数に対する比で決めること
        /// </summary>
        [Fact]
        public void Get_優勢な候補_生の本数ではなく期待本数に対する比で優劣を決める()
        {
            var l_結果 = Scaffolder.Get_優勢な候補([V_構築_scaffold候補(4, 200UL, 0.05D), V_構築_scaffold候補(6, 12UL, 0.95D)], p_優勢閾値: 0.8M, p_最小証拠数: 10UL);
            Assert.NotNull(l_結果);
            Assert.Equal(6, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 本数が少なくても、その場所で期待される本数に見合っていれば採ること
        /// </summary>
        [Fact]
        public void Get_優勢な候補_期待本数に見合っていれば本数が少なくても採用する()
        {
            var l_結果 = Scaffolder.Get_優勢な候補([V_構築_scaffold候補(4, 10UL, 0.95D)], p_優勢閾値: 0.8M, p_最小証拠数: 10UL);
            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 期待本数がほぼ 0 の候補があるときは比ではなく生の支持数で比べ、拮抗していれば採用しないことを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_期待本数が小さい候補があれば生の支持数で拮抗を採用しない()
        {
            var l_候補 = new List<Scaffold候補>
            {
                V_構築_scaffold候補(4, 20UL, 119.4D, 0.2D),
                V_構築_scaffold候補(6, 14UL, 0.0D, 0.0D)
            };
            Assert.Null(Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 3UL));
        }

        /// <summary>
        /// 期待本数が小さい候補があるときは、生の支持数で優勢なら採用することを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_期待本数が小さい候補があれば生の支持数で優勢な候補を採用する()
        {
            var l_候補 = new List<Scaffold候補>
            {
                V_構築_scaffold候補(4, 100UL, 0.9D, 50D),
                V_構築_scaffold候補(6, 12UL, 0.05D, 0.5D)
            };
            var l_結果 = Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 10UL);
            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 生の支持数で判定するとき、最良が支持数の合計に対して優勢でなければ採用しないことを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_生の支持数で3候補が拮抗するときは採用しない()
        {
            var l_候補 = new List<Scaffold候補>
            {
                V_構築_scaffold候補(4, 17UL, 0.0D, 0.0D),
                V_構築_scaffold候補(6, 19UL, 0.0D, 0.0D),
                V_構築_scaffold候補(8, 15UL, 2.57D, 5.8D)
            };
            Assert.Null(Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 3UL));
        }

        /// <summary>
        /// 期待本数がすべて 1 以上なら、生の支持数ではなく従来どおり比で判定することを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_期待本数が下限以上なら比で判定する()
        {
            var l_候補 = new List<Scaffold候補>
            {
                V_構築_scaffold候補(4, 30UL, 6.0D, 5.0D),
                V_構築_scaffold候補(6, 5UL, 0.5D, 10.0D)
            };
            var l_結果 = Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 3UL);
            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 支持数の下限未満の候補は、期待本数が小さくても判定に入れないことを検証する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_支持数の下限未満の候補は判定に入れない()
        {
            var l_候補 = new List<Scaffold候補>
            {
                V_構築_scaffold候補(4, 20UL, 5.0D, 4.0D),
                V_構築_scaffold候補(6, 2UL, 0.0D, 0.0D)
            };
            var l_結果 = Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 3UL);
            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 生の支持数で判定するかどうかの判定を検証する
        /// </summary>
        [Fact]
        public void Is生の支持数で判定_期待本数が下限未満の候補があるときだけ真になる()
        {
            Assert.True(Scaffolder.Is生の支持数で判定([V_構築_scaffold候補(4, 20UL, 119.4D, 0.2D), V_構築_scaffold候補(6, 14UL, 0.0D, 0.0D)], 3UL));
            Assert.False(Scaffolder.Is生の支持数で判定([V_構築_scaffold候補(4, 30UL, 6.0D, 5.0D), V_構築_scaffold候補(6, 5UL, 0.5D, 10.0D)], 3UL));
            Assert.False(Scaffolder.Is生の支持数で判定([V_構築_scaffold候補(4, 20UL, 5.0D, 4.0D), V_構築_scaffold候補(6, 2UL, 0.0D, 0.0D)], 3UL));
        }

        /// <summary>
        /// 比で最良の候補の期待本数が足りていて、期待本数の小さい競合の支持数が最良より少なければ、今までどおり比で判定する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_期待本数の小さい競合が最良より少ない支持なら比で判定する()
        {
            var l_候補 = new List<Scaffold候補>
            {
                V_構築_scaffold候補(4, 21UL, 6.0D, 3.5D),
                V_構築_scaffold候補(6, 7UL, 0.0D, 0.0D)
            };
            Assert.False(Scaffolder.Is生の支持数で判定(l_候補, 3UL));
            var l_結果 = Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 3UL);
            Assert.NotNull(l_結果);
            Assert.Equal(4, l_結果!.Value.A_行き先);
        }

        /// <summary>
        /// 比で最良の候補の期待本数が足りていても、期待本数の小さい競合が最良以上の支持数を持てば、生の支持数で判定する
        /// </summary>
        [Fact]
        public void Get_優勢な候補_期待本数の小さい競合が最良以上の支持なら生の支持数で判定する()
        {
            var l_候補 = new List<Scaffold候補>
            {
                V_構築_scaffold候補(8, 15UL, 2.57D, 5.8D),
                V_構築_scaffold候補(4, 17UL, 0.0D, 0.0D)
            };
            Assert.True(Scaffolder.Is生の支持数で判定(l_候補, 3UL));
            Assert.Null(Scaffolder.Get_優勢な候補(l_候補, p_優勢閾値: 0.8M, p_最小証拠数: 3UL));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// scaffold の候補を組み立てる
        /// </summary>
        /// <param name="p_行き先">繋ぐ先の頂点</param>
        /// <param name="p_支持数">支持したペアの数</param>
        /// <param name="p_期待比">理想本数に対する比</param>
        /// <param name="p_期待本数">理想本数モデルによる期待本数 (既定は NaN、つまり理想本数モデルが無い)</param>
        /// <returns>scaffold の候補</returns>
        private static Scaffold候補 V_構築_scaffold候補(int p_行き先, ulong p_支持数, double p_期待比, double p_期待本数 = double.NaN)
        {
            return new Scaffold候補(p_行き先, p_支持数, 300, p_期待比, p_期待本数);
        }

        #endregion
    }
}
