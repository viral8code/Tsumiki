using Tsumiki.Cores.Scaffolding;
using Tsumiki.Core;
using Tsumiki.Models.ContigBuilding;
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

        /// <summary>
        /// contig 1 → 2 → 3 の鎖で contig 2 が多コピーなら、両側の継ぎ目を両向きとも外す
        /// </summary>
        [Fact]
        public void V_外す_多コピーを通り抜ける辺_通り抜ける辺を外す()
        {
            var l_確定辺 = V_構築_鎖の確定辺();

            var l_外した数 = Scaffolder.V_外す_多コピーを通り抜ける辺(l_確定辺, x => x == 2, (_, _) => false);

            Assert.Equal(2, l_外した数);
            Assert.All(l_確定辺, x => Assert.Null(x));
        }

        /// <summary>
        /// 多コピーの contig の片端にだけ確定辺があるなら、残す
        /// </summary>
        [Fact]
        public void V_外す_多コピーを通り抜ける辺_片端だけなら残す()
        {
            var l_確定辺 = new (int A_行き先, int A_ギャップ長)?[8];
            l_確定辺[2] = (4, 0);
            l_確定辺[5] = (3, 0);

            var l_外した数 = Scaffolder.V_外す_多コピーを通り抜ける辺(l_確定辺, x => x == 2, (_, _) => false);

            Assert.Equal(0, l_外した数);
            Assert.Equal((4, 0), l_確定辺[2]);
            Assert.Equal((3, 0), l_確定辺[5]);
        }

        /// <summary>
        /// 通り抜けられる contig が多コピーでなければ、残す
        /// </summary>
        [Fact]
        public void V_外す_多コピーを通り抜ける辺_多コピーでなければ残す()
        {
            var l_確定辺 = V_構築_鎖の確定辺();

            var l_外した数 = Scaffolder.V_外す_多コピーを通り抜ける辺(l_確定辺, _ => false, (_, _) => false);

            Assert.Equal(0, l_外した数);
            Assert.Equal(V_構築_鎖の確定辺(), l_確定辺);
        }

        /// <summary>
        /// 両隣 (contig 1 の出口の頂点 2 から contig 3 の入口の頂点 6) を直接結ぶ支持があれば、多コピーの contig 2 を通り抜ける辺を残す
        /// </summary>
        [Fact]
        public void V_外す_多コピーを通り抜ける辺_両隣を直接結ぶ支持があれば残す()
        {
            var l_確定辺 = V_構築_鎖の確定辺();

            var l_外した数 = Scaffolder.V_外す_多コピーを通り抜ける辺(l_確定辺, x => x == 2, (x, y) => x == 2 && y == 6);

            Assert.Equal(0, l_外した数);
            Assert.Equal(V_構築_鎖の確定辺(), l_確定辺);
        }

        /// <summary>
        /// 多コピーの contig を通らない別の鎖 (contig 4 → 5) は残す
        /// </summary>
        [Fact]
        public void V_外す_多コピーを通り抜ける辺_ほかの鎖は残す()
        {
            var l_確定辺 = new (int A_行き先, int A_ギャップ長)?[12];
            V_構築_鎖の確定辺().CopyTo(l_確定辺, 0);
            l_確定辺[8] = (10, 10);
            l_確定辺[11] = (9, 10);

            var l_外した数 = Scaffolder.V_外す_多コピーを通り抜ける辺(l_確定辺, x => x == 2, (_, _) => false);

            Assert.Equal(2, l_外した数);
            Assert.All(l_確定辺[..8], x => Assert.Null(x));
            Assert.Equal((10, 10), l_確定辺[8]);
            Assert.Equal((9, 10), l_確定辺[11]);
        }

        /// <summary>
        /// 先頭の unitig の逆向きの頂点と、末尾の unitig の walk の向きの頂点がどちらも多コピー未解決の contig だけを返す<br/>
        /// contig 2 (unitig 2 → 3): 頂点 5 (unitig 2 の逆鎖) と 6 (unitig 3 の順鎖) が未解決 → 多コピー。contig 3 (unitig 4): 頂点 8 だけ未解決 → 片端なので多コピーでない
        /// </summary>
        [Fact]
        public void Get_多コピーのcontig_両端が未解決のcontigを返す()
        {
            var l_結果 = Scaffolder.Get_多コピーのcontig(V_構築_配置(), new HashSet<int> { 5, 6, 8 });

            Assert.Equal([2], l_結果);
        }

        /// <summary>
        /// 1 つの unitig でできた contig は、その unitig の順鎖と逆鎖の頂点がどちらも未解決なら多コピー
        /// </summary>
        [Fact]
        public void Get_多コピーのcontig_単独のunitigは両向きが未解決なら返す()
        {
            var l_結果 = Scaffolder.Get_多コピーのcontig(V_構築_配置(), new HashSet<int> { 8, 9, 10 });

            Assert.Equal([3], l_結果);
        }

        /// <summary>
        /// walk で逆鎖として使われた unitig は、向きを入れ替えて見る<br/>
        /// contig 1 = unitig 1 (逆鎖で使用) だけ: 末尾は walk の向きの頂点 3、先頭はその逆の頂点 2
        /// </summary>
        [Fact]
        public void Get_多コピーのcontig_逆鎖で使われたunitigは向きを入れ替えて見る()
        {
            var l_配置 = new Dictionary<int, Unitig配置>
            {
                [1] = new Unitig配置(1, false, 0, 2, true),
                [2] = new Unitig配置(1, false, 1, 2, false)
            };

            Assert.Equal([1], Scaffolder.Get_多コピーのcontig(l_配置, new HashSet<int> { 2, 4 }));
            Assert.Empty(Scaffolder.Get_多コピーのcontig(l_配置, new HashSet<int> { 3, 4 }));
        }

        /// <summary>
        /// 未解決の頂点が無ければ、多コピーの contig は無い
        /// </summary>
        [Fact]
        public void Get_多コピーのcontig_未解決の頂点が無ければ空()
        {
            Assert.Empty(Scaffolder.Get_多コピーのcontig(V_構築_配置(), new HashSet<int>()));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// contig 1 → 2 → 3 の鎖の確定辺 (両向き) を組み立てる
        /// </summary>
        /// <returns>頂点ごとの確定辺</returns>
        private static (int A_行き先, int A_ギャップ長)?[] V_構築_鎖の確定辺()
        {
            var l_確定辺 = new (int A_行き先, int A_ギャップ長)?[8];
            l_確定辺[2] = (4, 0);
            l_確定辺[5] = (3, 0);
            l_確定辺[4] = (6, 0);
            l_確定辺[7] = (5, 0);
            return l_確定辺;
        }

        /// <summary>
        /// unitig 1 (contig 1)、unitig 2 と 3 (contig 2)、unitig 4 (contig 3)、unitig 5 (contig 4) の配置を組み立てる
        /// </summary>
        /// <returns>unitig の配置</returns>
        private static Dictionary<int, Unitig配置> V_構築_配置()
        {
            return new Dictionary<int, Unitig配置>
            {
                [1] = new Unitig配置(1, false, 0, 1, false),
                [2] = new Unitig配置(2, false, 0, 2, false),
                [3] = new Unitig配置(2, false, 1, 2, false),
                [4] = new Unitig配置(3, false, 0, 1, false),
                [5] = new Unitig配置(4, false, 0, 1, false)
            };
        }

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
