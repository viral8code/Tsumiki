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
        /// walk 順は unitig の walk 順の位置の順に並ぶことを検証する
        /// </summary>
        [Fact]
        public void Get_contig別のwalk順_walk順の位置の順に並ぶ()
        {
            var l_結果 = Scaffolder.Get_contig別のwalk順(V_構築_配置(false));
            Assert.Equal(new[] { 1, 2, 3 }, l_結果[1]);
            Assert.Equal(new[] { 4 }, l_結果[2]);
        }

        /// <summary>
        /// 出口側の端の unitig は距離 0 で、従来の末端頂点と同じ頂点を返すことを検証する
        /// </summary>
        [Fact]
        public void Is成功_Contig端の近くの頂点_出口側の端のunitigは距離0で従来の頂点を返す()
        {
            var l_結果 = V_Is端の近く(V_構築_配置(false), 3, p_Is出口側: true, p_距離の上限: 0);
            Assert.True(l_結果.A_成功);
            Assert.Equal(2, l_結果.A_頂点);
            Assert.Equal(0, l_結果.A_距離);
        }

        /// <summary>
        /// 内側の unitig は、上限以内なら端までの距離を (unitig 長 − 重なり長) の和で求めることを検証する
        /// </summary>
        [Fact]
        public void Is成功_Contig端の近くの頂点_内側のunitigは上限以内なら距離を和で求める()
        {
            var l_配置 = V_構築_配置(false);
            var l_結果1 = V_Is端の近く(l_配置, 2, p_Is出口側: true, p_距離の上限: 30);
            Assert.True(l_結果1.A_成功);
            Assert.Equal(2, l_結果1.A_頂点);
            Assert.Equal(26, l_結果1.A_距離);
            Assert.False(V_Is端の近く(l_配置, 2, p_Is出口側: true, p_距離の上限: 20).A_成功);

            var l_結果2 = V_Is端の近く(l_配置, 1, p_Is出口側: true, p_距離の上限: 40);
            Assert.True(l_結果2.A_成功);
            Assert.Equal(2, l_結果2.A_頂点);
            Assert.Equal(32, l_結果2.A_距離);
        }

        /// <summary>
        /// 実効の向きが逆鎖のとき、出口側は walk の前側へ進み、頂点は逆鎖の頂点になることを検証する
        /// </summary>
        [Fact]
        public void Is成功_Contig端の近くの頂点_出口側で逆鎖の内側のunitigは前側へ進む()
        {
            var l_結果 = V_Is端の近く(V_構築_配置(false), -2, p_Is出口側: true, p_距離の上限: 30);
            Assert.True(l_結果.A_成功);
            Assert.Equal(3, l_結果.A_頂点);
            Assert.Equal(16, l_結果.A_距離);
        }

        /// <summary>
        /// 入口側の内側の unitig は、実効の向きに応じて前側または後ろ側へ進むことを検証する
        /// </summary>
        [Fact]
        public void Is成功_Contig端の近くの頂点_入口側は実効の向きに応じて進む向きを決める()
        {
            var l_配置 = V_構築_配置(false);
            var l_順鎖 = V_Is端の近く(l_配置, 2, p_Is出口側: false, p_距離の上限: 30);
            Assert.True(l_順鎖.A_成功);
            Assert.Equal(2, l_順鎖.A_頂点);
            Assert.Equal(16, l_順鎖.A_距離);

            var l_逆鎖 = V_Is端の近く(l_配置, -2, p_Is出口側: false, p_距離の上限: 30);
            Assert.True(l_逆鎖.A_成功);
            Assert.Equal(3, l_逆鎖.A_頂点);
            Assert.Equal(26, l_逆鎖.A_距離);
        }

        /// <summary>
        /// contig が逆相補で書き出されたときは、頂点の向きが反転することを検証する
        /// </summary>
        [Fact]
        public void Is成功_Contig端の近くの頂点_逆相補で書き出されたcontigでは頂点の向きが反転する()
        {
            var l_配置 = V_構築_配置(true);
            var l_端 = V_Is端の近く(l_配置, 3, p_Is出口側: true, p_距離の上限: 0);
            Assert.True(l_端.A_成功);
            Assert.Equal(3, l_端.A_頂点);

            var l_内側 = V_Is端の近く(l_配置, 2, p_Is出口側: true, p_距離の上限: 30);
            Assert.True(l_内側.A_成功);
            Assert.Equal(3, l_内側.A_頂点);
            Assert.Equal(26, l_内側.A_距離);
        }

        /// <summary>
        /// 配置が無い unitig、または長さが欠けている unitig は偽になることを検証する
        /// </summary>
        [Fact]
        public void Is成功_Contig端の近くの頂点_配置やunitig長が欠けていれば偽を返す()
        {
            var l_配置 = V_構築_配置(false);
            Assert.False(V_Is端の近く(l_配置, 9, p_Is出口側: true, p_距離の上限: 100).A_成功);

            var l_長さが欠けている = V_構築_unitig長();
            _ = l_長さが欠けている.Remove(3);
            Assert.False(V_Is端の近く(l_配置, 2, p_Is出口側: true, p_距離の上限: 100, p_unitig長: l_長さが欠けている).A_成功);
        }

        /// <summary>
        /// 同じ contig の内側の unitig を指したペアは、端に対応づけられず、両端が端ならそのペアが辺になることを検証する
        /// </summary>
        [Fact]
        public void Get_対称化した辺_内側のunitigを指したペアは既知長に距離を足して辺に数える()
        {
            var l_ペア経路 = new Dictionary<(int, int), List<int>>
            {
                [(2, 4)] = [100, 110]
            };
            var l_内部 = 0;
            var l_未配置 = 0;
            var l_結果 = Scaffolder.Get_対称化した辺(V_構築_配置(false), l_ペア経路, Scaffolder.Get_contig別のwalk順(V_構築_配置(false)), V_構築_unitig長(), 4, 30, new HashSet<int>(), ref l_内部, ref l_未配置);
            Assert.Equal(0, l_内部);
            Assert.Equal(0, l_未配置);
            Assert.Equal(2UL, l_結果[(2, 4)].A_支持数);
            Assert.Equal(new[] { 126, 136 }, l_結果[(2, 4)].A_既知長標本.ToArray());
            Assert.Equal(2UL, l_結果[(5, 3)].A_支持数);
            Assert.Equal(new[] { 126, 136 }, l_結果[(5, 3)].A_既知長標本.ToArray());
        }

        /// <summary>
        /// 内側の unitig を端に対応づけない contig に属する unitig は、距離 0 のときだけ数えることを検証する
        /// </summary>
        [Fact]
        public void Get_対称化した辺_除外したcontigの内側のunitigは数えない()
        {
            var l_内部 = 0;
            var l_未配置 = 0;
            var l_ペア経路 = new Dictionary<(int, int), List<int>>
            {
                [(2, 4)] = [100, 110]
            };
            var l_結果 = Scaffolder.Get_対称化した辺(V_構築_配置(false), l_ペア経路, Scaffolder.Get_contig別のwalk順(V_構築_配置(false)), V_構築_unitig長(), 4, 30, new HashSet<int> { 1 }, ref l_内部, ref l_未配置);
            Assert.False(l_結果.ContainsKey((2, 4)));
            Assert.Equal(1, l_内部);

            var l_端のペア = new Dictionary<(int, int), List<int>>
            {
                [(3, 4)] = [100]
            };
            var l_内部2 = 0;
            var l_未配置2 = 0;
            var l_結果2 = Scaffolder.Get_対称化した辺(V_構築_配置(false), l_端のペア, Scaffolder.Get_contig別のwalk順(V_構築_配置(false)), V_構築_unitig長(), 4, 30, new HashSet<int> { 1 }, ref l_内部2, ref l_未配置2);
            Assert.Equal(1UL, l_結果2[(2, 4)].A_支持数);
            Assert.Equal(new[] { 100 }, l_結果2[(2, 4)].A_既知長標本.ToArray());
            Assert.Equal(0, l_内部2);
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

        /// <summary>
        /// 手で組んだ unitig の配置を作る。contig 1 は unitig 1, 2, 3 (長さ 20, 10, 30) の順鎖、contig 2 は unitig 4 (長さ 40) だけ
        /// </summary>
        /// <param name="p_IsContig逆相補">contig 1 を逆相補で書き出したことにするか</param>
        /// <returns>unitig ID → 配置</returns>
        private static Dictionary<int, Unitig配置> V_構築_配置(bool p_IsContig逆相補)
        {
            return new Dictionary<int, Unitig配置>
            {
                [1] = new Unitig配置(1, p_IsContig逆相補, 0, 3, false),
                [2] = new Unitig配置(1, p_IsContig逆相補, 1, 3, false),
                [3] = new Unitig配置(1, p_IsContig逆相補, 2, 3, false),
                [4] = new Unitig配置(2, false, 0, 1, false)
            };
        }

        /// <summary>
        /// 手で組んだ unitig 長 (unitig 1, 2, 3, 4 の順に 20, 10, 30, 40)
        /// </summary>
        /// <returns>unitig ID → unitig 長</returns>
        private static Dictionary<int, int> V_構築_unitig長()
        {
            return new Dictionary<int, int> { [1] = 20, [2] = 10, [3] = 30, [4] = 40 };
        }

        /// <summary>
        /// 端の近くの判定を、重なり長 4 で呼び出す
        /// </summary>
        /// <param name="p_配置">unitig の配置</param>
        /// <param name="p_符号付きunitigID">符号付き unitig ID</param>
        /// <param name="p_Is出口側">出口側か</param>
        /// <param name="p_距離の上限">距離の上限</param>
        /// <param name="p_unitig長">unitig 長、既定は V_構築_unitig長 の値</param>
        /// <returns>判定の結果</returns>
        private static (bool A_成功, int A_頂点, int A_距離) V_Is端の近く(Dictionary<int, Unitig配置> p_配置, int p_符号付きunitigID, bool p_Is出口側, int p_距離の上限, Dictionary<int, int>? p_unitig長 = null)
        {
            var l_成功 = Scaffolder.Is成功_Contig端の近くの頂点(p_配置, Scaffolder.Get_contig別のwalk順(p_配置), p_unitig長 ?? V_構築_unitig長(), 4, p_符号付きunitigID, p_Is出口側, p_距離の上限, out var l_頂点, out var l_距離);
            return (l_成功, l_頂点, l_距離);
        }

        #endregion
    }
}
