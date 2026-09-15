using Tsumiki.Core;
using Tsumiki.Cores.UnitigBuilding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 入口と出口が 2 本に限らない反復を、通り抜けたリードの並びで解く検証
    /// </summary>
    public class MultiWayRepeatResolutionTests
    {
        #region 定数

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 8;

        /// <summary>
        /// 優勢閾値
        /// </summary>
        private const decimal 優勢閾値 = 0.8M;

        /// <summary>
        /// 最小証拠数
        /// </summary>
        private const ulong 最小証拠数 = 10UL;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 入口 3・出口 3 の反復を、並びが示す 3 組に解く
        /// </summary>
        [Fact]
        public void V_入口3出口3の反復を並びどおりに解く()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(3, 3, 3001);
            var l_反復 = 頂点(1);
            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([入口(0), l_反復, 出口(3, 1)], 20UL);
            l_索引.V_追加([入口(1), l_反復, 出口(3, 2)], 20UL);
            l_索引.V_追加([入口(2), l_反復, 出口(3, 0)], 20UL);

            var l_解決数 = V_解決(l_グラフ, l_unitig配列, l_索引);

            Assert.Equal(1, l_解決数);
            V_検証_一本道(l_グラフ, 入口(0), 出口(3, 1));
            V_検証_一本道(l_グラフ, 入口(1), 出口(3, 2));
            V_検証_一本道(l_グラフ, 入口(2), 出口(3, 0));
        }

        /// <summary>
        /// 入口 2・出口 3 で、支持の無い出口 (行き止まり) は元の反復に残す
        /// </summary>
        [Fact]
        public void V_入口2出口3で支持の無い出口は元の反復に残す()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(2, 3, 3101);
            var l_反復 = 頂点(1);
            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([入口(0), l_反復, 出口(2, 0)], 20UL);
            l_索引.V_追加([入口(1), l_反復, 出口(2, 1)], 20UL);

            var l_解決数 = V_解決(l_グラフ, l_unitig配列, l_索引);

            Assert.Equal(1, l_解決数);
            var l_複製 = RepeatGraphFixture.Get_唯一の次(l_グラフ, 入口(0));
            var l_元 = RepeatGraphFixture.Get_唯一の次(l_グラフ, 入口(1));
            Assert.NotEqual(l_複製, l_元);

            // 入口が 1 本も残らないので、決着した組の 1 つは支持の無い出口と一緒に元の反復へ残る
            var l_残った側 = l_グラフ.A_出辺[l_複製].Count == 1 ? l_元 : l_複製;
            var l_移った側 = l_残った側 == l_元 ? l_複製 : l_元;
            Assert.Single(l_グラフ.A_出辺[l_移った側]);
            Assert.Equal(2, l_グラフ.A_出辺[l_残った側].Count);
            Assert.Contains(出口(2, 2), l_グラフ.A_出辺[l_残った側]);
        }

        /// <summary>
        /// 1 組だけ決着したら、その組だけを複製へ移し、残りは元の反復に残す
        /// </summary>
        [Fact]
        public void V_1組だけ決着したらその組だけを移す()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(3, 3, 3201);
            var l_反復 = 頂点(1);
            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([入口(0), l_反復, 出口(3, 0)], 20UL);
            l_索引.V_追加([入口(1), l_反復, 出口(3, 1)], 10UL);
            l_索引.V_追加([入口(1), l_反復, 出口(3, 2)], 10UL);
            l_索引.V_追加([入口(2), l_反復, 出口(3, 1)], 10UL);
            l_索引.V_追加([入口(2), l_反復, 出口(3, 2)], 10UL);

            var l_解決数 = V_解決(l_グラフ, l_unitig配列, l_索引);

            Assert.Equal(1, l_解決数);
            V_検証_一本道(l_グラフ, 入口(0), 出口(3, 0));
            Assert.Equal(2, l_グラフ.Get_入次数(l_反復));
            Assert.Equal(2, l_グラフ.A_出辺[l_反復].Count);
        }

        /// <summary>
        /// 支持が五分五分なら触らない
        /// </summary>
        [Fact]
        public void V_支持が五分五分なら触らない()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(2, 2, 3301);
            var l_反復 = 頂点(1);
            var l_索引 = new ReadPathIndex();
            foreach (var l_入口 in new[] { 入口(0), 入口(1) })
            {
                foreach (var l_出口 in new[] { 出口(2, 0), 出口(2, 1) })
                {
                    l_索引.V_追加([l_入口, l_反復, l_出口], 10UL);
                }
            }
            var l_頂点数 = l_グラフ.A_出辺.Count;

            Assert.Equal(0, V_解決(l_グラフ, l_unitig配列, l_索引));
            Assert.Equal(l_頂点数, l_グラフ.A_出辺.Count);
        }

        /// <summary>
        /// 並びとペアが別の対応付けを示すなら触らない
        /// </summary>
        [Fact]
        public void V_並びとペアが食い違えば触らない()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(2, 2, 3401);
            var l_反復 = 頂点(1);
            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([入口(0), l_反復, 出口(2, 0)], 20UL);
            l_索引.V_追加([入口(1), l_反復, 出口(2, 1)], 20UL);
            Dictionary<(int, int), ulong> l_ペア連結 = new() { [(入口(0), 出口(2, 1))] = 30UL, [(入口(1), 出口(2, 0))] = 30UL };
            var l_頂点数 = l_グラフ.A_出辺.Count;

            Assert.Equal(0, l_グラフ.V_解決_短い反復(l_unitig配列, [], l_ペア連結, 500, 優勢閾値, 最小証拠数, p_経路索引: l_索引));
            Assert.Equal(l_頂点数, l_グラフ.A_出辺.Count);
        }

        /// <summary>
        /// 1 組だけ決着して入口と出口が 1 本ずつ残るとき、その残る組が他の入口と出口を争っているなら反復を解かない
        /// </summary>
        /// <remarks>
        /// 残る出口へ、複製へ移す入口からも抜けたリードがある (出口が 2 コピーある) 形<br/>
        /// 解くと、決着していない残る組を walk が一本道として通り抜ける
        /// </remarks>
        [Fact]
        public void V_残る1組が他の入口と出口を争っていれば解かない()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(2, 2, 4101);
            var l_反復 = 頂点(1);
            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([入口(0), l_反復, 出口(2, 1)], 11UL);
            l_索引.V_追加([入口(1), l_反復, 出口(2, 0)], 9UL);
            l_索引.V_追加([入口(1), l_反復, 出口(2, 1)], 6UL);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(入口(0), 出口(2, 1))] = 18UL,
                [(入口(1), 出口(2, 0))] = 27UL,
                [(入口(1), 出口(2, 1))] = 5UL,
            };
            var l_頂点数 = l_グラフ.A_出辺.Count;

            Assert.Equal(0, l_グラフ.V_解決_短い反復(l_unitig配列, [], l_ペア連結, 500, 優勢閾値, 最小証拠数, p_経路索引: l_索引));
            Assert.Equal(l_頂点数, l_グラフ.A_出辺.Count);
        }

        /// <summary>
        /// ペアだけで全組が決着しても、通り抜けたリードの並びがその組と争っているなら反復を解かない
        /// </summary>
        [Fact]
        public void V_ペアで決着した組を並びが争っていれば解かない()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(2, 2, 4201);
            var l_反復 = 頂点(1);
            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([入口(0), l_反復, 出口(2, 1)], 11UL);
            l_索引.V_追加([入口(1), l_反復, 出口(2, 0)], 9UL);
            l_索引.V_追加([入口(1), l_反復, 出口(2, 1)], 6UL);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(入口(0), 出口(2, 1))] = 24UL,
                [(入口(1), 出口(2, 0))] = 27UL,
            };
            var l_頂点数 = l_グラフ.A_出辺.Count;

            Assert.Equal(0, l_グラフ.V_解決_短い反復(l_unitig配列, [], l_ペア連結, 500, 優勢閾値, 最小証拠数, p_経路索引: l_索引));
            Assert.Equal(l_頂点数, l_グラフ.A_出辺.Count);
        }

        /// <summary>
        /// 入口側の分岐と出口側の分岐が一本道で結ばれた反復の鎖を、鎖ごと複製して解く
        /// </summary>
        [Fact]
        public void V_複数unitigにまたがる反復を鎖ごと解く()
        {
            var l_反復1 = RepeatGraphFixture.Get_乱数配列(16, 3501);
            var l_中間 = l_反復1[^(k長 - 1)..] + RepeatGraphFixture.Get_乱数配列(10, 3502);
            var l_反復2 = l_中間[^(k長 - 1)..] + RepeatGraphFixture.Get_乱数配列(10, 3503);
            var l_入1 = RepeatGraphFixture.Get_乱数配列(11, 3504) + "A" + l_反復1[..(k長 - 1)];
            var l_入2 = RepeatGraphFixture.Get_乱数配列(11, 3505) + "C" + l_反復1[..(k長 - 1)];
            var l_出1 = l_反復2[^(k長 - 1)..] + "G" + RepeatGraphFixture.Get_乱数配列(11, 3506);
            var l_出2 = l_反復2[^(k長 - 1)..] + "T" + RepeatGraphFixture.Get_乱数配列(11, 3507);
            var (l_グラフ, l_unitig配列) = RepeatGraphFixture.Get_グラフ(k長, l_反復1, l_中間, l_反復2, l_入1, l_入2, l_出1, l_出2);
            Assert.Equal(2, l_グラフ.Get_入次数(頂点(1)));
            Assert.Equal([頂点(2)], l_グラフ.A_出辺[頂点(1)]);
            Assert.Equal([頂点(3)], l_グラフ.A_出辺[頂点(2)]);
            Assert.Equal(2, l_グラフ.A_出辺[頂点(3)].Count);

            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([頂点(4), 頂点(1), 頂点(2), 頂点(3), 頂点(7)], 20UL);
            l_索引.V_追加([頂点(5), 頂点(1), 頂点(2), 頂点(3), 頂点(6)], 20UL);

            Assert.Equal(1, V_解決(l_グラフ, l_unitig配列, l_索引));

            foreach (var (l_入口, l_出口) in new[] { (頂点(4), 頂点(7)), (頂点(5), 頂点(6)) })
            {
                var l_1 = RepeatGraphFixture.Get_唯一の次(l_グラフ, l_入口);
                var l_2 = RepeatGraphFixture.Get_唯一の次(l_グラフ, l_1);
                var l_3 = RepeatGraphFixture.Get_唯一の次(l_グラフ, l_2);
                Assert.Equal([l_出口], l_グラフ.A_出辺[l_3]);
                Assert.Equal(l_反復2, l_unitig配列[l_3]);
                Assert.Equal(1, l_グラフ.Get_入次数(l_1));
            }
        }

        /// <summary>
        /// A-R-B-R-C のように同じ反復が 2 回現れる形を、A-R-B と B-R-C に解く
        /// </summary>
        [Fact]
        public void V_同じ反復が2回現れる並びを解く()
        {
            var l_反復 = RepeatGraphFixture.Get_乱数配列(16, 3601);
            var l_A = RepeatGraphFixture.Get_乱数配列(11, 3602) + "A" + l_反復[..(k長 - 1)];
            var l_B = l_反復[^(k長 - 1)..] + "G" + RepeatGraphFixture.Get_乱数配列(11, 3603) + "C" + l_反復[..(k長 - 1)];
            var l_C = l_反復[^(k長 - 1)..] + "T" + RepeatGraphFixture.Get_乱数配列(11, 3604);
            var (l_グラフ, l_unitig配列) = RepeatGraphFixture.Get_グラフ(k長, l_反復, l_A, l_B, l_C);
            Assert.Equal(2, l_グラフ.Get_入次数(頂点(1)));
            Assert.Equal(2, l_グラフ.A_出辺[頂点(1)].Count);

            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([頂点(2), 頂点(1), 頂点(3)], 20UL);
            l_索引.V_追加([頂点(3), 頂点(1), 頂点(4)], 20UL);

            Assert.Equal(1, V_解決(l_グラフ, l_unitig配列, l_索引));
            V_検証_一本道(l_グラフ, 頂点(2), 頂点(3));
            V_検証_一本道(l_グラフ, 頂点(3), 頂点(4));
        }

        /// <summary>
        /// この k の証拠が 1 件も無い反復は、前段 k の確定経路の並びに従って解く
        /// </summary>
        [Fact]
        public void V_このkの証拠が無い反復は引き継ぎ経路で解く()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(2, 2, 3701);
            var l_反復 = 頂点(1);
            var l_引き継ぎ = new ReadPathIndex();
            l_引き継ぎ.V_追加([入口(0), l_反復, 出口(2, 1)], 1UL);
            l_引き継ぎ.V_追加([入口(1), l_反復, 出口(2, 0)], 1UL);

            Assert.Equal(1, l_グラフ.V_解決_短い反復(l_unitig配列, [], new Dictionary<(int, int), ulong>(), 500, 優勢閾値, 最小証拠数, p_経路索引: new ReadPathIndex(), p_引き継ぎ経路索引: l_引き継ぎ));
            V_検証_一本道(l_グラフ, 入口(0), 出口(2, 1));
            V_検証_一本道(l_グラフ, 入口(1), 出口(2, 0));
        }

        /// <summary>
        /// この k のリードが前段 k の経路と別の対応付けを示すなら、この k のリードに従う
        /// </summary>
        [Fact]
        public void V_このkのリードは引き継ぎ経路より優先する()
        {
            var (l_グラフ, l_unitig配列) = Get_放射状の反復(2, 2, 3901);
            var l_反復 = 頂点(1);
            var l_索引 = new ReadPathIndex();
            l_索引.V_追加([入口(0), l_反復, 出口(2, 0)], 20UL);
            l_索引.V_追加([入口(1), l_反復, 出口(2, 1)], 20UL);
            var l_引き継ぎ = new ReadPathIndex();
            l_引き継ぎ.V_追加([入口(0), l_反復, 出口(2, 1)], 1UL);
            l_引き継ぎ.V_追加([入口(1), l_反復, 出口(2, 0)], 1UL);

            Assert.Equal(1, l_グラフ.V_解決_短い反復(l_unitig配列, [], new Dictionary<(int, int), ulong>(), 500, 優勢閾値, 最小証拠数, p_経路索引: l_索引, p_引き継ぎ経路索引: l_引き継ぎ));
            V_検証_一本道(l_グラフ, 入口(0), 出口(2, 0));
            V_検証_一本道(l_グラフ, 入口(1), 出口(2, 1));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// unitig ID の順鎖の頂点番号
        /// </summary>
        /// <param name="p_ID"></param>
        /// <returns></returns>
        private static int 頂点(int p_ID)
        {
            return ContigMaker.Get_頂点番号(p_ID);
        }

        /// <summary>
        /// 放射状の反復の i 番目の入口 (ID 2 から)
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        private static int 入口(int i)
        {
            return 頂点(2 + i);
        }

        /// <summary>
        /// 放射状の反復の j 番目の出口 (入口の後ろの ID から)
        /// </summary>
        /// <param name="p_入口数"></param>
        /// <param name="j"></param>
        /// <returns></returns>
        private static int 出口(int p_入口数, int j)
        {
            return 頂点(2 + p_入口数 + j);
        }

        /// <summary>
        /// 反復 (ID 1) に入口と出口を指定数だけ付けたグラフ
        /// </summary>
        /// <param name="p_入口数"></param>
        /// <param name="p_出口数"></param>
        /// <param name="p_乱数種"></param>
        /// <returns></returns>
        private static (UnitigGraph A_グラフ, List<string> A_unitig配列) Get_放射状の反復(int p_入口数, int p_出口数, int p_乱数種)
        {
            // 入口どうし・出口どうしは反復に接する 1 塩基を変えて、末端の k-mer を重ねない
            var l_反復 = RepeatGraphFixture.Get_乱数配列(16, p_乱数種);
            List<string> l_unitig群 = [l_反復];
            for (var i = 0; i < p_入口数; i++)
            {
                l_unitig群.Add(RepeatGraphFixture.Get_乱数配列(11, p_乱数種 + 10 + i) + "ACG"[i] + l_反復[..(k長 - 1)]);
            }
            for (var j = 0; j < p_出口数; j++)
            {
                l_unitig群.Add(l_反復[^(k長 - 1)..] + "GTA"[j] + RepeatGraphFixture.Get_乱数配列(11, p_乱数種 + 20 + j));
            }
            var (l_グラフ, l_unitig配列) = RepeatGraphFixture.Get_グラフ(k長, [.. l_unitig群]);
            Assert.Equal(p_入口数, l_グラフ.Get_入次数(頂点(1)));
            Assert.Equal(p_出口数, l_グラフ.A_出辺[頂点(1)].Count);
            return (l_グラフ, l_unitig配列);
        }

        /// <summary>
        /// 並びの索引だけを証拠にして反復を解く
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_索引"></param>
        /// <returns></returns>
        private static int V_解決(UnitigGraph p_グラフ, List<string> p_unitig配列, ReadPathIndex p_索引)
        {
            return p_グラフ.V_解決_短い反復(p_unitig配列, [], new Dictionary<(int, int), ulong>(), 500, 優勢閾値, 最小証拠数, p_経路索引: p_索引);
        }

        /// <summary>
        /// 入口から反復のコピー 1 つだけを通って出口へ行く一本道になっていることを確かめる
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_入口"></param>
        /// <param name="p_出口"></param>
        private static void V_検証_一本道(UnitigGraph p_グラフ, int p_入口, int p_出口)
        {
            var l_コピー = p_グラフ.A_出辺[p_入口].Single(x => p_グラフ.A_出辺[x].Contains(p_出口));
            Assert.Equal(1, p_グラフ.Get_入次数(l_コピー));
            Assert.Equal([p_出口], p_グラフ.A_出辺[l_コピー]);
        }

        #endregion
    }
}
