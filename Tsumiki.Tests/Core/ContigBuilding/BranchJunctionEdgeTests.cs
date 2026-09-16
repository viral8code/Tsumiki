using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// contig が分岐のある継ぎ目で通った辺を取り出す処理の検証
    /// </summary>
    public class BranchJunctionEdgeTests
    {
        #region 定数

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 8;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 入次数 2 の頂点へ入る継ぎ目だけを (k+1)-mer で返し、一本道の継ぎ目は返さない
        /// </summary>
        [Fact]
        public void V_入次数2の頂点へ入る継ぎ目だけを返す()
        {
            var l_反復 = RepeatGraphFixture.Get_乱数配列(16, 5101);
            var l_入口1 = RepeatGraphFixture.Get_乱数配列(11, 5102) + "A" + l_反復[..(k長 - 1)];
            var l_入口2 = RepeatGraphFixture.Get_乱数配列(11, 5103) + "C" + l_反復[..(k長 - 1)];
            var l_出口 = l_反復[^(k長 - 1)..] + "G" + RepeatGraphFixture.Get_乱数配列(11, 5104);
            var (l_グラフ, l_unitig配列) = RepeatGraphFixture.Get_グラフ(k長, l_反復, l_入口1, l_入口2, l_出口);
            Assert.Equal(2, l_グラフ.Get_入次数(Get_頂点(1)));
            Assert.Equal([Get_頂点(4)], l_グラフ.A_出辺[Get_頂点(1)]);

            var l_辺 = ContigMaker.Get_分岐の継ぎ目(l_グラフ, l_unitig配列, [Get_頂点(2), Get_頂点(1), Get_頂点(4)], k長).ToList();

            Assert.Equal([l_入口1[^k長..] + l_反復[k長 - 1]], l_辺);
        }

        /// <summary>
        /// 出次数 2 の頂点から出る継ぎ目を返す
        /// </summary>
        [Fact]
        public void V_出次数2の頂点から出る継ぎ目を返す()
        {
            var l_反復 = RepeatGraphFixture.Get_乱数配列(16, 5201);
            var l_入口 = RepeatGraphFixture.Get_乱数配列(11, 5202) + "A" + l_反復[..(k長 - 1)];
            var l_出口1 = l_反復[^(k長 - 1)..] + "G" + RepeatGraphFixture.Get_乱数配列(11, 5203);
            var l_出口2 = l_反復[^(k長 - 1)..] + "T" + RepeatGraphFixture.Get_乱数配列(11, 5204);
            var (l_グラフ, l_unitig配列) = RepeatGraphFixture.Get_グラフ(k長, l_反復, l_入口, l_出口1, l_出口2);
            Assert.Equal(2, l_グラフ.A_出辺[Get_頂点(1)].Count);

            var l_辺 = ContigMaker.Get_分岐の継ぎ目(l_グラフ, l_unitig配列, [Get_頂点(2), Get_頂点(1), Get_頂点(4)], k長).ToList();

            Assert.Equal([l_反復[^k長..] + l_出口2[k長 - 1]], l_辺);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// unitig ID の順鎖の頂点番号
        /// </summary>
        /// <param name="p_ID"></param>
        /// <returns></returns>
        private static int Get_頂点(int p_ID)
        {
            return ContigMaker.Get_頂点番号(p_ID);
        }

        #endregion
    }
}
