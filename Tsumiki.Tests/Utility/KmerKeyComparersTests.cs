using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer の値と配置候補の比較器が、既定の比較と同じ同値性と列挙順を保つことの検証
    /// </summary>
    public class KmerKeyComparersTests
    {
        #region 公開メソッド

        /// <summary>
        /// UInt128組比較器_等しい値は等しくハッシュも一致する
        /// </summary>
        [Fact]
        public void V_UInt128組比較器_等しい値は等しくハッシュも一致する()
        {
            var l_x = Get_組(new UInt128(1UL, 2UL), new UInt128(3UL, 4UL));
            var l_y = Get_組(new UInt128(1UL, 2UL), new UInt128(3UL, 4UL));
            Assert.True(UInt128組比較器.C_既定.Equals(l_x, l_y));
            Assert.Equal(UInt128組比較器.C_既定.GetHashCode(l_x), UInt128組比較器.C_既定.GetHashCode(l_y));
        }

        /// <summary>
        /// UInt128組比較器_成分が1つだけ違えば等しくない
        /// </summary>
        /// <param name="p_上位上"></param>
        /// <param name="p_上位下"></param>
        /// <param name="p_下位上"></param>
        /// <param name="p_下位下"></param>
        [Theory]
        [InlineData(0UL, 2UL, 3UL, 4UL)]
        [InlineData(1UL, 0UL, 3UL, 4UL)]
        [InlineData(1UL, 2UL, 0UL, 4UL)]
        [InlineData(1UL, 2UL, 3UL, 0UL)]
        public void V_UInt128組比較器_成分が1つだけ違えば等しくない(ulong p_上位上, ulong p_上位下, ulong p_下位上, ulong p_下位下)
        {
            var l_x = Get_組(new UInt128(1UL, 2UL), new UInt128(3UL, 4UL));
            var l_y = Get_組(new UInt128(p_上位上, p_上位下), new UInt128(p_下位上, p_下位下));
            Assert.False(UInt128組比較器.C_既定.Equals(l_x, l_y));
        }

        /// <summary>
        /// UInt128組比較器_上位だけ異なる値は等しくない
        /// </summary>
        [Fact]
        public void V_UInt128組比較器_上位だけ異なる値は等しくない()
        {
            var l_x = Get_組(new UInt128(1UL << 63, 0UL), new UInt128(0UL, 1UL));
            var l_y = Get_組(new UInt128(1UL << 62, 0UL), new UInt128(0UL, 1UL));
            Assert.False(UInt128組比較器.C_既定.Equals(l_x, l_y));
        }

        /// <summary>
        /// UInt128組比較器_乱数のキーで既定と同じ列挙順になる
        /// </summary>
        [Fact]
        public void V_UInt128組比較器_乱数のキーで既定と同じ列挙順になる()
        {
            var l_乱数 = new Random(1);
            var l_キー群 = Enumerable.Range(0, 100_000).Select(_ => Get_組(Get_乱数128(l_乱数), Get_乱数128(l_乱数))).ToList();
            var l_既定 = new Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>();
            var l_新 = new Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>(UInt128組比較器.C_既定);
            for (var i = 0; i < 60_000; i++)
            {
                l_既定[l_キー群[i]] = (ulong)i;
                l_新[l_キー群[i]] = (ulong)i;
            }

            for (var i = 0; i < 60_000; i += 3)
            {
                _ = l_既定.Remove(l_キー群[i]);
                _ = l_新.Remove(l_キー群[i]);
            }

            for (var i = 60_000; i < l_キー群.Count; i++)
            {
                l_既定[l_キー群[i]] = (ulong)i;
                l_新[l_キー群[i]] = (ulong)i;
            }

            Assert.Equal(l_既定.Keys.ToList(), l_新.Keys.ToList());
            Assert.Equal(l_既定.Values.ToList(), l_新.Values.ToList());
        }

        /// <summary>
        /// 配置候補比較器_等しい候補は等しくハッシュも一致する
        /// </summary>
        [Fact]
        public void V_配置候補比較器_等しい候補は等しくハッシュも一致する()
        {
            Assert.True(配置候補比較器.C_既定.Equals((3, true, -7), (3, true, -7)));
            Assert.Equal(配置候補比較器.C_既定.GetHashCode((3, true, -7)), 配置候補比較器.C_既定.GetHashCode((3, true, -7)));
        }

        /// <summary>
        /// 配置候補比較器_成分が1つだけ違えば等しくない
        /// </summary>
        /// <param name="p_配列番号"></param>
        /// <param name="p_Is逆鎖"></param>
        /// <param name="p_対角線"></param>
        [Theory]
        [InlineData(4, true, -7)]
        [InlineData(3, false, -7)]
        [InlineData(3, true, -6)]
        public void V_配置候補比較器_成分が1つだけ違えば等しくない(int p_配列番号, bool p_Is逆鎖, int p_対角線)
        {
            Assert.False(配置候補比較器.C_既定.Equals((3, true, -7), (p_配列番号, p_Is逆鎖, p_対角線)));
        }

        /// <summary>
        /// 配置候補比較器_乱数の候補で既定と同じ列挙順になる
        /// </summary>
        [Fact]
        public void V_配置候補比較器_乱数の候補で既定と同じ列挙順になる()
        {
            var l_乱数 = new Random(1);
            var l_キー群 = Enumerable.Range(0, 100_000).Select(_ => (l_乱数.Next(0, 200), l_乱数.Next(2) == 0, l_乱数.Next(-500, 500))).ToList();
            var l_既定 = new Dictionary<(int A_配列番号, bool A_Is逆鎖, int A_対角線), int>();
            var l_新 = new Dictionary<(int A_配列番号, bool A_Is逆鎖, int A_対角線), int>(配置候補比較器.C_既定);
            for (var i = 0; i < 60_000; i++)
            {
                l_既定[l_キー群[i]] = i;
                l_新[l_キー群[i]] = i;
            }

            for (var i = 0; i < 60_000; i += 3)
            {
                _ = l_既定.Remove(l_キー群[i]);
                _ = l_新.Remove(l_キー群[i]);
            }

            for (var i = 60_000; i < l_キー群.Count; i++)
            {
                l_既定[l_キー群[i]] = i;
                l_新[l_キー群[i]] = i;
            }

            Assert.Equal(l_既定.Keys.ToList(), l_新.Keys.ToList());
            Assert.Equal(l_既定.Values.ToList(), l_新.Values.ToList());
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 同値性の検証に使う、上位と下位の組を作る
        /// </summary>
        /// <param name="p_上位">上位</param>
        /// <param name="p_下位">下位</param>
        /// <returns>値の組</returns>
        private static (UInt128 A_上位, UInt128 A_下位) Get_組(UInt128 p_上位, UInt128 p_下位)
        {
            return (p_上位, p_下位);
        }

        /// <summary>
        /// 乱数から 128 ビットの値を作る
        /// </summary>
        /// <param name="p_乱数">乱数</param>
        /// <returns>128 ビットの値</returns>
        private static UInt128 Get_乱数128(Random p_乱数)
        {
            return new UInt128((ulong)p_乱数.NextInt64(), (ulong)p_乱数.NextInt64());
        }

        #endregion
    }
}
