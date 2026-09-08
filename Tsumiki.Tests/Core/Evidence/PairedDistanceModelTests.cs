using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアの隣接証拠を期待本数との比で測るためのモデル。
    /// 観測本数をそのまま固定の下限と比べると幾何的な偏りを拾うため、
    /// 期待位置数の算出と、裾に強いギャップ長推定を固定する。
    /// </summary>
    public class PairedDistanceModelTests
    {
        /// <summary>中央 400、おおよそ 350-450 に広がるフラグメント長分布。</summary>
        private static List<int> Get_分布(int p_件数 = 2000)
        {
            var l_乱数 = new Random(4649);
            return [.. Enumerable.Range(0, p_件数).Select(_ => 400 + l_乱数.Next(-50, 51))];
        }

        private static PairedDistanceModel Get_モデル() => new(Get_分布(), p_リード長: 100);

        [Fact]
        public void Get_一貫した支持_FindsTheGapFromACleanCluster()
        {
            var l_モデル = Get_モデル();
            // ギャップ 100 なら既知長は 400-100 = 300 付近に集まる。
            List<int> l_標本 = [295, 300, 302, 298, 305, 300];

            var (l_本数, l_ギャップ) = l_モデル.Get_一貫した支持(l_標本);

            Assert.Equal(6, l_本数);
            Assert.InRange(l_ギャップ, 90, 110);
        }

        /// <summary>
        /// 誤マップ由来の裾が過半を占めても、峰から出るギャップ長と本数が
        /// 変わらないこと。中央値で測るとここが壊れる。
        /// </summary>
        [Fact]
        public void Get_一貫した支持_IsNotMovedByAHeavyTailOfMismappedPairs()
        {
            var l_モデル = Get_モデル();
            List<int> l_峰 = [295, 300, 302, 298, 305, 300];
            List<int> l_裾 = [2300, 4100, 6907, 8000, 11269, 15000, 3050, 5200];

            var (l_本数, l_ギャップ) = l_モデル.Get_一貫した支持([.. l_峰, .. l_裾]);

            Assert.Equal(6, l_本数);
            Assert.InRange(l_ギャップ, 90, 110);
        }

        [Fact]
        public void Get_一貫した支持_ReportsNoSupportForAnEmptySample()
        {
            Assert.Equal(0, Get_モデル().Get_一貫した支持([]).A_本数);
        }

        /// <summary>
        /// 期待位置数は、接合点の両側が短いほど、ギャップが広いほど小さくなること。
        /// 観測本数を固定の下限と比べてはいけない理由そのもの。
        /// </summary>
        [Fact]
        public void Get_期待位置数_ShrinksWithShortFlanksAndWideGaps()
        {
            var l_モデル = Get_モデル();

            var l_短い辺 = l_モデル.Get_期待位置数(200, 200, 0);
            var l_長い辺 = l_モデル.Get_期待位置数(50_000, 50_000, 0);
            Assert.True(l_長い辺 > l_短い辺);

            var l_広いギャップ = l_モデル.Get_期待位置数(50_000, 50_000, 200);
            Assert.True(l_長い辺 > l_広いギャップ);
        }

        /// <summary>
        /// 接合点から1フラグメント長ぶんの窓しか寄与しないので、それより長い
        /// 配列では期待位置数は増えない。
        /// </summary>
        [Fact]
        public void Get_期待位置数_SaturatesBeyondTheFragmentLength()
        {
            var l_モデル = Get_モデル();

            Assert.Equal(
                l_モデル.Get_期待位置数(50_000, 50_000, 0),
                l_モデル.Get_期待位置数(5_000, 5_000, 0),
                6);
        }

        /// <summary>
        /// フラグメント長を超えるギャップは跨げないので期待は 0 になること。
        /// </summary>
        [Fact]
        public void Get_期待位置数_IsZeroBeyondTheFragmentLength()
        {
            Assert.Equal(0, Get_モデル().Get_期待位置数(50_000, 50_000, 1_000), 6);
        }

        [Fact]
        public void Get_期待位置数_単一_IsZeroForASequenceShorterThanTheFragment()
        {
            var l_モデル = Get_モデル();

            Assert.Equal(0, l_モデル.Get_期待位置数_単一(100), 6);
            Assert.True(l_モデル.Get_期待位置数_単一(10_000) > 9_000);
        }

        [Fact]
        public void Model_IsUnusableWithoutSamples()
        {
            Assert.False(new PairedDistanceModel([], p_リード長: 100).A_使えるか);
            Assert.True(Get_モデル().A_使えるか);
        }
    }
}
