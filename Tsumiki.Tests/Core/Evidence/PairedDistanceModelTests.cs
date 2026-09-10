using Tsumiki.Cores.Evidence;
using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアの隣接証拠を期待本数との比で測るためのモデル
    /// </summary>
    /// <remarks>
    /// 観測本数をそのまま固定の下限と比べると幾何的な偏りを拾うため、
    /// 期待位置数の算出と、裾に強いギャップ長推定を固定する
    /// </remarks>
    public class PairedDistanceModelTests
    {
        /// <summary>
        /// 中央 400、おおよそ 350-450 に広がるフラグメント長分布
        /// </summary>
        private static List<int> Get_分布(int p_件数 = 2000)
        {
            var l_乱数 = new Random(4649);
            return [.. Enumerable.Range(0, p_件数).Select(_ => 400 + l_乱数.Next(-50, 51))];
        }

        /// <summary>
        /// 検証に使う断片長のモデル
        /// </summary>
        /// <returns>断片長のモデル</returns>
        private static PairedDistanceModel Get_モデル() => new(Get_分布(), p_リード長: 100);

        /// <summary>
        /// きれいなクラスタからギャップ長を見つけられることを確かめる
        /// </summary>
        [Fact]
        public void Get_一貫した支持_きれいなクラスタからギャップ長を見つける()
        {
            var l_モデル = Get_モデル();
            // ギャップ 100 なら既知長は 400-100 = 300 付近に集まる
            List<int> l_標本 = [295, 300, 302, 298, 305, 300];

            var (l_本数, l_ギャップ) = l_モデル.Get_一貫した支持(l_標本);

            Assert.Equal(6, l_本数);
            Assert.InRange(l_ギャップ, 90, 110);
        }

        /// <summary>
        /// 誤マップ由来の裾が過半を占めても、峰から出るギャップ長と本数が
        /// 変わらないこと
        /// </summary>
        /// <remarks>
        /// 中央値で測るとここが壊れる
        /// </remarks>
        [Fact]
        public void Get_一貫した支持_誤マップの重い裾があってもギャップ長と本数は変わらない()
        {
            var l_モデル = Get_モデル();
            List<int> l_峰 = [295, 300, 302, 298, 305, 300];
            List<int> l_裾 = [2300, 4100, 6907, 8000, 11269, 15000, 3050, 5200];

            var (l_本数, l_ギャップ) = l_モデル.Get_一貫した支持([.. l_峰, .. l_裾]);

            Assert.Equal(6, l_本数);
            Assert.InRange(l_ギャップ, 90, 110);
        }

        /// <summary>
        /// 標本が空なら支持本数は 0 になることを確かめる
        /// </summary>
        [Fact]
        public void Get_一貫した支持_標本が空なら支持本数は0になる()
        {
            Assert.Equal(0, Get_モデル().Get_一貫した支持([]).A_本数);
        }

        /// <summary>
        /// 期待位置数は、接合点の両側が短いほど、ギャップが広いほど小さくなること
        /// </summary>
        /// <remarks>
        /// 観測本数を固定の下限と比べてはいけない理由そのもの
        /// </remarks>
        [Fact]
        public void Get_期待位置数_両側が短いほどギャップが広いほど小さくなる()
        {
            var l_モデル = Get_モデル();

            var l_短い辺 = l_モデル.Get_期待位置数(200, 200, 0);
            var l_長い辺 = l_モデル.Get_期待位置数(50_000, 50_000, 0);
            Assert.True(l_長い辺 > l_短い辺);

            var l_広いギャップ = l_モデル.Get_期待位置数(50_000, 50_000, 200);
            Assert.True(l_長い辺 > l_広いギャップ);
        }

        /// <summary>
        /// 接合点から 1 フラグメント長ぶんの窓しか寄与しないので、それより長い
        /// 配列では期待位置数は増えない
        /// </summary>
        [Fact]
        public void Get_期待位置数_フラグメント長を超えると増えなくなる()
        {
            var l_モデル = Get_モデル();

            Assert.Equal(
                l_モデル.Get_期待位置数(50_000, 50_000, 0),
                l_モデル.Get_期待位置数(5_000, 5_000, 0),
                6);
        }

        /// <summary>
        /// フラグメント長を超えるギャップは跨げないので期待は 0 になること
        /// </summary>
        [Fact]
        public void Get_期待位置数_フラグメント長を超えるギャップは0になる()
        {
            Assert.Equal(0, Get_モデル().Get_期待位置数(50_000, 50_000, 1_000), 6);
        }

        /// <summary>
        /// フラグメントより短い配列では期待位置数が 0 になることを確かめる
        /// </summary>
        [Fact]
        public void Get_期待位置数_単一_フラグメントより短い配列では0になる()
        {
            var l_モデル = Get_モデル();

            Assert.Equal(0, l_モデル.Get_期待位置数_単一(100), 6);
            Assert.True(l_モデル.Get_期待位置数_単一(10_000) > 9_000);
        }

        /// <summary>
        /// モデルは標本が無ければ使えないことを確かめる
        /// </summary>
        [Fact]
        public void モデルは標本が無ければ使えない()
        {
            Assert.False(new PairedDistanceModel([], p_リード長: 100).A_使えるか);
            Assert.True(Get_モデル().A_使えるか);
        }
    }
}
