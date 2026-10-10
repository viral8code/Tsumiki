using Tsumiki.Cores.Evidence;
using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペア証拠を生カウントではなく期待本数との比で測るための較正器の検証
    /// </summary>
    public class EvidenceCalibratorTests
    {
        #region 公開メソッド

        /// <summary>
        /// リード長が不明な較正器は使えないこと
        /// </summary>
        [Fact]
        public void V_較正器_リード長が不明なら使えない()
        {
            var l_較正器 = 証拠較正器.Get_較正器(Get_同一unitig標本(), p_リード長: null, [1_000L, 2_000L]);
            Assert.False(l_較正器.A_Is使用可能);
        }

        /// <summary>
        /// 標本が空な較正器は使えないこと
        /// </summary>
        [Fact]
        public void V_較正器_標本が空なら使えない()
        {
            var l_較正器 = 証拠較正器.Get_較正器([], p_リード長: 100, [1_000L, 2_000L]);
            Assert.False(l_較正器.A_Is使用可能);
        }

        /// <summary>
        /// 全 unitig がフラグメントより短い較正器は使えないこと
        /// </summary>
        [Fact]
        public void V_較正器_全unitigがフラグメントより短いと使えない()
        {
            var l_較正器 = 証拠較正器.Get_較正器(Get_同一unitig標本(), p_リード長: 100, [50L, 80L]);
            Assert.False(l_較正器.A_Is使用可能);
        }

        /// <summary>
        /// 較正器が使えない場合の正規化済み支持はゼロになること
        /// </summary>
        [Fact]
        public void V_正規化済み支持_較正器が使えないとゼロ()
        {
            var l_較正器 = 証拠較正器.Get_較正器([], p_リード長: 100, [1_000L, 2_000L]);
            Assert.Equal(0D, l_較正器.Get_正規化済み支持(p_観測本数: 10UL, p_基準長: 1_000L, p_比較長: 1_000L, p_ギャップ長: 0));
        }

        /// <summary>
        /// 較正器の核心: 同じ観測本数でも、短い unitig への辺は正規化後の値が大きくなる (=期待が小さいので少ない観測でも強い支持とみなされる)
        /// </summary>
        [Fact]
        public void V_正規化済み支持_同じ観測本数でも隣接unitigが短いほど大きい()
        {
            var l_標本 = Get_同一unitig標本();
            var l_較正器 = 証拠較正器.Get_較正器(l_標本, p_リード長: 100, [100_000L, 150L, 50_000L]);
            var l_短い辺への支持 = l_較正器.Get_正規化済み支持(p_観測本数: 8UL, p_基準長: 100_000L, p_比較長: 150L, p_ギャップ長: 0);
            var l_長い辺への支持 = l_較正器.Get_正規化済み支持(p_観測本数: 8UL, p_基準長: 100_000L, p_比較長: 50_000L, p_ギャップ長: 0);
            Assert.True(l_較正器.A_Is使用可能);
            Assert.True(l_短い辺への支持 > l_長い辺への支持, $"expected short-flank ratio ({l_短い辺への支持}) to exceed long-flank ratio ({l_長い辺への支持}) for the same raw count");
        }

        /// <summary>
        /// 観測本数がちょうど期待本数どおりなら比はおよそ 1.0 になること
        /// </summary>
        [Fact]
        public void V_正規化済み支持_観測本数が理想どおりならおよそ1()
        {
            var l_標本 = Get_同一unitig標本(20_000);
            IReadOnlyList<long> l_unitig長一覧 = [50_000L, 50_000L, 50_000L];
            var l_較正器 = 証拠較正器.Get_較正器(l_標本, p_リード長: 100, l_unitig長一覧);
            var l_モデル = new PairedDistanceModel(l_標本, p_リード長: 100);
            var l_期待位置数 = l_モデル.Get_期待位置数(50_000L, 50_000L, 0);
            var l_密度相当の観測本数 = (ulong)Math.Round(l_標本.Count / (3D * l_モデル.Get_期待位置数_単一(50_000L)) * l_期待位置数);
            var l_比 = l_較正器.Get_正規化済み支持(l_密度相当の観測本数, 50_000L, 50_000L, 0);
            Assert.InRange(l_比, 0.9D, 1.1D);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 中央 400 付近に集まるフラグメント長の標本
        /// </summary>
        /// <param name="p_件数"></param>
        /// <returns></returns>
        private static List<int> Get_同一unitig標本(int p_件数 = 2_000)
        {
            var l_乱数 = new Random(20_260_908);
            return[..Enumerable.Range(0, p_件数).Select(_ => 400 + l_乱数.Next(-50, 51))];
        }

        #endregion
    }
}
