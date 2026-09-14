using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// カバレッジが右へ長く裾を引くスペクトルでの解析の検証
    /// </summary>
    public class KmerHistogramSkewTests
    {
        #region 公開メソッド

        /// <summary>
        /// 山 (最頻値) が平均より大きく下に来ても、ゲノムサイズを水増ししない
        /// </summary>
        /// <remarks>
        /// GC の偏りが強いライブラリでは単一コピーのカバレッジが対数正規分布のように広がる
        /// </remarks>
        [Fact]
        public void V_裾の長いスペクトルでもゲノムサイズを水増ししない()
        {
            const long l_真のゲノムサイズ = 5_000_000L;
            var l_ヒストグラム = V_構築_対数正規スペクトル(p_中央値: 120D, p_σ: 0.45D, p_ゲノムサイズ: l_真のゲノムサイズ, p_エラー係数: 10_000_000L);

            var l_解析 = KmerHistogram.Get_解析結果(l_ヒストグラム);

            Assert.NotNull(l_解析);

            // 延べ数を山の位置で割る推定がこの形では大きく外れることを前提として確かめておく
            Assert.True(l_解析.A_ゲノム由来の延べ数 / (double)l_解析.A_ピーク出現回数 > 1.3D * l_真のゲノムサイズ);
            Assert.InRange(l_解析.A_推定ゲノムサイズ, (long)(0.85D * l_真のゲノムサイズ), (long)(1.25D * l_真のゲノムサイズ));
            Assert.InRange(l_解析.A_単一コピー基準値, 100D, 130D);
            Assert.True(l_解析.A_単一コピー上限 > 1.5D * l_解析.A_単一コピー基準値);
        }

        /// <summary>
        /// 単一コピー上限未満は 1、以上は基準値との比を丸めた 2 以上のコピー数になる
        /// </summary>
        [Fact]
        public void V_期待コピー数は上限未満で1になる()
        {
            Assert.Equal(1, KmerHistogram.Get_期待コピー数(199D, 100D, 200D));
            Assert.Equal(2, KmerHistogram.Get_期待コピー数(200D, 100D, 200D));
            Assert.Equal(4, KmerHistogram.Get_期待コピー数(410D, 100D, 200D));
            Assert.Equal(1, KmerHistogram.Get_期待コピー数(500D, 0D, 0D));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// エラー由来の裾と、対数正規分布状に広がった単一コピーの山を重ねたヒストグラムを作る
        /// </summary>
        /// <param name="p_中央値"></param>
        /// <param name="p_σ">対数の標準偏差</param>
        /// <param name="p_ゲノムサイズ"></param>
        /// <param name="p_エラー係数"></param>
        /// <returns></returns>
        private static Dictionary<ulong, long> V_構築_対数正規スペクトル(double p_中央値, double p_σ, long p_ゲノムサイズ, long p_エラー係数)
        {
            var l_ヒストグラム = new Dictionary<ulong, long>();
            var l_μ = Math.Log(p_中央値);
            for (var l_出現回数 = 1UL; l_出現回数 <= 2_000UL; l_出現回数++)
            {
                var c = (double)l_出現回数;
                var l_エラー = (long)(p_エラー係数 / (c * c));
                var l_ゲノム = (long)(p_ゲノムサイズ * Math.Exp(-Math.Pow(Math.Log(c) - l_μ, 2D) / (2D * p_σ * p_σ)) / (c * p_σ * Math.Sqrt(2D * Math.PI)));
                l_ヒストグラム[l_出現回数] = l_エラー + l_ゲノム;
            }
            return l_ヒストグラム;
        }

        #endregion
    }
}
