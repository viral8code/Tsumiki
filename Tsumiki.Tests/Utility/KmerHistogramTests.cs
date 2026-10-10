using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer スペクトルから谷と山を求める処理の検証
    /// </summary>
    public class KmerHistogramTests
    {
        #region 定数

        /// <summary>
        /// 項目 1 10  2 0  3 5
        /// </summary>
        private const string C_項目_1_10__2_0__3_5 = "1:10, 2:0, 3:5";

        /// <summary>
        /// 項目 1 10  2 5
        /// </summary>
        private const string C_項目_1_10__2_5 = "1:10, 2:5";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 空のヒストグラムでは null を返す
        /// </summary>
        [Fact]
        public void V_空のヒストグラムではnullを返す()
        {
            Assert.Null(KmerHistogram.Get_推奨カットオフ(new Dictionary<ulong, long>()));
        }

        /// <summary>
        /// 典型的な二峰性スペクトルで谷を見つける
        /// </summary>
        [Fact]
        public void V_典型的な二峰性スペクトルで谷を見つける()
        {
            Dictionary<ulong, long> l_histogram = new()
            {
                [1UL] = 10_000L,
                [2UL] = 3_000L,
                [3UL] = 500L,
                [4UL] = 800L,
                [5UL] = 5_000L,
                [30UL] = 20_000L,
                [31UL] = 19_000L,
            };
            var l_analysis = KmerHistogram.Get_解析結果(l_histogram);
            Assert.NotNull(l_analysis);
            Assert.Equal(3UL, l_analysis.A_谷);
            Assert.Equal(30UL, l_analysis.A_ピーク出現回数);
        }

        /// <summary>
        /// 推奨カットオフは谷そのものではない
        /// </summary>
        [Fact]
        public void V_エラーが既に支配していなければ谷より手前で止める()
        {
            Dictionary<ulong, long> l_histogram = new()
            {
                [1UL] = 10_000L,
                [2UL] = 3_000L,
                [3UL] = 500L,
                [4UL] = 800L,
                [5UL] = 5_000L,
                [30UL] = 20_000L,
                [31UL] = 19_000L,
            };
            Assert.Equal(2UL, KmerHistogram.Get_推奨カットオフ(l_histogram));
        }

        /// <summary>
        /// エラー由来の k-mer が桁違いに多い (高カバレッジ) 場合は、集合がエラーに埋め尽くされないところまでカットオフを上げること
        /// </summary>
        [Fact]
        public void V_低頻度エラーが集合を支配する場合はカットオフを上げる()
        {
            const int l_truePeak = 50;
            const long l_trueGenomeSize = 6_000_000L;
            var l_histogram = V_構築_現実的スペクトル(l_truePeak, l_trueGenomeSize, p_エラー係数: 60_000_000L);
            l_histogram[2UL] = 9_000_000L;
            var l_suggestion = KmerHistogram.Get_推奨カットオフ(l_histogram);
            Assert.NotNull(l_suggestion);
            Assert.True(l_suggestion > 2, $"cutoff should have been raised above 2, but was {l_suggestion}");
            var l_analysis = KmerHistogram.Get_解析結果(l_histogram);
            Assert.NotNull(l_analysis);
            Assert.True(l_suggestion <= l_analysis.A_谷, $"cutoff {l_suggestion} exceeded the valley {l_analysis.A_谷}");
        }

        /// <summary>
        /// 単調減少するヒストグラムでは null を返す
        /// </summary>
        [Fact]
        public void V_単調減少するヒストグラムではnullを返す()
        {
            Dictionary<ulong, long> l_histogram = new()
            {
                [1UL] = 100L,
                [2UL] = 50L,
                [3UL] = 10L,
                [4UL] = 1L,
            };
            Assert.Null(KmerHistogram.Get_推奨カットオフ(l_histogram));
        }

        /// <summary>
        /// 出現回数が 2 までしか無いヒストグラムはスペクトルとして成立しておらず、谷も山も判定できない
        /// </summary>
        [Fact]
        public void V_退化した2区分のヒストグラムではnullを返す()
        {
            Dictionary<ulong, long> l_histogram = new()
            {
                [1UL] = 5L,
                [2UL] = 500L,
            };
            Assert.Null(KmerHistogram.Get_推奨カットオフ(l_histogram));
        }

        /// <summary>
        /// 解析上の谷が 1 でも、推奨値としては 2 を下回らないこと
        /// </summary>
        [Fact]
        public void V_出現回数1に谷がある場合は下限まで引き上げる()
        {
            Dictionary<ulong, long> l_histogram = new()
            {
                [1UL] = 100L,
                [2UL] = 300L,
                [3UL] = 900L,
                [4UL] = 2_000L,
                [5UL] = 4_000L,
                [6UL] = 5_000L,
                [7UL] = 3_000L,
            };
            var l_analysis = KmerHistogram.Get_解析結果(l_histogram);
            Assert.NotNull(l_analysis);
            Assert.Equal(1UL, l_analysis.A_谷);
            Assert.Equal(KmerHistogram.C_推奨カットオフの下限, KmerHistogram.Get_推奨カットオフ(l_histogram));
        }

        /// <summary>
        /// エラー由来の裾と単一コピーの山が重なった、実データに近い連続的なスペクトル
        /// </summary>
        [Fact]
        public void V_連続的な二峰性スペクトルから谷_山_ゲノムサイズを報告する()
        {
            const int l_truePeak = 30;
            const long l_trueGenomeSize = 6_000_000L;
            var l_histogram = V_構築_現実的スペクトル(l_truePeak, l_trueGenomeSize, p_エラー係数: 10_000_000L);
            var l_analysis = KmerHistogram.Get_解析結果(l_histogram);
            Assert.NotNull(l_analysis);
            Assert.InRange(l_analysis.A_谷, 8UL, 22UL);
            Assert.InRange(l_analysis.A_ピーク出現回数, 27UL, 33UL);
            Assert.InRange(l_analysis.A_推定ゲノムサイズ, (long)(l_trueGenomeSize * 0.8D), (long)(l_trueGenomeSize * 1.3D));
        }

        /// <summary>
        /// アダプタ配列やコンタミ由来の、桁違いに出現回数の多い k-mer が混ざっていてもゲノムサイズ推定が壊れないこと
        /// </summary>
        [Fact]
        public void V_極端な外れ値があってもゲノムサイズ推定は膨れない()
        {
            const int l_truePeak = 30;
            const long l_trueGenomeSize = 6_000_000L;
            var l_histogram = V_構築_現実的スペクトル(l_truePeak, l_trueGenomeSize, p_エラー係数: 10_000_000L);
            var l_baseline = KmerHistogram.Get_解析結果(l_histogram);
            l_histogram[1_000_000UL] = 50L;
            var l_withOutliers = KmerHistogram.Get_解析結果(l_histogram);
            Assert.NotNull(l_baseline);
            Assert.NotNull(l_withOutliers);
            Assert.Equal(l_baseline.A_推定ゲノムサイズ, l_withOutliers.A_推定ゲノムサイズ);
        }

        /// <summary>
        /// 要約は、ヒストグラム自身が持つ最大キーで止まる
        /// </summary>
        [Fact]
        public void V_要約はヒストグラム自身の最大キーで止まる()
        {
            Dictionary<ulong, long> l_histogram = new()
            {
                [1UL] = 10L,
                [3UL] = 5L,
            };
            var l_summary = KmerHistogram.Get_要約(l_histogram, p_表示上限: 10UL);
            Assert.Equal(C_項目_1_10__2_0__3_5, l_summary);
        }

        /// <summary>
        /// 要約は、ヒストグラムがさらに続いていても表示上限で切り詰める
        /// </summary>
        [Fact]
        public void V_要約は表示上限を超えると切り詰める()
        {
            Dictionary<ulong, long> l_histogram = new()
            {
                [1UL] = 10L,
                [2UL] = 5L,
                [3UL] = 2L,
                [4UL] = 1L,
            };
            var l_summary = KmerHistogram.Get_要約(l_histogram, p_表示上限: 2UL);
            Assert.Equal(C_項目_1_10__2_5, l_summary);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// エラー由来の裾 (出現回数の二乗に反比例して減衰) と、単一コピーの山 (平均 p_ピーク の正規分布状) を重ね合わせた、実データに近い形のヒストグラムを作る
        /// </summary>
        /// <param name="p_ピーク"></param>
        /// <param name="p_ゲノムサイズ"></param>
        /// <param name="p_エラー係数"></param>
        /// <returns></returns>
        private static Dictionary<ulong, long> V_構築_現実的スペクトル(int p_ピーク, long p_ゲノムサイズ, long p_エラー係数)
        {
            var l_histogram = new Dictionary<ulong, long>();
            var l_sd = Math.Sqrt(p_ピーク);
            for (var l_count = 1UL; l_count <= (ulong)(p_ピーク * 2); l_count++)
            {
                var l_error = (long)(p_エラー係数 / (double)(l_count * l_count));
                var l_genome = (long)(p_ゲノムサイズ * Math.Exp(-Math.Pow((double)l_count - p_ピーク, 2D) / (2D * l_sd * l_sd)) / (l_sd * Math.Sqrt(2D * Math.PI)));
                l_histogram[l_count] = l_error + l_genome;
            }

            return l_histogram;
        }

        #endregion
    }
}
