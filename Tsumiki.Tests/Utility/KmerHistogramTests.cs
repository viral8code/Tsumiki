using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    public class KmerHistogramTests
    {
        [Fact]
        public void SuggestCutoff_EmptyHistogram_ReturnsNull()
        {
            Assert.Null(KmerHistogram.Get_推奨カットオフ(new Dictionary<ulong, long>()));
        }

        [Fact]
        public void Analyse_ClassicBimodalSpectrum_FindsValley()
        {
            // エラー由来の山(count=1,2)、谷(count=3)、真のゲノム由来の山(count~30)。
            Dictionary<ulong, long> histogram = new()
            {
                [1] = 10_000,
                [2] = 3_000,
                [3] = 500,
                [4] = 800,
                [5] = 5_000,
                [30] = 20_000,
                [31] = 19_000,
            };

            var analysis = KmerHistogram.Get_解析結果(histogram);

            Assert.NotNull(analysis);
            Assert.Equal(3UL, analysis.A_谷);
            Assert.Equal(30UL, analysis.A_ピーク出現回数);
        }

        /// <summary>
        /// 推奨カットオフは谷そのものではない。谷はエラー由来の曲線と
        /// ゲノム由来の曲線が交わる点なので、そこで切るとゲノム側の左裾を
        /// 削ってしまう。エラーが集合を支配しない範囲でできるだけ低く返す。
        ///
        /// 上のスペクトルなら、出現回数2以上を残せば 48,300 種類で
        /// 推定ゲノムサイズ 40,623 の 1.19 倍に収まるため、谷(3)まで
        /// 上げる必要はない。
        /// </summary>
        [Fact]
        public void SuggestCutoff_StopsBelowTheValley_WhenErrorsAlreadyDoNotDominate()
        {
            Dictionary<ulong, long> histogram = new()
            {
                [1] = 10_000,
                [2] = 3_000,
                [3] = 500,
                [4] = 800,
                [5] = 5_000,
                [30] = 20_000,
                [31] = 19_000,
            };

            Assert.Equal(2UL, KmerHistogram.Get_推奨カットオフ(histogram));
        }

        /// <summary>
        /// エラー由来の k-mer が桁違いに多い(高カバレッジ)場合は、
        /// 集合がエラーに埋め尽くされないところまでカットオフを上げること。
        /// 品質は変わらないがメモリが減る。
        /// </summary>
        [Fact]
        public void SuggestCutoff_RaisesTheCutoff_WhenLowCountErrorsDominateTheSet()
        {
            const int truePeak = 50;
            const long trueGenomeSize = 6_000_000;
            var histogram = BuildRealisticSpectrum(truePeak, trueGenomeSize, p_エラー係数: 60_000_000);
            // 出現回数2のエラー k-mer を、ゲノムの種類数を超える規模で載せる。
            histogram[2] = 9_000_000;

            var suggestion = KmerHistogram.Get_推奨カットオフ(histogram);

            Assert.NotNull(suggestion);
            Assert.True(suggestion > 2, $"cutoff should have been raised above 2, but was {suggestion}");
            // 谷を超えて上げてはいけない(そこから先はゲノム由来しか残っていない)。
            var analysis = KmerHistogram.Get_解析結果(histogram);
            Assert.NotNull(analysis);
            Assert.True(suggestion <= analysis.A_谷, $"cutoff {suggestion} exceeded the valley {analysis.A_谷}");
        }

        [Fact]
        public void SuggestCutoff_MonotonicDecrease_ReturnsNull()
        {
            Dictionary<ulong, long> histogram = new()
            {
                [1] = 100,
                [2] = 50,
                [3] = 10,
                [4] = 1,
            };

            Assert.Null(KmerHistogram.Get_推奨カットオフ(histogram));
        }

        /// <summary>
        /// 出現回数が2までしか無いヒストグラムはスペクトルとして成立しておらず、
        /// 谷も山も判定できない。推測で値を返すより、判定不能を返して
        /// 既定値を維持させるほうが安全。
        /// </summary>
        [Fact]
        public void SuggestCutoff_DegenerateTwoBucketHistogram_ReturnsNull()
        {
            Dictionary<ulong, long> histogram = new()
            {
                [1] = 5,
                [2] = 500,
            };

            Assert.Null(KmerHistogram.Get_推奨カットオフ(histogram));
        }

        /// <summary>
        /// 解析上の谷が1でも、推奨値としては2を下回らないこと。
        /// 出現回数1の k-mer はどのカバレッジ帯でもほぼ全てエラー由来であり、
        /// 残すとメモリを食ったうえでグラフが偽の枝だらけになる。
        /// </summary>
        [Fact]
        public void SuggestCutoff_ValleyAtCountOne_IsRaisedToTheFloor()
        {
            // count=1 が最小(エラーがほとんど無いデータ)で、そこから
            // 単一コピーの山へ立ち上がるスペクトル。
            Dictionary<ulong, long> histogram = new()
            {
                [1] = 100,
                [2] = 300,
                [3] = 900,
                [4] = 2_000,
                [5] = 4_000,
                [6] = 5_000,
                [7] = 3_000,
            };

            var analysis = KmerHistogram.Get_解析結果(histogram);

            Assert.NotNull(analysis);
            Assert.Equal(1UL, analysis.A_谷);
            Assert.Equal(KmerHistogram.推奨カットオフの下限, KmerHistogram.Get_推奨カットオフ(histogram));
        }

        /// <summary>
        /// エラー由来の裾と単一コピーの山が重なった、実データに近い連続的な
        /// スペクトル。谷・山・ゲノムサイズがまとめて取れること。
        ///
        /// 素朴な「最初に頻度が増えた位置」だけを見る実装は、谷の底が平らな
        /// 実データでノイズに引きずられて答えがぶれた(同じ検体の 100x で
        /// 6 と 11 の両方が出た)。底の最小値を取り直すことで安定させている。
        /// </summary>
        [Fact]
        public void Analyse_ContinuousBimodalSpectrum_ReportsValleyPeakAndGenomeSize()
        {
            const int truePeak = 30;
            const long trueGenomeSize = 6_000_000;
            var histogram = BuildRealisticSpectrum(truePeak, trueGenomeSize, p_エラー係数: 10_000_000);

            var analysis = KmerHistogram.Get_解析結果(histogram);

            Assert.NotNull(analysis);
            // 谷はエラーの裾とゲノムの山の交点付近に来る。
            Assert.InRange(analysis.A_谷, 8UL, 22UL);
            Assert.InRange(analysis.A_ピーク出現回数, 27UL, 33UL);
            // カットオフを超えて残るエラー k-mer のぶんだけ上振れするが、
            // 真の値のオーダーは取れていなければならない。
            Assert.InRange(analysis.A_推定ゲノムサイズ, (long)(trueGenomeSize * 0.8), (long)(trueGenomeSize * 1.3));
        }

        /// <summary>
        /// アダプタ配列やコンタミ由来の、桁違いに出現回数の多い k-mer が
        /// 混ざっていてもゲノムサイズ推定が壊れないこと。素直に延べ数へ
        /// 足し込むと、たった数十種類でゲノムサイズが何倍にも膨れる。
        /// </summary>
        [Fact]
        public void Analyse_ExtremeOutlierCounts_DoNotInflateTheGenomeSizeEstimate()
        {
            const int truePeak = 30;
            const long trueGenomeSize = 6_000_000;
            var histogram = BuildRealisticSpectrum(truePeak, trueGenomeSize, p_エラー係数: 10_000_000);
            var baseline = KmerHistogram.Get_解析結果(histogram);

            // 100万回出現する k-mer を50種類混ぜる(延べ 5000 万)。
            histogram[1_000_000] = 50;
            var withOutliers = KmerHistogram.Get_解析結果(histogram);

            Assert.NotNull(baseline);
            Assert.NotNull(withOutliers);
            Assert.Equal(baseline.A_推定ゲノムサイズ, withOutliers.A_推定ゲノムサイズ);
        }

        /// <summary>
        /// エラー由来の裾(出現回数の二乗に反比例して減衰)と、単一コピーの
        /// 山(平均 p_ピーク の正規分布状)を重ね合わせた、実データに近い形の
        /// ヒストグラムを作る。
        /// </summary>
        private static Dictionary<ulong, long> BuildRealisticSpectrum(int p_ピーク, long p_ゲノムサイズ, long p_エラー係数)
        {
            var histogram = new Dictionary<ulong, long>();
            var sd = Math.Sqrt(p_ピーク);
            for (var count = 1UL; count <= (ulong)(p_ピーク * 2); count++)
            {
                var error = (long)(p_エラー係数 / (double)(count * count));
                var genome = (long)(p_ゲノムサイズ
                    * Math.Exp(-Math.Pow((double)count - p_ピーク, 2) / (2 * sd * sd))
                    / (sd * Math.Sqrt(2 * Math.PI)));
                histogram[count] = error + genome;
            }
            return histogram;
        }

        [Fact]
        public void FormatSummary_StopsAtHistogramsOwnMaxKey_EvenIfMaxCountAllowsMore()
        {
            Dictionary<ulong, long> histogram = new()
            {
                [1] = 10,
                [3] = 5,
            };

            var summary = KmerHistogram.Get_要約(histogram, p_表示上限: 10);

            Assert.Equal("1:10, 2:0, 3:5", summary);
        }

        [Fact]
        public void FormatSummary_TruncatesAtMaxCount_WhenHistogramExtendsFurther()
        {
            Dictionary<ulong, long> histogram = new()
            {
                [1] = 10,
                [2] = 5,
                [3] = 2,
                [4] = 1,
            };

            var summary = KmerHistogram.Get_要約(histogram, p_表示上限: 2);

            Assert.Equal("1:10, 2:5", summary);
        }
    }
}
