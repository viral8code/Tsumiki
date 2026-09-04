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
        public void SuggestCutoff_ClassicBimodalSpectrum_FindsValley()
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

            var suggestion = KmerHistogram.Get_推奨カットオフ(histogram);

            Assert.Equal(3UL, suggestion);
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

        [Fact]
        public void SuggestCutoff_RisesImmediatelyAtCountTwo_ReturnsOne()
        {
            Dictionary<ulong, long> histogram = new()
            {
                [1] = 5,
                [2] = 500,
            };

            Assert.Equal(1UL, KmerHistogram.Get_推奨カットオフ(histogram));
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
