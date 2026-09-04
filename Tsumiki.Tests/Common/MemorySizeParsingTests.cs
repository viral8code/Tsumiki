using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// -mem のサイズ指定("2G" など)の解釈を固定する。
    /// 接尾辞は2進接頭辞(1K = 1024)。メモリ量の指定なので 1000 刻みより
    /// 1024 刻みのほうが直感に合う。接尾辞が無い場合は MB とみなす。
    /// </summary>
    public class MemorySizeParsingTests
    {
        [Theory]
        [InlineData("2G", 2L * 1024 * 1024 * 1024)]
        [InlineData("2g", 2L * 1024 * 1024 * 1024)]
        [InlineData("2GB", 2L * 1024 * 1024 * 1024)]
        [InlineData("512M", 512L * 1024 * 1024)]
        [InlineData("512m", 512L * 1024 * 1024)]
        [InlineData("1024K", 1024L * 1024)]
        [InlineData("1T", 1024L * 1024 * 1024 * 1024)]
        [InlineData("1.5G", (long)(1.5 * 1024 * 1024 * 1024))]
        [InlineData("0.5G", 512L * 1024 * 1024)]
        public void ParseMemorySize_ReadsSuffixedSizes(string text, long expected)
        {
            Assert.Equal(expected, Util.V_変換_メモリサイズ(text));
        }

        /// <summary>
        /// 接尾辞なしは MB。単なる数値で指定したときに「バイト」と解釈すると
        /// 現実的にありえない小ささになるため。
        /// </summary>
        [Theory]
        [InlineData("768", 768L * 1024 * 1024)]
        [InlineData("2048", 2048L * 1024 * 1024)]
        public void ParseMemorySize_BareNumberMeansMegabytes(string text, long expected)
        {
            Assert.Equal(expected, Util.V_変換_メモリサイズ(text));
        }

        [Theory]
        [InlineData("")]
        [InlineData("   ")]
        [InlineData("abc")]
        [InlineData("-1G")]
        [InlineData("0")]
        [InlineData("G")]
        public void ParseMemorySize_RejectsInvalidInput(string text)
        {
            _ = Assert.Throws<ArgumentException>(() => Util.V_変換_メモリサイズ(text));
        }

        [Fact]
        public void Parameters_AcceptsSuffixedMemoryBudget()
        {
            var param = new Parameters { A_メモリ予算 = "2G" };
            Assert.Equal(2L * 1024 * 1024 * 1024, param.A_メモリ予算バイト数);
        }

        /// <summary>
        /// 表示は読みやすい形に戻ること(パラメータの一覧表示で使う)。
        /// </summary>
        [Theory]
        [InlineData(2L * 1024 * 1024 * 1024, "2 GB")]
        [InlineData(768L * 1024 * 1024, "768 MB")]
        [InlineData(1536L * 1024 * 1024, "1.5 GB")]
        public void FormatMemorySize_RoundTripsToAReadableForm(long bytes, string expected)
        {
            Assert.Equal(expected, Util.Get_表示用メモリサイズ(bytes));
        }
    }
}
