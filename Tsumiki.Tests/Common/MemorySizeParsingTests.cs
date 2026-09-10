using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// -mem のサイズ指定 ("2G" など) の解釈を固定する
    /// </summary>
    /// <remarks>
    /// 接尾辞は 2 進接頭辞 (1 K = 1024)<br/>
    /// メモリ量の指定なので 1000 刻みより
    /// 1024 刻みのほうが直感に合う<br/>
    /// 接尾辞が無い場合は MB とみなす
    /// </remarks>
    public class MemorySizeParsingTests
    {
        /// <summary>
        /// 接尾辞付きのサイズ指定を読み取る
        /// </summary>
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
        public void 接尾辞付きサイズを読み取る(string p_文字列, long p_期待値)
        {
            Assert.Equal(p_期待値, Util.V_変換_メモリサイズ(p_文字列));
        }

        /// <summary>
        /// 接尾辞なしは MB
        /// </summary>
        /// <remarks>
        /// 単なる数値で指定したときに「バイト」と解釈すると
        /// 現実的にありえない小ささになるため
        /// </remarks>
        [Theory]
        [InlineData("768", 768L * 1024 * 1024)]
        [InlineData("2048", 2048L * 1024 * 1024)]
        public void 接尾辞なしはメガバイト扱い(string p_文字列, long p_期待値)
        {
            Assert.Equal(p_期待値, Util.V_変換_メモリサイズ(p_文字列));
        }

        /// <summary>
        /// 不正な入力を拒否する
        /// </summary>
        [Theory]
        [InlineData("")]
        [InlineData("   ")]
        [InlineData("abc")]
        [InlineData("-1G")]
        [InlineData("0")]
        [InlineData("G")]
        public void 不正な入力を拒否する(string p_文字列)
        {
            _ = Assert.Throws<ArgumentException>(() => Util.V_変換_メモリサイズ(p_文字列));
        }

        /// <summary>
        /// Parameters が接尾辞付きのメモリ予算を受け付ける
        /// </summary>
        [Fact]
        public void メモリ予算に接尾辞付き指定を受け付ける()
        {
            var l_パラメータ = new Parameters { A_メモリ予算 = "2G" };
            Assert.Equal(2L * 1024 * 1024 * 1024, l_パラメータ.A_メモリ予算バイト数);
        }

        /// <summary>
        /// 表示は読みやすい形に戻ること (パラメータの一覧表示で使う)
        /// </summary>
        [Theory]
        [InlineData(2L * 1024 * 1024 * 1024, "2 GB")]
        [InlineData(768L * 1024 * 1024, "768 MB")]
        [InlineData(1536L * 1024 * 1024, "1.5 GB")]
        public void 読みやすい形式へ変換する(long p_バイト数, string p_期待値)
        {
            Assert.Equal(p_期待値, Util.Get_表示用メモリサイズ(p_バイト数));
        }
    }
}
