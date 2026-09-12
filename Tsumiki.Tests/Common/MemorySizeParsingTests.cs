using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// -mem のサイズ指定 ("2G" など) の解釈を固定する
    /// </summary>
    /// <remarks>
    /// 接尾辞は 2 進接頭辞 (1 K = 1024) <br/>
    /// メモリ量の指定なので 1000 刻みより 1024 刻みのほうが直感に合う<br/>
    /// 接尾辞が無い場合は MB とみなす
    /// </remarks>
    public class MemorySizeParsingTests
    {
        #region 公開メソッド

        /// <summary>
        /// 接尾辞付きのサイズ指定を読み取る
        /// </summary>
        /// <param name="p_文字列"></param>
        /// <param name="p_期待値"></param>
        [Theory]
        [InlineData("2G", 2L * 1_024L * 1_024L * 1_024L)]
        [InlineData("2g", 2L * 1_024L * 1_024L * 1_024L)]
        [InlineData("2GB", 2L * 1_024L * 1_024L * 1_024L)]
        [InlineData("512M", 512L * 1_024L * 1_024L)]
        [InlineData("512m", 512L * 1_024L * 1_024L)]
        [InlineData("1024K", 1_024L * 1_024L)]
        [InlineData("1T", 1_024L * 1_024L * 1_024L * 1_024L)]
        [InlineData("1.5G", (long)(1.5D * 1_024D * 1_024D * 1_024D))]
        [InlineData("0.5G", 512L * 1_024L * 1_024L)]
        public void V_接尾辞付きサイズを読み取る(string p_文字列, long p_期待値)
        {
            Assert.Equal(p_期待値, Util.V_変換_メモリサイズ(p_文字列));
        }

        /// <summary>
        /// 接尾辞なしは MB
        /// </summary>
        /// <remarks>
        /// 単なる数値で指定したときに「バイト」と解釈すると現実的にありえない小ささになるため
        /// </remarks>
        /// <param name="p_文字列"></param>
        /// <param name="p_期待値"></param>
        [Theory]
        [InlineData("768", 768L * 1_024L * 1_024L)]
        [InlineData("2048", 2_048L * 1_024L * 1_024L)]
        public void V_接尾辞なしはメガバイト扱い(string p_文字列, long p_期待値)
        {
            Assert.Equal(p_期待値, Util.V_変換_メモリサイズ(p_文字列));
        }

        /// <summary>
        /// 不正な入力を拒否する
        /// </summary>
        /// <param name="p_文字列"></param>
        [Theory]
        [InlineData("")]
        [InlineData("   ")]
        [InlineData("abc")]
        [InlineData("-1G")]
        [InlineData("0")]
        [InlineData("G")]
        public void V_不正な入力を拒否する(string p_文字列)
        {
            _ = Assert.Throws<ArgumentException>(() => Util.V_変換_メモリサイズ(p_文字列));
        }

        /// <summary>
        /// Parameters が接尾辞付きのメモリ予算を受け付ける
        /// </summary>
        [Fact]
        public void V_メモリ予算に接尾辞付き指定を受け付ける()
        {
            var l_パラメータ = new Parameters { A_メモリ予算 = "2G" };
            Assert.Equal(2L * 1_024L * 1_024L * 1_024L, l_パラメータ.A_メモリ予算バイト数);
        }

        /// <summary>
        /// 表示は読みやすい形に戻ること (パラメータの一覧表示で使う)
        /// </summary>
        /// <param name="p_バイト数"></param>
        /// <param name="p_期待値"></param>
        [Theory]
        [InlineData(2L * 1_024L * 1_024L * 1_024L, "2 GB")]
        [InlineData(768L * 1_024L * 1_024L, "768 MB")]
        [InlineData(1_536L * 1_024L * 1_024L, "1.5 GB")]
        public void V_読みやすい形式へ変換する(long p_バイト数, string p_期待値)
        {
            Assert.Equal(p_期待値, Util.Get_表示用メモリサイズ(p_バイト数));
        }

        #endregion

    }
}
