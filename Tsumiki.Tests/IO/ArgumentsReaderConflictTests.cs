using Tsumiki.Commons;
using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// 同時に指定すると意味が通らない組を、始める前に止めることの検証
    /// </summary>
    public class ArgumentsReaderConflictTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 300
        /// </summary>
        private const string C_項目_300 = "300";

        #endregion

        #region 内部変数

        /// <summary>
        /// 実在するだけのダミーのリードのパス
        /// </summary>
        private readonly string[] _ダミーリードパス = [Path.GetTempFileName(), Path.GetTempFileName(), Path.GetTempFileName(), Path.GetTempFileName()];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 片付ける
        /// </summary>
        public void Dispose()
        {
            foreach (var l_パス in this._ダミーリードパス)
            {
                File.Delete(l_パス);
            }
        }

        /// <summary>
        /// 相反する組は理由付きで止め、片方だけなら通す
        /// </summary>
        /// <param name="p_追加">入力に足す引数</param>
        /// <param name="p_Is止める">止めるべきか</param>
        /// <param name="p_理由に含む">止めるときの理由に含まれるはずのキー</param>
        [Theory]
        [InlineData(new[] { "-inmem", "-mem", "4G" }, true, "-mem")]
        [InlineData(new[] { "-mem", "4G", "-inmem" }, true, "-mem")]
        [InlineData(new[] { "-inmem", "-rs" }, true, "-rs")]
        [InlineData(new[] { "-inmem" }, false, "")]
        [InlineData(new[] { "-inmem", "-prep" }, true, "-prep")]
        [InlineData(new[] { "-prep", "-rs" }, false, "")]
        [InlineData(new[] { "-mem", "4G", "-rs" }, false, "")]
        [InlineData(new[] { "-i", "300" }, false, "")]
        public void V_相反する組だけを止める(string[] p_追加, bool p_Is止める, string p_理由に含む)
        {
            string[] l_引数 = [Consts.引数キー.リード1のパス, this._ダミーリードパス[0], Consts.引数キー.リード2のパス, this._ダミーリードパス[1], ..p_追加];
            if (!p_Is止める)
            {
                _ = ArgumentsReader.Get_実行時引数(l_引数);
                return;
            }

            var l_例外 = Assert.Throws<ArgumentException>(() => ArgumentsReader.Get_実行時引数(l_引数));
            Assert.Contains(Consts.引数キー.オンメモリ, l_例外.Message, StringComparison.Ordinal);
            Assert.Contains(p_理由に含む, l_例外.Message, StringComparison.Ordinal);
        }

        /// <summary>
        /// インサートサイズの指定は、ペアのライブラリが 2 つ以上あると止める
        /// </summary>
        [Fact]
        public void V_インサートサイズの指定は複数のペアライブラリと併用できない()
        {
            var l_順リード = $"{this._ダミーリードパス[0]},{this._ダミーリードパス[2]}";
            var l_逆リード = $"{this._ダミーリードパス[1]},{this._ダミーリードパス[3]}";
            var l_例外 = Assert.Throws<ArgumentException>(() => ArgumentsReader.Get_実行時引数([Consts.引数キー.リード1のパス, l_順リード, Consts.引数キー.リード2のパス, l_逆リード, Consts.引数キー.インサートサイズ, C_項目_300]));
            Assert.Contains(Consts.引数キー.インサートサイズ, l_例外.Message, StringComparison.Ordinal);
            _ = ArgumentsReader.Get_実行時引数([
                Consts.引数キー.リード1のパス,
                this._ダミーリードパス[0],
                Consts.引数キー.リード2のパス,
                this._ダミーリードパス[1],
                Consts.引数キー.シングルのパス,
                this._ダミーリードパス[2],
                Consts.引数キー.インサートサイズ,
                C_項目_300
            ]);
        }

        #endregion
    }
}
