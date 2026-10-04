using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// FASTA の読み込み (折り返した配列・空の配列・途中で終わるファイル) の確認
    /// </summary>
    public class FastaReaderTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 一時ファイルのパス
        /// </summary>
        private readonly string _パス = Path.GetTempFileName();

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 1 行の配列も、複数行に折り返した配列も、1 本につなげて読む
        /// </summary>
        [Fact]
        public void 折り返した配列もつなげて読む()
        {
            File.WriteAllText(this._パス, ">a desc\nACGT\n>b\nAAAA\nCCCC\r\nGG\n\n>c\nTTT\n\n\n");

            var l_全件 = FastaReader.Get_全エントリ(this._パス);

            Assert.Equal([("a desc", "ACGT"), ("b", "AAAACCCCGG"), ("c", "TTT")], l_全件);
        }

        /// <summary>
        /// 配列の無い見出しは空の配列として読む
        /// </summary>
        [Fact]
        public void 配列の無い見出しは空の配列()
        {
            File.WriteAllText(this._パス, ">a\n>b\nAC\n>c\n");

            Assert.Equal([("a", ""), ("b", "AC"), ("c", "")], FastaReader.Get_全エントリ(this._パス));
        }

        /// <summary>
        /// 空行しか無いファイルは、止まらずに例外にする
        /// </summary>
        [Fact]
        public void 空行しか無いファイルは止まらずに例外()
        {
            File.WriteAllText(this._パス, "\n\n");

            _ = Assert.Throws<InvalidDataException>(() => FastaReader.Get_全エントリ(this._パス));
        }

        /// <summary>
        /// 一時ファイルを片付ける
        /// </summary>
        public void Dispose()
        {
            File.Delete(this._パス);
            GC.SuppressFinalize(this);
        }

        #endregion
    }
}
