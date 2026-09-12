using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// 壊れた FASTQ で範囲外アクセスや無限ループにせず、どこが不正かを言って止まることを固定する
    /// </summary>
    public class FastqReaderTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public FastqReaderTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_fastq_reader_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 正しく整形されたレコードを読み込めることを検証する
        /// </summary>
        [Fact]
        public void Get_次のリード_正しく整形されたレコードを読み込む()
        {
            var l_パス = this.Get_書き出し先("@r1\nACGTACGT\n+\nIIIIIIII\n");
            using var l_読み込み = new FastqReader(l_パス);

            Assert.True(l_読み込み.Get_続きがあるか());
            var l_リード = l_読み込み.Get_次のリード_軽量();
            Assert.Equal("@r1", l_リード.A_ID);
            Assert.Equal("ACGTACGT", l_リード.A_生リード);
        }

        /// <summary>
        /// クオリティ長が配列長と異なる場合は例外を投げることを検証する
        /// </summary>
        [Fact]
        public void Get_次のリード_クオリティ長が配列長と異なる場合は例外を投げる()
        {
            var l_パス = this.Get_書き出し先("@r1\nACGTACGT\n+\nIIII\n");
            using var l_読み込み = new FastqReader(l_パス);

            var l_例外 = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード_軽量());
            Assert.Contains("@r1", l_例外.Message);
        }

        /// <summary>
        /// レコードの途中でファイルが途切れている場合は例外を投げることを検証する
        /// </summary>
        [Fact]
        public void Get_次のリード_レコードの途中でファイルが途切れている場合は例外を投げる()
        {
            var l_パス = this.Get_書き出し先("@r1\nACGTACGT\n+\n");
            using var l_読み込み = new FastqReader(l_パス);

            _ = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード_軽量());
        }

        /// <summary>
        /// 曖昧塩基を扱う経路でも同様に長さの不一致を検査することを検証する
        /// </summary>
        [Fact]
        public void Get_次のリード_曖昧塩基を扱う経路でも同様に検査する()
        {
            var l_パス = this.Get_書き出し先("@r1\nACGTACGT\n+\nIIIII\n");
            using var l_読み込み = new FastqReader(l_パス);

            _ = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード());
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 中身を書き出した一時ファイルのパスを返す
        /// </summary>
        /// <param name="p_内容">書き出す中身</param>
        /// <returns>書き出したパス</returns>
        private string Get_書き出し先(string p_内容)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, Guid.NewGuid().ToString("N") + ".fastq");
            File.WriteAllText(l_パス, p_内容);
            return l_パス;
        }

        #endregion

    }
}
