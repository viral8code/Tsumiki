using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// FASTQ の書き出し (非圧縮と gzip 圧縮) の検証
    /// </summary>
    public class FastqWriterTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki fastq writer tests
        /// </summary>
        private const string C_項目_tsumiki_fastq_writer_tests = "tsumiki_fastq_writer_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// gzip の先頭 1 バイト目
        /// </summary>
        private const byte C_gzip先頭1 = 0x1F;

        /// <summary>
        /// gzip の先頭 2 バイト目
        /// </summary>
        private const byte C_gzip先頭2 = 0x8B;

        /// <summary>
        /// FASTQ のヘッダ行の先頭文字
        /// </summary>
        private const char C_ヘッダ記号 = '@';

        #endregion

        #region 内部変数

        /// <summary>
        /// 作業用ディレクトリ
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の作業ディレクトリを用意する
        /// </summary>
        public FastqWriterTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_fastq_writer_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 圧縮して書いた 3 リードは、先頭が gzip の印で、FastqReader で読むと同じ内容になる
        /// </summary>
        [Fact]
        public void V_圧縮して書くと先頭がgzipの印で_読むと同じ内容になる()
        {
            var l_パス = this.Get_パス("out.fq.gz");
            var l_書いた = new[]
            {
                ("@r0", "ACGT", "IIII"),
                ("@r1", "GGCA", "FFFF"),
                ("@r2", "TTAA", "!!!!"),
            };
            using (var l_書込 = new FastqWriter(l_パス, p_Is圧縮: true))
            {
                foreach (var (l_ID, l_配列, l_クオリティ) in l_書いた)
                {
                    l_書込.V_書き込み(l_ID, l_配列, l_クオリティ);
                }
            }

            var l_先頭 = File.ReadAllBytes(l_パス).Take(2).ToArray();
            Assert.Equal([C_gzip先頭1, C_gzip先頭2], l_先頭);

            using var l_読込 = new FastqReader(l_パス);
            var l_読んだ = new List<(string A_ID, string A_配列, string A_クオリティ)>();
            while (l_読込.Has続き())
            {
                l_読んだ.Add(l_読込.Get_次のレコード());
            }

            Assert.Equal(l_書いた, l_読んだ);
        }

        /// <summary>
        /// 圧縮の指定を省くと、従来どおり非圧縮の平文で、先頭は '@' になる
        /// </summary>
        [Fact]
        public void V_圧縮を省くと非圧縮で先頭はアットマーク()
        {
            var l_パス = this.Get_パス("plain.fq");
            using (var l_書込 = new FastqWriter(l_パス))
            {
                l_書込.V_書き込み("@r0", "ACGT", "IIII");
            }

            var l_バイト = File.ReadAllBytes(l_パス);
            Assert.Equal((byte)C_ヘッダ記号, l_バイト[0]);
            var l_期待 = string.Join(Environment.NewLine, "@r0", "ACGT", "+", "IIII") + Environment.NewLine;
            Assert.Equal(l_期待, File.ReadAllText(l_パス));
        }

        /// <summary>
        /// 作業ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }

            GC.SuppressFinalize(this);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 作業ディレクトリ内のパスを返す
        /// </summary>
        /// <param name="p_名前">ファイル名</param>
        /// <returns>パス</returns>
        private string Get_パス(string p_名前)
        {
            return Path.Combine(this._作業ディレクトリ, p_名前);
        }

        #endregion
    }
}
