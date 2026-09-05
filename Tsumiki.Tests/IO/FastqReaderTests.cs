using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// 壊れた FASTQ で範囲外アクセスや無限ループにせず、
    /// どこが不正かを言って止まることを固定する。
    /// </summary>
    public class FastqReaderTests : IDisposable
    {
        private readonly string _tempDir;

        public FastqReaderTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_fastq_reader_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private string Get_書き出し先(string p_内容)
        {
            var l_パス = Path.Combine(this._tempDir, Guid.NewGuid().ToString("N") + ".fastq");
            File.WriteAllText(l_パス, p_内容);
            return l_パス;
        }

        [Fact]
        public void Get_次のリード_ReadsAWellFormedRecord()
        {
            var l_パス = this.Get_書き出し先("@r1\nACGTACGT\n+\nIIIIIIII\n");
            using var l_読み込み = new FastqReader(l_パス);

            Assert.True(l_読み込み.Get_続きがあるか());
            var l_リード = l_読み込み.Get_次のリード_軽量();
            Assert.Equal("@r1", l_リード.A_ID);
            Assert.Equal("ACGTACGT", l_リード.A_生リード);
        }

        [Fact]
        public void Get_次のリード_ThrowsWhenQualityLengthDiffersFromSequence()
        {
            var l_パス = this.Get_書き出し先("@r1\nACGTACGT\n+\nIIII\n");
            using var l_読み込み = new FastqReader(l_パス);

            var l_例外 = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード_軽量());
            Assert.Contains("@r1", l_例外.Message);
        }

        [Fact]
        public void Get_次のリード_ThrowsWhenTheFileIsTruncatedMidRecord()
        {
            var l_パス = this.Get_書き出し先("@r1\nACGTACGT\n+\n");
            using var l_読み込み = new FastqReader(l_パス);

            _ = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード_軽量());
        }

        [Fact]
        public void Get_次のリード_ChecksTheAmbiguousBasePathToo()
        {
            var l_パス = this.Get_書き出し先("@r1\nACGTACGT\n+\nIIIII\n");
            using var l_読み込み = new FastqReader(l_パス);

            _ = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード());
        }
    }
}
