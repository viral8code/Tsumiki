using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// リードファイルからの代表リード長の抽出。k 長の自動選択の入力になる。
    /// </summary>
    public class ReadLengthSnifferTests : IDisposable
    {
        private readonly string _tempDir;

        public ReadLengthSnifferTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_readlength_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private string WriteFastq(string name, params int[] readLengths)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new StreamWriter(path);
            var id = 1;
            foreach (var length in readLengths)
            {
                writer.WriteLine($"@read{id++}");
                writer.WriteLine(new string('A', length));
                writer.WriteLine("+");
                writer.WriteLine(new string('I', length));
            }
            return path;
        }

        [Fact]
        public void GetReadLength_UniformReads_ReturnsThatLength()
        {
            var path = this.WriteFastq("uniform.fq", 150, 150, 150, 150);

            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(path));
        }

        /// <summary>
        /// トリミング済みのデータではリード長がばらつく。平均や最大ではなく
        /// 中央値を使うことで、少数の極端に短いリードに引きずられない。
        /// </summary>
        [Fact]
        public void GetReadLength_TrimmedReads_ReturnsTheMedianRatherThanTheMeanOrMax()
        {
            var path = this.WriteFastq("trimmed.fq", 35, 40, 148, 150, 150, 150, 151);

            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(path));
        }

        [Fact]
        public void GetReadLength_StopsAfterTheSampleLimit()
        {
            // 先頭2本だけを見れば 200 が中央値になる。ファイル全体を見ると 50。
            var lengths = new[] { 200, 200 }.Concat(Enumerable.Repeat(50, 100)).ToArray();
            var path = this.WriteFastq("limited.fq", lengths);

            Assert.Equal(200, ReadLengthSniffer.Get_代表リード長(path, p_標本上限: 2));
        }

        [Fact]
        public void GetReadLength_EmptyFile_ReturnsNull()
        {
            var path = Path.Combine(this._tempDir, "empty.fq");
            File.WriteAllText(path, string.Empty);

            Assert.Null(ReadLengthSniffer.Get_代表リード長(path));
        }

        /// <summary>
        /// k は「どちらのリードからも k-mer が取れる」必要があるため、
        /// ペアで長さが違う場合は短いほうに合わせる。
        /// </summary>
        [Fact]
        public void GetReadLength_PairedFilesWithDifferentLengths_ReturnsTheShorterOne()
        {
            var path1 = this.WriteFastq("pair.1.fq", 150, 150, 150);
            var path2 = this.WriteFastq("pair.2.fq", 100, 100, 100);

            Assert.Equal(100, ReadLengthSniffer.Get_代表リード長(path1, path2));
        }

        [Fact]
        public void GetReadLength_SingleEnd_IgnoresTheMissingSecondFile()
        {
            var path1 = this.WriteFastq("single.fq", 150, 150, 150);

            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(path1, null));
        }
    }
}
