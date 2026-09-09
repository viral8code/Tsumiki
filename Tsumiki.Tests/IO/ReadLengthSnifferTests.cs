using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// リードファイルからの代表リード長の抽出
    /// </summary>
    /// <remarks>
    /// k 長の自動選択の入力になる
    /// </remarks>
    public class ReadLengthSnifferTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public ReadLengthSnifferTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_readlength_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        /// <summary>
        /// リードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_リード長群">書き出すリードの長さ</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_FASTQ(string p_ファイル名, params int[] p_リード長群)
        {
            var l_パス = Path.Combine(this._tempDir, p_ファイル名);
            using var l_writer = new StreamWriter(l_パス);
            var l_id = 1;
            foreach (var l_長さ in p_リード長群)
            {
                l_writer.WriteLine($"@read{l_id++}");
                l_writer.WriteLine(new string('A', l_長さ));
                l_writer.WriteLine("+");
                l_writer.WriteLine(new string('I', l_長さ));
            }
            return l_パス;
        }

        /// <summary>
        /// すべて同じ長さのリードでは、その長さを代表リード長として返すことを検証する
        /// </summary>
        [Fact]
        public void GetReadLength_すべて同じ長さのリードではその長さを返す()
        {
            var path = this.V_書き出し_FASTQ("uniform.fq", 150, 150, 150, 150);

            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(path));
        }

        /// <summary>
        /// トリミング済みのデータではリード長がばらつく
        /// </summary>
        /// <remarks>
        /// 平均や最大ではなく
        /// 中央値を使うことで、少数の極端に短いリードに引きずられない
        /// </remarks>
        [Fact]
        public void GetReadLength_トリミング済みデータでは平均や最大ではなく中央値を返す()
        {
            var path = this.V_書き出し_FASTQ("trimmed.fq", 35, 40, 148, 150, 150, 150, 151);

            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(path));
        }

        /// <summary>
        /// 標本上限に達したらそこで読み取りを打ち切ることを検証する
        /// </summary>
        [Fact]
        public void GetReadLength_標本上限に達したら打ち切る()
        {
            // 先頭 2 本だけを見れば 200 が中央値になる
            // ファイル全体を見ると 50
            var lengths = new[] { 200, 200 }.Concat(Enumerable.Repeat(50, 100)).ToArray();
            var path = this.V_書き出し_FASTQ("limited.fq", lengths);

            Assert.Equal(200, ReadLengthSniffer.Get_代表リード長(path, p_標本上限: 2));
        }

        /// <summary>
        /// 空のファイルでは null を返すことを検証する
        /// </summary>
        [Fact]
        public void GetReadLength_空ファイルではnullを返す()
        {
            var path = Path.Combine(this._tempDir, "empty.fq");
            File.WriteAllText(path, string.Empty);

            Assert.Null(ReadLengthSniffer.Get_代表リード長(path));
        }

        /// <summary>
        /// k は「どちらのリードからも k-mer が取れる」必要があるため、
        /// ペアで長さが違う場合は短いほうに合わせる
        /// </summary>
        [Fact]
        public void GetReadLength_ペアで長さが違う場合は短いほうを返す()
        {
            var path1 = this.V_書き出し_FASTQ("pair.1.fq", 150, 150, 150);
            var path2 = this.V_書き出し_FASTQ("pair.2.fq", 100, 100, 100);

            Assert.Equal(100, ReadLengthSniffer.Get_代表リード長(path1, path2));
        }

        /// <summary>
        /// シングルエンドでは、無い方のファイルを無視することを検証する
        /// </summary>
        [Fact]
        public void GetReadLength_シングルエンドでは無い方のファイルを無視する()
        {
            var path1 = this.V_書き出し_FASTQ("single.fq", 150, 150, 150);

            Assert.Equal(150, ReadLengthSniffer.Get_代表リード長(path1, null));
        }
    }
}
