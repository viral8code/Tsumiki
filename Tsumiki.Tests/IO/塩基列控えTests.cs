using Tsumiki.IO;
using Tsumiki.Tests.Utility;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// 塩基列の控えから返す塩基列が、FASTQ をそのまま読んだときと同じであることの検証
    /// </summary>
    [Collection(中間データ置き場の集まり.名前)]
    public class 塩基列控えTests : IDisposable
    {
        private readonly string _作業ディレクトリ;

        public 塩基列控えTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_bases_cache_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        public void Dispose()
        {
            塩基列控え.V_全消去();
            塩基列控え.A_置き場所 = null;
            中間データ置き場.A_Is有効 = false;
            中間データ置き場.V_消去();
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        // 詰められるもの・詰められないもの (N、小文字、空、端数の長さ、長いもの) を混ぜる
        private static readonly string[] 塩基列群 =
        [
            "ACGT",
            "A",
            "ACGTA",
            "NNACGTNN",
            "acgtACGT",
            "",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            new string('G', 5000) + "C",
            "ACGTRYKM",
        ];

        private string V_書出_FASTQ(string p_名前, IEnumerable<string> p_配列群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_名前);
            using var l_書き込み = new StreamWriter(l_パス);
            var l_番号 = 0;
            foreach (var l_配列 in p_配列群)
            {
                l_書き込み.WriteLine($"@r{l_番号++}");
                l_書き込み.WriteLine(l_配列.Length == 0 ? "A" : l_配列);
                l_書き込み.WriteLine("+");
                l_書き込み.WriteLine(new string('I', l_配列.Length == 0 ? 1 : l_配列.Length));
            }
            return l_パス;
        }

        [Theory]
        [InlineData(false)]
        [InlineData(true)]
        public void 二回目は控えから同じ塩基列を返す(bool p_Isオンメモリ)
        {
            中間データ置き場.A_Is有効 = p_Isオンメモリ;
            塩基列控え.A_置き場所 = Path.Combine(this._作業ディレクトリ, "bases");
            var l_パス = this.V_書出_FASTQ("reads.fq", 塩基列群);
            var l_期待 = new FastqReaderの塩基列(l_パス).ToList();

            var l_一回目 = 塩基列控え.Get_塩基列(l_パス).ToList();
            Assert.True(塩基列控え.Has控え(l_パス));
            var l_二回目 = 塩基列控え.Get_塩基列(l_パス).ToList();

            Assert.Equal(l_期待, l_一回目);
            Assert.Equal(l_期待, l_二回目);
        }

        [Fact]
        public void 置き場所が無ければ控えを作らない()
        {
            var l_パス = this.V_書出_FASTQ("reads.fq", 塩基列群);

            var l_結果 = 塩基列控え.Get_塩基列(l_パス).ToList();

            Assert.Equal(new FastqReaderの塩基列(l_パス).ToList(), l_結果);
            Assert.False(塩基列控え.Has控え(l_パス));
        }

        [Fact]
        public void 途中で読むのをやめたら控えを登録しない()
        {
            塩基列控え.A_置き場所 = Path.Combine(this._作業ディレクトリ, "bases");
            var l_パス = this.V_書出_FASTQ("reads.fq", 塩基列群);

            _ = 塩基列控え.Get_塩基列(l_パス).Take(2).ToList();

            Assert.False(塩基列控え.Has控え(l_パス));
            Assert.Equal(new FastqReaderの塩基列(l_パス).ToList(), 塩基列控え.Get_塩基列(l_パス).ToList());
        }

        [Theory]
        [InlineData(false)]
        [InlineData(true)]
        public void 置き場経由で書き直したら控えを捨てる(bool p_Isオンメモリ)
        {
            中間データ置き場.A_Is有効 = p_Isオンメモリ;
            塩基列控え.A_置き場所 = Path.Combine(this._作業ディレクトリ, "bases");
            var l_パス = Path.Combine(this._作業ディレクトリ, "inter.fq");
            V_書出_置き場(l_パス, "ACGT");
            _ = 塩基列控え.Get_塩基列(l_パス).ToList();
            Assert.True(塩基列控え.Has控え(l_パス));

            V_書出_置き場(l_パス, "GGGG");

            Assert.False(塩基列控え.Has控え(l_パス));
            Assert.Equal(["GGGG"], 塩基列控え.Get_塩基列(l_パス).ToList());
        }

        [Fact]
        public void ディスク上で直接書き換わっても古い控えを返さない()
        {
            塩基列控え.A_置き場所 = Path.Combine(this._作業ディレクトリ, "bases");
            var l_パス = this.V_書出_FASTQ("reads.fq", ["ACGT"]);
            _ = 塩基列控え.Get_塩基列(l_パス).ToList();

            _ = this.V_書出_FASTQ("reads.fq", ["ACGTACGT", "TT"]);

            Assert.Equal(["ACGTACGT", "TT"], 塩基列控え.Get_塩基列(l_パス).ToList());
        }

        [Fact]
        public void 全消去で控えのファイルも消える()
        {
            var l_置き場所 = Path.Combine(this._作業ディレクトリ, "bases");
            塩基列控え.A_置き場所 = l_置き場所;
            var l_パス = this.V_書出_FASTQ("reads.fq", 塩基列群);
            _ = 塩基列控え.Get_塩基列(l_パス).ToList();
            Assert.NotEmpty(Directory.GetFiles(l_置き場所));

            塩基列控え.V_全消去();

            Assert.Empty(Directory.GetFiles(l_置き場所));
            Assert.False(塩基列控え.Has控え(l_パス));
        }

        private static void V_書出_置き場(string p_パス, string p_配列)
        {
            using var l_書き込み = new StreamWriter(中間データ置き場.Get_書込ストリーム(p_パス));
            l_書き込み.WriteLine("@r");
            l_書き込み.WriteLine(p_配列);
            l_書き込み.WriteLine("+");
            l_書き込み.WriteLine(new string('I', p_配列.Length));
        }

        // 控えを通さずに FASTQ から読んだ塩基列 (比較の基準)
        private sealed class FastqReaderの塩基列(string p_パス) : IEnumerable<string>
        {
            public IEnumerator<string> GetEnumerator()
            {
                using var l_読み込み = new FastqReader(p_パス);
                while (l_読み込み.Has続き())
                {
                    yield return l_読み込み.Get_次のレコード().A_配列;
                }
            }

            System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
            {
                return this.GetEnumerator();
            }
        }
    }
}
