using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Tests.Utility;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// 塩基列の控えから返す塩基列が、FASTQ をそのまま読んだときと同じであることの検証
    /// </summary>
    [Collection(中間データ置き場の集まり.C_名前)]
    public class 塩基列控えTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki bases cache tests
        /// </summary>
        private const string C_項目_tsumiki_bases_cache_tests = "tsumiki_bases_cache_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 項目 bases
        /// </summary>
        private const string C_項目_bases = "bases";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// ファイル名 inter fq
        /// </summary>
        private const string C_ファイル名_inter_fq = "inter.fq";

        /// <summary>
        /// 塩基配列 GGGG
        /// </summary>
        private const string C_塩基配列_GGGG = "GGGG";

        /// <summary>
        /// 塩基配列 ACGTACGT
        /// </summary>
        private const string C_塩基配列_ACGTACGT = "ACGTACGT";

        /// <summary>
        /// 塩基配列 TT
        /// </summary>
        private const string C_塩基配列_TT = "TT";

        /// <summary>
        /// 項目 A
        /// </summary>
        private const string C_項目_A = "A";

        /// <summary>
        /// FASTQ 品質区切り
        /// </summary>
        private const string C_FASTQ品質区切り = "+";

        /// <summary>
        /// FASTQ 見出し
        /// </summary>
        private const string C_FASTQ見出し = "@r";

        /// <summary>
        /// 塩基列群
        /// </summary>
        private static readonly string[] C_塩基列群 = [
            "ACGT",
            "A",
            "ACGTA",
            "NNACGTNN",
            "acgtACGT",
            string.Empty,
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            new string('G', 5000) + "C",
            "ACGTRYKM",
        ];

        #endregion

        #region 内部変数

        /// <summary>
        /// 作業ディレクトリ
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 塩基列控え Tests
        /// </summary>
        public 塩基列控えTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_bases_cache_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// Dispose
        /// </summary>
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

        /// <summary>
        /// 二回目は控えから同じ塩基列を返す
        /// </summary>
        /// <param name="p_Isオンメモリ"></param>
        [Theory]
        [InlineData(false)]
        [InlineData(true)]
        public void V_二回目は控えから同じ塩基列を返す(bool p_Isオンメモリ)
        {
            中間データ置き場.A_Is有効 = p_Isオンメモリ;
            塩基列控え.A_置き場所 = Path.Combine(this._作業ディレクトリ, C_項目_bases);
            var l_パス = this.V_書出_FASTQ(C_ファイル名_reads_fq, C_塩基列群);
            var l_期待 = new FastqReaderの塩基列(l_パス).ToList();
            var l_一回目 = 塩基列控え.Get_塩基列(l_パス).ToList();
            Assert.True(塩基列控え.Has控え(l_パス));
            var l_二回目 = 塩基列控え.Get_塩基列(l_パス).ToList();
            Assert.Equal(l_期待, l_一回目);
            Assert.Equal(l_期待, l_二回目);
        }

        /// <summary>
        /// 置き場所が無ければ控えを作らない
        /// </summary>
        [Fact]
        public void V_置き場所が無ければ控えを作らない()
        {
            var l_パス = this.V_書出_FASTQ(C_ファイル名_reads_fq, C_塩基列群);
            var l_結果 = 塩基列控え.Get_塩基列(l_パス).ToList();
            Assert.Equal(new FastqReaderの塩基列(l_パス).ToList(), l_結果);
            Assert.False(塩基列控え.Has控え(l_パス));
        }

        /// <summary>
        /// 途中で読むのをやめたら控えを登録しない
        /// </summary>
        [Fact]
        public void V_途中で読むのをやめたら控えを登録しない()
        {
            塩基列控え.A_置き場所 = Path.Combine(this._作業ディレクトリ, C_項目_bases);
            var l_パス = this.V_書出_FASTQ(C_ファイル名_reads_fq, C_塩基列群);
            _ = 塩基列控え.Get_塩基列(l_パス).Take(2).ToList();
            Assert.False(塩基列控え.Has控え(l_パス));
            Assert.Equal(new FastqReaderの塩基列(l_パス).ToList(), 塩基列控え.Get_塩基列(l_パス).ToList());
        }

        /// <summary>
        /// 置き場経由で書き直したら控えを捨てる
        /// </summary>
        /// <param name="p_Isオンメモリ"></param>
        [Theory]
        [InlineData(false)]
        [InlineData(true)]
        public void V_置き場経由で書き直したら控えを捨てる(bool p_Isオンメモリ)
        {
            中間データ置き場.A_Is有効 = p_Isオンメモリ;
            塩基列控え.A_置き場所 = Path.Combine(this._作業ディレクトリ, C_項目_bases);
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_inter_fq);
            V_書出_置き場(l_パス, Consts.塩基文字);
            _ = 塩基列控え.Get_塩基列(l_パス).ToList();
            Assert.True(塩基列控え.Has控え(l_パス));
            V_書出_置き場(l_パス, C_塩基配列_GGGG);
            Assert.False(塩基列控え.Has控え(l_パス));
            Assert.Equal([C_塩基配列_GGGG], 塩基列控え.Get_塩基列(l_パス).ToList());
        }

        /// <summary>
        /// ディスク上で直接書き換わっても古い控えを返さない
        /// </summary>
        [Fact]
        public void V_ディスク上で直接書き換わっても古い控えを返さない()
        {
            塩基列控え.A_置き場所 = Path.Combine(this._作業ディレクトリ, C_項目_bases);
            var l_パス = this.V_書出_FASTQ(C_ファイル名_reads_fq, [Consts.塩基文字]);
            _ = 塩基列控え.Get_塩基列(l_パス).ToList();
            _ = this.V_書出_FASTQ(C_ファイル名_reads_fq, [C_塩基配列_ACGTACGT, C_塩基配列_TT]);
            Assert.Equal([C_塩基配列_ACGTACGT, C_塩基配列_TT], 塩基列控え.Get_塩基列(l_パス).ToList());
        }

        /// <summary>
        /// 全消去で控えのファイルも消える
        /// </summary>
        [Fact]
        public void V_全消去で控えのファイルも消える()
        {
            var l_置き場所 = Path.Combine(this._作業ディレクトリ, C_項目_bases);
            塩基列控え.A_置き場所 = l_置き場所;
            var l_パス = this.V_書出_FASTQ(C_ファイル名_reads_fq, C_塩基列群);
            _ = 塩基列控え.Get_塩基列(l_パス).ToList();
            Assert.NotEmpty(Directory.GetFiles(l_置き場所));
            塩基列控え.V_全消去();
            Assert.Empty(Directory.GetFiles(l_置き場所));
            Assert.False(塩基列控え.Has控え(l_パス));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 書出_FASTQ
        /// </summary>
        /// <param name="p_名前"></param>
        /// <param name="p_配列群"></param>
        /// <returns></returns>
        private string V_書出_FASTQ(string p_名前, IEnumerable<string> p_配列群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_名前);
            using var l_書き込み = new StreamWriter(l_パス);
            var l_番号 = 0;
            foreach (var l_配列 in p_配列群)
            {
                l_書き込み.WriteLine($"@r{l_番号++}");
                l_書き込み.WriteLine(l_配列.Length == 0 ? C_項目_A : l_配列);
                l_書き込み.WriteLine(C_FASTQ品質区切り);
                l_書き込み.WriteLine(new string('I', l_配列.Length == 0 ? 1 : l_配列.Length));
            }

            return l_パス;
        }

        /// <summary>
        /// 書出_置き場
        /// </summary>
        /// <param name="p_パス"></param>
        /// <param name="p_配列"></param>
        private static void V_書出_置き場(string p_パス, string p_配列)
        {
            using var l_書き込み = new StreamWriter(中間データ置き場.Get_書込ストリーム(p_パス));
            l_書き込み.WriteLine(C_FASTQ見出し);
            l_書き込み.WriteLine(p_配列);
            l_書き込み.WriteLine(C_FASTQ品質区切り);
            l_書き込み.WriteLine(new string('I', p_配列.Length));
        }

        #endregion
    }

    /// <summary>
    /// 控えを通さずに FASTQ から塩基列を列挙する比較用データ
    /// </summary>
    /// <param name="p_パス"></param>
    internal sealed class FastqReaderの塩基列(string p_パス) : IEnumerable<string>
    {
        #region 公開メソッド

        /// <summary>
        /// GetEnumerator
        /// </summary>
        /// <returns></returns>
        public IEnumerator<string> GetEnumerator()
        {
            using var l_読み込み = new FastqReader(p_パス);
            while (l_読み込み.Has続き())
            {
                yield return l_読み込み.Get_次のレコード().A_配列;
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// GetEnumerator
        /// </summary>
        /// <returns></returns>
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return this.GetEnumerator();
        }

        #endregion
    }
}
