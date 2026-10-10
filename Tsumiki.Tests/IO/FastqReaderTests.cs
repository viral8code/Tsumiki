using Tsumiki.Commons;
using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// 壊れた FASTQ で範囲外アクセスや無限ループにせず、どこが不正かを言って止まることを固定する
    /// </summary>
    public class FastqReaderTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki fastq reader tests
        /// </summary>
        private const string C_項目_tsumiki_fastq_reader_tests = "tsumiki_fastq_reader_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 書式 r1 ACGTACGT   IIIIIIII
        /// </summary>
        private const string C_書式_r1_ACGTACGT___IIIIIIII = "@r1\nACGTACGT\n+\nIIIIIIII\n";

        /// <summary>
        /// 項目 r1
        /// </summary>
        private const string C_項目_r1 = "@r1";

        /// <summary>
        /// 塩基配列 ACGTACGT
        /// </summary>
        private const string C_塩基配列_ACGTACGT = "ACGTACGT";

        /// <summary>
        /// 書式 r1 ACGTACGT   IIII
        /// </summary>
        private const string C_書式_r1_ACGTACGT___IIII = "@r1\nACGTACGT\n+\nIIII\n";

        /// <summary>
        /// 書式 r1 ACGTACGT
        /// </summary>
        private const string C_書式_r1_ACGTACGT = "@r1\nACGTACGT\n+\n";

        /// <summary>
        /// 書式 r1 ACGTACGT   IIIII
        /// </summary>
        private const string C_書式_r1_ACGTACGT___IIIII = "@r1\nACGTACGT\n+\nIIIII\n";

        /// <summary>
        /// 書式 a 1 ACGT   IIII  b 1 CCC  c 1 GGGG   IIII
        /// </summary>
        private const string C_書式_a_1_ACGT___IIII__b_1_CCC__c_1_GGGG___IIII = "@a/1\nACGT\n+\nIIII\n@b/1\nCCCC\n+\nIIII\n@c/1\nGGGG\n+\nIIII\n";

        /// <summary>
        /// 書式 a 2 TTTT   IIII  b 2 AAAA   IIII
        /// </summary>
        private const string C_書式_a_2_TTTT___IIII__b_2_AAAA___IIII = "@a/2\nTTTT\n+\nIIII\n@b/2\nAAAA\n+\nIIII\n";

        /// <summary>
        /// 塩基配列 TTTT
        /// </summary>
        private const string C_塩基配列_TTTT = "TTTT";

        /// <summary>
        /// 塩基配列 CCCC
        /// </summary>
        private const string C_塩基配列_CCCC = "CCCC";

        /// <summary>
        /// 塩基配列 AAAA
        /// </summary>
        private const string C_塩基配列_AAAA = "AAAA";

        /// <summary>
        /// 塩基配列 GGGG
        /// </summary>
        private const string C_塩基配列_GGGG = "GGGG";

        /// <summary>
        /// ファイル名 fastq
        /// </summary>
        private const string C_ファイル名_fastq = ".fastq";

        #endregion

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
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_fastq_reader_tests + Guid.NewGuid().ToString(C_GUID書式));
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
            var l_パス = this.Get_書き出し先(C_書式_r1_ACGTACGT___IIIIIIII);
            using var l_読み込み = new FastqReader(l_パス);
            Assert.True(l_読み込み.Has続き());
            var l_リード = l_読み込み.Get_次のリード_軽量();
            Assert.Equal(C_項目_r1, l_リード.A_ID);
            Assert.Equal(C_塩基配列_ACGTACGT, l_リード.A_生リード);
        }

        /// <summary>
        /// クオリティ長が配列長と異なる場合は例外を投げることを検証する
        /// </summary>
        [Fact]
        public void Get_次のリード_クオリティ長が配列長と異なる場合は例外を投げる()
        {
            var l_パス = this.Get_書き出し先(C_書式_r1_ACGTACGT___IIII);
            using var l_読み込み = new FastqReader(l_パス);
            var l_例外 = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード_軽量());
            Assert.Contains(C_項目_r1, l_例外.Message);
        }

        /// <summary>
        /// レコードの途中でファイルが途切れている場合は例外を投げることを検証する
        /// </summary>
        [Fact]
        public void Get_次のリード_レコードの途中でファイルが途切れている場合は例外を投げる()
        {
            var l_パス = this.Get_書き出し先(C_書式_r1_ACGTACGT);
            using var l_読み込み = new FastqReader(l_パス);
            _ = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード_軽量());
        }

        /// <summary>
        /// 曖昧塩基を扱う経路でも同様に長さの不一致を検査することを検証する
        /// </summary>
        [Fact]
        public void Get_次のリード_曖昧塩基を扱う経路でも同様に検査する()
        {
            var l_パス = this.Get_書き出し先(C_書式_r1_ACGTACGT___IIIII);
            using var l_読み込み = new FastqReader(l_パス);
            _ = Assert.Throws<InvalidDataException>(() => l_読み込み.Get_次のリード());
        }

        /// <summary>
        /// ペアの 2 ファイルを同じ順に組にして流し、片方が先に尽きたら残りを相方なしで流すことを検証する
        /// </summary>
        [Fact]
        public void Get_ペア塩基列_同じ順に組にし余りは相方なしで流す()
        {
            var l_パス1 = this.Get_書き出し先(C_書式_a_1_ACGT___IIII__b_1_CCC__c_1_GGGG___IIII);
            var l_パス2 = this.Get_書き出し先(C_書式_a_2_TTTT___IIII__b_2_AAAA___IIII);
            var l_組群 = FastqReader.Get_ペア塩基列(l_パス1, l_パス2).ToList();
            Assert.Equal([(Consts.塩基文字, C_塩基配列_TTTT), (C_塩基配列_CCCC, C_塩基配列_AAAA), (C_塩基配列_GGGG, string.Empty)], l_組群);
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
            var l_パス = Path.Combine(this._作業ディレクトリ, Guid.NewGuid().ToString(C_GUID書式) + C_ファイル名_fastq);
            File.WriteAllText(l_パス, p_内容);
            return l_パス;
        }

        #endregion
    }
}
