using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 出した配列の各位置がリードに裏付けられているかの検査
    /// </summary>
    public class ReadSupportCheckerTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki support tests
        /// </summary>
        private const string C_項目_tsumiki_support_tests = "tsumiki_support_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 項目 SEQ1
        /// </summary>
        private const string C_項目_SEQ1 = "SEQ1";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// ファイル名 asm fasta
        /// </summary>
        private const string C_ファイル名_asm_fasta = "asm.fasta";

        /// <summary>
        /// この検証で使う r-mer 長
        /// </summary>
        private const int C_r長 = 31;

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
        public ReadSupportCheckerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_support_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_スレッド数 = 1
            };
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

            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// リードどおりの配列なら支持のない位置は出ないことを確かめる
        /// </summary>
        [Fact]
        public void Get_検査結果_リードどおりの配列なら支持のない位置は出ない()
        {
            var l_真値 = Get_乱数配列(2_000, p_種: 20_260_913);
            var l_FASTA = this.Get_FASTA((C_項目_SEQ1, l_真値));
            var l_リード = this.Get_リード(C_ファイル名_reads_fq, l_真値);
            var l_結果 = ReadSupportChecker.Get_検査結果(l_FASTA, [(l_リード, string.Empty)], C_r長);
            Assert.NotNull(l_結果);
            Assert.Equal(0L, l_結果!.Value.A_支持のない位置数);
            Assert.Empty(l_結果.Value.A_区間);
        }

        /// <summary>
        /// 無関係な 2 本を繋いだ接合を 1 つの区間として指すことを確かめる
        /// </summary>
        [Fact]
        public void Get_検査結果_無関係な2本を繋いだ接合を1つの区間として指す()
        {
            var l_左 = Get_乱数配列(1_000, p_種: 20_260_914);
            var l_右 = Get_乱数配列(1_000, p_種: 20_260_915);
            var l_リード = this.Get_リード(C_ファイル名_reads_fq, l_左, l_右);
            var l_FASTA = this.Get_FASTA((C_項目_SEQ1, l_左 + l_右));
            var l_結果 = ReadSupportChecker.Get_検査結果(l_FASTA, [(l_リード, string.Empty)], C_r長);
            Assert.NotNull(l_結果);
            var l_区間 = Assert.Single(l_結果!.Value.A_区間);
            Assert.Equal(C_項目_SEQ1, l_区間.A_配列ID);
            Assert.Equal(C_r長 - 1, l_結果.Value.A_支持のない位置数);
            Assert.Equal(l_左.Length - C_r長 + 2, l_区間.A_開始);
            Assert.Equal(l_左.Length + C_r長 - 1, l_区間.A_終了);
        }

        /// <summary>
        /// ギャップの N は支持を問わないことを確かめる
        /// </summary>
        [Fact]
        public void Get_検査結果_ギャップのNは支持を問わない()
        {
            var l_真値 = Get_乱数配列(2_000, p_種: 20_260_916);
            var l_リード = this.Get_リード(C_ファイル名_reads_fq, l_真値);
            var l_FASTA = this.Get_FASTA((C_項目_SEQ1, l_真値[..1_000] + new string('N', 50) + l_真値[1_000..]));
            var l_結果 = ReadSupportChecker.Get_検査結果(l_FASTA, [(l_リード, string.Empty)], C_r長);
            Assert.NotNull(l_結果);
            Assert.Equal(0L, l_結果!.Value.A_支持のない位置数);
            Assert.Empty(l_結果.Value.A_区間);
        }

        /// <summary>
        /// r が長すぎる場合は調べないことを確かめる
        /// </summary>
        [Fact]
        public void Get_検査結果_rが長すぎる場合は調べない()
        {
            var l_真値 = Get_乱数配列(500, p_種: 20_260_917);
            var l_FASTA = this.Get_FASTA((C_項目_SEQ1, l_真値));
            var l_リード = this.Get_リード(C_ファイル名_reads_fq, l_真値);
            Assert.Null(ReadSupportChecker.Get_検査結果(l_FASTA, [(l_リード, string.Empty)], p_r長: 65));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数.Next(4)]));
        }

        /// <summary>
        /// 元の配列を 100 bp のリードで隙間なく覆った FASTQ を作る
        /// </summary>
        /// <param name="p_名前">ファイル名</param>
        /// <param name="p_元">元になる配列</param>
        /// <returns>書き出したパス</returns>
        private string Get_リード(string p_名前, params string[] p_元)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_名前);
            using var l_書き込み = new FastqWriter(l_パス);
            var l_ID = 0;
            foreach (var l_配列 in p_元)
            {
                for (var i = 0; i + 100 <= l_配列.Length; i += 10)
                {
                    l_書き込み.V_書き込み($"@r{l_ID++}", l_配列.Substring(i, 100), new string('I', 100));
                }
            }

            return l_パス;
        }

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_全件">書き出す名前と配列</param>
        /// <returns>書き出したパス</returns>
        private string Get_FASTA(params (string A_ID, string A_配列)[] p_全件)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_asm_fasta);
            using var l_書き込み = new FastaWriter(l_パス);
            foreach (var (l_ID, l_配列)in p_全件)
            {
                l_書き込み.V_書き込み(l_ID, l_配列);
            }

            return l_パス;
        }

        #endregion
    }
}
