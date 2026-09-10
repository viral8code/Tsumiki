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
    /// <remarks>
    /// 誤って繋いだ接合は両側それぞれが正しい配列なので局所の量では
    /// 見えず、繋ぎ目を跨ぐ r-mer の不在だけがそれを示す
    /// </remarks>
    public class ReadSupportCheckerTests : IDisposable
    {
        /// <summary>
        /// この検証で使う r-mer 長
        /// </summary>
        private const int R = 31;

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public ReadSupportCheckerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_support_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
            ConfigurationManager.A_実行時引数 = new Parameters { A_スレッド数 = 1 };
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
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string Get_乱数配列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }

        /// <summary>
        /// 元の配列を 100 bp のリードで隙間なく覆った FASTQ を作る
        /// </summary>
        /// <param name="p_名前">ファイル名</param>
        /// <param name="p_元">元になる配列</param>
        /// <returns>書き出したパス</returns>
        private string Get_リード(string p_名前, params string[] p_元)
        {
            var l_パス = Path.Combine(this._tempDir, p_名前);
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
            var l_パス = Path.Combine(this._tempDir, "asm.fasta");
            using var l_書き込み = new FastaWriter(l_パス);
            foreach (var (l_ID, l_配列) in p_全件)
            {
                l_書き込み.V_書き込み(l_ID, l_配列);
            }
            return l_パス;
        }

        /// <summary>
        /// リードどおりの配列なら支持のない位置は出ないことを確かめる
        /// </summary>
        [Fact]
        public void Get_検査結果_リードどおりの配列なら支持のない位置は出ない()
        {
            var l_真値 = Get_乱数配列(2000, p_種: 20260913);
            var l_FASTA = this.Get_FASTA(("SEQ1", l_真値));
            var l_リード = this.Get_リード("reads.fq", l_真値);

            var l_結果 = ReadSupportChecker.Get_検査結果(l_FASTA, l_リード, null, R);

            Assert.NotNull(l_結果);
            Assert.Equal(0, l_結果!.Value.A_支持のない位置数);
            Assert.Empty(l_結果.Value.A_区間);
        }

        /// <summary>
        /// 無関係な 2 本を繋いだ接合を 1 つの区間として指すことを確かめる
        /// </summary>
        [Fact]
        public void Get_検査結果_無関係な2本を繋いだ接合を1つの区間として指す()
        {
            var l_左 = Get_乱数配列(1000, p_種: 20260914);
            var l_右 = Get_乱数配列(1000, p_種: 20260915);

            // リードは左右それぞれからしか出ない
            // 繋いだ接合を読んだリードは無い
            var l_リード = this.Get_リード("reads.fq", l_左, l_右);
            var l_FASTA = this.Get_FASTA(("SEQ1", l_左 + l_右));

            var l_結果 = ReadSupportChecker.Get_検査結果(l_FASTA, l_リード, null, R);

            Assert.NotNull(l_結果);
            var l_区間 = Assert.Single(l_結果!.Value.A_区間);
            Assert.Equal("SEQ1", l_区間.A_配列ID);

            // 接合を跨ぐ r-mer は R-1 個
            // 覆う塩基は接合の両側 R-1 塩基ぶん
            Assert.Equal(R - 1, l_結果.Value.A_支持のない位置数);
            Assert.Equal(l_左.Length - R + 2, l_区間.A_開始);
            Assert.Equal(l_左.Length + R - 1, l_区間.A_終了);
        }

        /// <summary>
        /// ギャップの N は支持を問わないことを確かめる
        /// </summary>
        [Fact]
        public void Get_検査結果_ギャップのNは支持を問わない()
        {
            var l_真値 = Get_乱数配列(2000, p_種: 20260916);
            var l_リード = this.Get_リード("reads.fq", l_真値);
            var l_FASTA = this.Get_FASTA(
                ("SEQ1", l_真値[..1000] + new string('N', 50) + l_真値[1000..]));

            var l_結果 = ReadSupportChecker.Get_検査結果(l_FASTA, l_リード, null, R);

            Assert.NotNull(l_結果);
            Assert.Equal(0, l_結果!.Value.A_支持のない位置数);
            Assert.Empty(l_結果.Value.A_区間);
        }

        /// <summary>
        /// r が長すぎる場合は調べないことを確かめる
        /// </summary>
        [Fact]
        public void Get_検査結果_rが長すぎる場合は調べない()
        {
            var l_真値 = Get_乱数配列(500, p_種: 20260917);
            var l_FASTA = this.Get_FASTA(("SEQ1", l_真値));
            var l_リード = this.Get_リード("reads.fq", l_真値);

            Assert.Null(ReadSupportChecker.Get_検査結果(l_FASTA, l_リード, null, p_r長: 65));
        }
    }
}
