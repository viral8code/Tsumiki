using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 前段 k で確定した経路を、この k の分岐選択へ投影する仕組み (P2) の検証
    /// </summary>
    public class PathCarryOverTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki pathcarryover tests
        /// </summary>
        private const string C_項目_tsumiki_pathcarryover_tests = "tsumiki_pathcarryover_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 contigs fasta
        /// </summary>
        private const string C_ファイル名_contigs_fasta = "contigs.fasta";

        /// <summary>
        /// 塩基配列 GGGTTAG
        /// </summary>
        private const string C_塩基配列_GGGTTAG = "GGGTTAG";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// 塩基配列 CATCGAA
        /// </summary>
        private const string C_塩基配列_CATCGAA = "CATCGAA";

        /// <summary>
        /// ファイル名 contigs2 fasta
        /// </summary>
        private const string C_ファイル名_contigs2_fasta = "contigs2.fasta";

        /// <summary>
        /// ファイル名 unitigs fasta
        /// </summary>
        private const string C_ファイル名_unitigs_fasta = "unitigs.fasta";

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k長 = 8;

        /// <summary>
        /// 分岐元の入口配列 (末尾 7 塩基が分岐先双方の先頭と重なる)
        /// </summary>
        private const string C_入口unitig = "TTTTTTTAAACCCG";

        /// <summary>
        /// 分岐先 1 (入口の末尾 7 塩基を共有する)
        /// </summary>
        private const string C_分岐先1unitig = "AAACCCGGGGTTAG";

        /// <summary>
        /// 分岐先 2 (入口の末尾 7 塩基を共有する、分岐先 1 とは中身が異なる)
        /// </summary>
        private const string C_分岐先2unitig = "AAACCCGCATCGAA";

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
        public PathCarryOverTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_pathcarryover_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
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
        }

        /// <summary>
        /// この k 自身の read/pair 支持が無く分岐を決められない場合でも、前段 k から引き継いだ経路が 正しい分岐対応 (入口 → 分岐先 1) だけを選ばせることを確かめる
        /// </summary>
        [Fact]
        public void V_このkの支持だけでは決められない分岐を経路引き継ぎが解決する()
        {
            var l_unitigパス = this.V_書き込み_unitigFASTA();
            var l_contigパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_contigs_fasta);
            var l_contig構築 = new ContigMaker(l_unitigパス);
            var l_引き継ぎ経路 = new[]
            {
                C_入口unitig + C_塩基配列_GGGTTAG
            };
            l_contig構築.V_結合_Contig(l_contigパス, p_優勢閾値: 0.8M, p_最小証拠数: 1UL, p_引き継ぎ経路群: l_引き継ぎ経路);
            var l_contig群 = FastaReader.Get_全エントリ(l_contigパス).Select(x => x.A_配列).ToList();
            var l_重なり長 = C_k長 - 1;
            var l_結合後の長さ = C_入口unitig.Length + C_分岐先1unitig.Length - l_重なり長;
            Assert.Equal(2, l_contig群.Count);
            Assert.Contains(l_contig群, x => x.Length == l_結合後の長さ);
            Assert.Contains(l_contig群, x => x.Length == C_分岐先2unitig.Length);
        }

        /// <summary>
        /// この k 自身の read 支持が経路引き継ぎと矛盾する場合、実測の支持が必ず勝つことを確かめる
        /// </summary>
        [Fact]
        public void V_実測支持と矛盾する経路引き継ぎは採用されない()
        {
            var l_unitigパス = this.V_書き込み_unitigFASTA();
            var l_contigパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_contigs_fasta);
            var l_リードパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_reads_fq);
            this.V_書き込み_FASTQ(l_リードパス, C_入口unitig + C_塩基配列_CATCGAA, p_本数: 5);
            var l_contig構築 = new ContigMaker(l_unitigパス);
            l_contig構築.V_マッピング_リード(l_リードパス);
            var l_引き継ぎ経路 = new[]
            {
                C_入口unitig + C_塩基配列_GGGTTAG
            };
            var l_contigパス2 = Path.Combine(this._作業ディレクトリ, C_ファイル名_contigs2_fasta);
            l_contig構築.V_結合_Contig(l_contigパス2, p_優勢閾値: 0.8M, p_最小証拠数: 1UL, p_引き継ぎ経路群: l_引き継ぎ経路);
            var l_contig群 = FastaReader.Get_全エントリ(l_contigパス2).Select(x => x.A_配列).ToList();
            var l_重なり長 = C_k長 - 1;
            var l_結合後の長さ = C_入口unitig.Length + C_分岐先2unitig.Length - l_重なり長;
            Assert.Equal(2, l_contig群.Count);
            Assert.Contains(l_contig群, x => x.Length == l_結合後の長さ);
            Assert.Contains(l_contig群, x => x.Length == C_分岐先1unitig.Length);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 入口・分岐先 1・分岐先 2 の 3 本からなる unitig FASTA を書き出す
        /// </summary>
        /// <returns>書き出したパス</returns>
        private string V_書き込み_unitigFASTA()
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_unitigs_fasta);
            using var l_書き込み = new FastaWriter(l_パス);
            l_書き込み.V_書き込み(1, C_入口unitig);
            l_書き込み.V_書き込み(2, C_分岐先1unitig);
            l_書き込み.V_書き込み(3, C_分岐先2unitig);
            return l_パス;
        }

        /// <summary>
        /// 同じ配列を指定本数だけ含む FASTQ を書き出す
        /// </summary>
        /// <param name="p_パス">書き出し先</param>
        /// <param name="p_配列">書き出す配列</param>
        /// <param name="p_本数">繰り返す本数</param>
        private void V_書き込み_FASTQ(string p_パス, string p_配列, int p_本数)
        {
            using var l_書き込み = new StreamWriter(p_パス);
            for (var i = 0; i < p_本数; i++)
            {
                l_書き込み.WriteLine($"@r{i}\n{p_配列}\n+\n{new string('I', p_配列.Length)}");
            }
        }

        #endregion
    }
}
