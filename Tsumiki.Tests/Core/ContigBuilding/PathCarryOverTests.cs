using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 前段 k で確定した経路を、この k の分岐選択へ投影する仕組み (P2) の検証
    /// </summary>
    /// <remarks>
    /// 決定的な受入基準:<br/>
    /// (1) 全 k-mer が既にこの k に存在していて集合の持ち越しだけでは何も変わらない場合でも、経路引き継ぎだけが正しい分岐対応を選べること<br/>
    /// (2) この k 自身の read/pair 支持が経路引き継ぎと矛盾する場合は、必ず実測の支持が勝つこと (引き継ぎが誤接続を持ち込まない)
    /// </remarks>
    public class PathCarryOverTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 8;

        /// <summary>
        /// 分岐元の入口配列 (末尾 7 塩基が分岐先双方の先頭と重なる)
        /// </summary>
        private const string 入口ユニティグ = "TTTTTTTAAACCCG";

        /// <summary>
        /// 分岐先 1 (入口の末尾 7 塩基を共有する)
        /// </summary>
        private const string 分岐先1ユニティグ = "AAACCCGGGGTTAG";

        /// <summary>
        /// 分岐先 2 (入口の末尾 7 塩基を共有する、分岐先 1 とは中身が異なる)
        /// </summary>
        private const string 分岐先2ユニティグ = "AAACCCGCATCGAA";

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
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_pathcarryover_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
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
        /// この k 自身の read/pair 支持が無く分岐を決められない場合でも、前段 k から引き継いだ経路が
        /// 正しい分岐対応 (入口 → 分岐先1) だけを選ばせることを確かめる
        /// </summary>
        [Fact]
        public void V_このkの支持だけでは決められない分岐を経路引き継ぎが解決する()
        {
            var l_unitigパス = this.V_書き込み_unitigFASTA();
            var l_contigパス = Path.Combine(this._作業ディレクトリ, "contigs.fasta");
            var l_contig構築 = new ContigMaker(l_unitigパス);

            // この k 自身の read/pair マッピングは一切行わない (real 支持ゼロ)
            // 前段 k で確定した経路 (入口+分岐先1 の全体配列) だけを引き継ぐ
            var l_引き継ぎ経路 = new[] { 入口ユニティグ + "GGGTTAG" };

            l_contig構築.V_結合_Contig(l_contigパス, p_優勢閾値: 0.8M, p_最小証拠数: 1UL, p_引き継ぎ経路群: l_引き継ぎ経路);

            var l_contig群 = FastaReader.Get_全エントリ(l_contigパス).Select(x => x.A_配列).ToList();

            // 入口 + 分岐先1 が 1 本に結合され、分岐先2 は単独で残るはず (2 本になる)
            // 出力の鎖の向き (順鎖/逆鎖どちらを正準として選ぶか) は walk の実装詳細なので、
            // 内容の部分一致ではなく長さで結合の有無を確かめる
            var l_重なり長 = k長 - 1;
            var l_結合後の長さ = 入口ユニティグ.Length + 分岐先1ユニティグ.Length - l_重なり長;
            Assert.Equal(2, l_contig群.Count);
            Assert.Contains(l_contig群, x => x.Length == l_結合後の長さ);
            Assert.Contains(l_contig群, x => x.Length == 分岐先2ユニティグ.Length);
        }

        /// <summary>
        /// この k 自身の read 支持が経路引き継ぎと矛盾する場合、実測の支持が必ず勝つことを確かめる
        /// </summary>
        /// <remarks>
        /// 低k由来の誤接続を高kが矛盾により棄却する、という P2 のもう一つの受入基準
        /// </remarks>
        [Fact]
        public void V_実測支持と矛盾する経路引き継ぎは採用されない()
        {
            var l_unitigパス = this.V_書き込み_unitigFASTA();
            var l_contigパス = Path.Combine(this._作業ディレクトリ, "contigs.fasta");
            var l_リードパス = Path.Combine(this._作業ディレクトリ, "reads.fq");

            // この k で実際に観測された (と仮定する) read は、入口 -> 分岐先2 を繰り返し裏付ける
            this.V_書き込み_FASTQ(l_リードパス, 入口ユニティグ + "CATCGAA", p_本数: 5);

            var l_contig構築 = new ContigMaker(l_unitigパス);
            l_contig構築.V_マッピング_リード(l_リードパス);

            // 経路引き継ぎは (誤って) 分岐先1 を示している、という矛盾した状況
            var l_引き継ぎ経路 = new[] { 入口ユニティグ + "GGGTTAG" };

            var l_contigパス2 = Path.Combine(this._作業ディレクトリ, "contigs2.fasta");
            l_contig構築.V_結合_Contig(l_contigパス2, p_優勢閾値: 0.8M, p_最小証拠数: 1UL, p_引き継ぎ経路群: l_引き継ぎ経路);

            var l_contig群 = FastaReader.Get_全エントリ(l_contigパス2).Select(x => x.A_配列).ToList();

            // 矛盾する引き継ぎを無視し、実測支持どおり 入口+分岐先2 が結合され、分岐先1 が単独で残る
            var l_重なり長 = k長 - 1;
            var l_結合後の長さ = 入口ユニティグ.Length + 分岐先2ユニティグ.Length - l_重なり長;
            Assert.Equal(2, l_contig群.Count);
            Assert.Contains(l_contig群, x => x.Length == l_結合後の長さ);
            Assert.Contains(l_contig群, x => x.Length == 分岐先1ユニティグ.Length);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 入口・分岐先1・分岐先2 の 3 本からなる unitig FASTA を書き出す
        /// </summary>
        /// <returns>書き出したパス</returns>
        private string V_書き込み_unitigFASTA()
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, "unitigs.fasta");
            using var l_書き込み = new FastaWriter(l_パス);
            l_書き込み.V_書き込み(1, 入口ユニティグ);
            l_書き込み.V_書き込み(2, 分岐先1ユニティグ);
            l_書き込み.V_書き込み(3, 分岐先2ユニティグ);
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
