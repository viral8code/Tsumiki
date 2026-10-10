using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// リードの貼り付けで、通り抜けた unitig の並びが入口と出口の対応を保ったまま数えられることの検証
    /// </summary>
    public class ReadPathIndexMappingTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki readpath tests
        /// </summary>
        private const string C_項目_tsumiki_readpath_tests = "tsumiki_readpath_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 項目 A
        /// </summary>
        private const string C_項目_A = "A";

        /// <summary>
        /// 項目 C
        /// </summary>
        private const string C_項目_C = "C";

        /// <summary>
        /// 項目 G
        /// </summary>
        private const string C_項目_G = "G";

        /// <summary>
        /// 項目 T
        /// </summary>
        private const string C_項目_T = "T";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// ファイル名 r1 fq
        /// </summary>
        private const string C_ファイル名_r1_fq = "r1.fq";

        /// <summary>
        /// ファイル名 r2 fq
        /// </summary>
        private const string C_ファイル名_r2_fq = "r2.fq";

        /// <summary>
        /// ファイル名 unitigs fasta
        /// </summary>
        private const string C_ファイル名_unitigs_fasta = "unitigs.fasta";

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k長 = 8;

        /// <summary>
        /// 読み取り窓の長さ
        /// </summary>
        private const int C_リード長 = 30;

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        /// <summary>
        /// A・C・R・B・D の順の unitig 配列
        /// </summary>
        private readonly string[] _unitig群;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// A-R-B と C-R-D が R を共有する構成を用意する
        /// </summary>
        public ReadPathIndexMappingTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_readpath_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
                A_スレッド数 = 1
            };
            var l_反復 = RepeatGraphFixture.Get_乱数配列(16, 2001);
            var l_A = RepeatGraphFixture.Get_乱数配列(11, 2002) + C_項目_A + l_反復[..(C_k長 - 1)];
            var l_C = RepeatGraphFixture.Get_乱数配列(11, 2003) + C_項目_C + l_反復[..(C_k長 - 1)];
            var l_B = l_反復[^(C_k長 - 1)..] + C_項目_G + RepeatGraphFixture.Get_乱数配列(11, 2004);
            var l_D = l_反復[^(C_k長 - 1)..] + C_項目_T + RepeatGraphFixture.Get_乱数配列(11, 2005);
            this._unitig群 = [l_A, l_C, l_反復, l_B, l_D];
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
        /// A-R-B と C-R-D を読んだリードからは、その 2 通りの並びだけが数えられ、A-R-D は 0 になる
        /// </summary>
        [Fact]
        public void V_反復を通り抜けた入口と出口の組だけが数えられる()
        {
            var l_contig構築 = new ContigMaker(this.V_書き込み_unitigFASTA());
            var l_リードパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_reads_fq);
            V_書き込み_FASTQ(l_リードパス, [..Get_窓群(this.Get_連結(0, 2, 3)), ..Get_窓群(this.Get_連結(1, 2, 4))]);
            l_contig構築.V_マッピング_リード(l_リードパス);
            var l_索引 = l_contig構築.Get_経路索引();
            int[] l_ARB = [Get_頂点(1), Get_頂点(3), Get_頂点(4)];
            int[] l_CRD = [Get_頂点(2), Get_頂点(3), Get_頂点(5)];
            int[] l_ARD = [Get_頂点(1), Get_頂点(3), Get_頂点(5)];
            Assert.Equal((ulong)Get_窓群(this.Get_連結(0, 2, 3)).Count, l_索引.Get_出現数(l_ARB));
            Assert.Equal((ulong)Get_窓群(this.Get_連結(1, 2, 4)).Count, l_索引.Get_出現数(l_CRD));
            Assert.Equal(0UL, l_索引.Get_出現数(l_ARD));
            Assert.Equal(l_索引.Get_出現数(l_ARB), l_索引.Get_出現数([Get_頂点(4) ^ 1, Get_頂点(3) ^ 1, Get_頂点(1) ^ 1]));
        }

        /// <summary>
        /// read1 と read2 (逆相補で読まれる) が同じ並びを読んだペアは、1 組を 1 回と数える
        /// </summary>
        [Fact]
        public void V_同じペアの両リードが同じ並びを読んでも二重に数えない()
        {
            var l_contig構築 = new ContigMaker(this.V_書き込み_unitigFASTA());
            var l_連結 = this.Get_連結(0, 2, 3);
            var l_窓群 = Get_窓群(l_連結);
            var l_順リードパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_r1_fq);
            var l_逆リードパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_r2_fq);
            V_書き込み_FASTQ(l_順リードパス, l_窓群);
            V_書き込み_FASTQ(l_逆リードパス, [..l_窓群.Select(Util.V_逆相補)]);
            l_contig構築.V_マッピング_ペアリード(l_順リードパス, l_逆リードパス);
            var l_索引 = l_contig構築.Get_経路索引();
            Assert.Equal((ulong)l_窓群.Count, l_索引.Get_出現数([Get_頂点(1), Get_頂点(3), Get_頂点(4)]));
        }

        /// <summary>
        /// 反復を複製して解いたあとも、並びは移った側の複製を指す
        /// </summary>
        [Fact]
        public void V_反復を複製したあとも並びが複製側の頂点を指す()
        {
            var l_contig構築 = new ContigMaker(this.V_書き込み_unitigFASTA());
            var l_リードパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_reads_fq);
            V_書き込み_FASTQ(l_リードパス, [..Get_窓群(this.Get_連結(0, 2, 3)), ..Get_窓群(this.Get_連結(1, 2, 4))]);
            l_contig構築.V_マッピング_リード(l_リードパス);
            var l_索引 = l_contig構築.Get_経路索引();
            var l_件数 = l_索引.Get_出現数([Get_頂点(1), Get_頂点(3), Get_頂点(4)]);
            var l_グラフ = l_contig構築.Get_グラフ();
            List<string> l_unitig配列 = [string.Empty, string.Empty, ..this._unitig群.SelectMany(x => new[] { x, Util.V_逆相補(x) })];
            var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig配列, [], new Dictionary<(int, int), ulong>(), p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL, p_経路索引: l_索引);
            Assert.Equal(1, l_解決数);
            var l_Aの次 = Assert.Single(l_グラフ.A_出辺[Get_頂点(1)]);
            var l_Cの次 = Assert.Single(l_グラフ.A_出辺[Get_頂点(2)]);
            Assert.NotEqual(l_Aの次, l_Cの次);
            Assert.Equal(l_件数, l_索引.Get_出現数([Get_頂点(1), l_Aの次, Get_頂点(4)]));
            Assert.Equal(l_件数, l_索引.Get_出現数([Get_頂点(2), l_Cの次, Get_頂点(5)]));
            Assert.Equal(0UL, l_索引.Get_出現数([Get_頂点(1), l_Cの次, Get_頂点(4)]));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// unitig ID の順鎖の頂点番号
        /// </summary>
        /// <param name="p_ID"></param>
        /// <returns></returns>
        private static int Get_頂点(int p_ID)
        {
            return ContigMaker.Get_頂点番号(p_ID);
        }

        /// <summary>
        /// 添字で選んだ unitig を k-1 の重なりで繋いだ配列
        /// </summary>
        /// <param name="p_添字群"></param>
        /// <returns></returns>
        private string Get_連結(params int[] p_添字群)
        {
            return this._unitig群[p_添字群[0]] + string.Concat(p_添字群.Skip(1).Select(x => this._unitig群[x][(C_k長 - 1)..]));
        }

        /// <summary>
        /// 3 unitig をすべて跨ぐ読み取り窓 (先頭の unitig の末尾 k-mer と末尾の unitig の先頭 k-mer を含むもの)
        /// </summary>
        /// <param name="p_連結"></param>
        /// <returns></returns>
        private static List<string> Get_窓群(string p_連結)
        {
            const int l_先頭の末尾kmer開始 = 19 - C_k長;
            const int l_末尾の先頭kmer終了 = 19 + 16 - (2 * (C_k長 - 1)) + C_k長;
            List<string> l_窓群 = [];
            for (var s = 0; s + C_リード長 <= p_連結.Length; s++)
            {
                if (s <= l_先頭の末尾kmer開始 && s + C_リード長 >= l_末尾の先頭kmer終了)
                {
                    l_窓群.Add(p_連結.Substring(s, C_リード長));
                }
            }

            return l_窓群;
        }

        /// <summary>
        /// 配列群を FASTQ に書き出す
        /// </summary>
        /// <param name="p_パス"></param>
        /// <param name="p_配列群"></param>
        private static void V_書き込み_FASTQ(string p_パス, IReadOnlyList<string> p_配列群)
        {
            using var l_書き込み = new StreamWriter(p_パス);
            for (var i = 0; i < p_配列群.Count; i++)
            {
                l_書き込み.WriteLine($"@r{i}\n{p_配列群[i]}\n+\n{new string('I', p_配列群[i].Length)}");
            }
        }

        /// <summary>
        /// unitig FASTA を書き出す
        /// </summary>
        /// <returns>書き出したパス</returns>
        private string V_書き込み_unitigFASTA()
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_unitigs_fasta);
            using var l_書き込み = new FastaWriter(l_パス);
            for (var i = 0; i < this._unitig群.Length; i++)
            {
                l_書き込み.V_書き込み(i + 1, this._unitig群[i]);
            }

            return l_パス;
        }

        #endregion
    }
}
