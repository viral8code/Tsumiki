using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 選んだアセンブリに無い配列を他の k から補う処理の検証
    /// </summary>
    public class AssemblyComplementerTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki complementer tests
        /// </summary>
        private const string C_項目_tsumiki_complementer_tests = "tsumiki_complementer_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 backbone fasta
        /// </summary>
        private const string C_ファイル名_backbone_fasta = "backbone.fasta";

        /// <summary>
        /// ファイル名 other fasta
        /// </summary>
        private const string C_ファイル名_other_fasta = "other.fasta";

        /// <summary>
        /// ファイル名 rescued fasta
        /// </summary>
        private const string C_ファイル名_rescued_fasta = "rescued.fasta";

        /// <summary>
        /// 項目 119
        /// </summary>
        private const string C_項目_119 = "119_";

        /// <summary>
        /// ファイル名 other1 fasta
        /// </summary>
        private const string C_ファイル名_other1_fasta = "other1.fasta";

        /// <summary>
        /// ファイル名 other2 fasta
        /// </summary>
        private const string C_ファイル名_other2_fasta = "other2.fasta";

        /// <summary>
        /// 塩基配列 NNNNN
        /// </summary>
        private const string C_塩基配列_NNNNN = "NNNNN";

        /// <summary>
        /// アンカー k 長
        /// </summary>
        private const int C_アンカーk長 = 19;

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 一時ディレクトリを作る
        /// </summary>
        public AssemblyComplementerTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_complementer_tests + Guid.NewGuid().ToString(C_GUID書式));
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
        /// 骨格に無い配列を他の k が持っていれば、その区間を骨格の後ろに足すことを確かめる
        /// </summary>
        [Fact]
        public void Get_補完_骨格に無い配列を他のkから足す()
        {
            var l_左 = Get_乱数配列(3_000, 701);
            var l_欠け = Get_乱数配列(1_500, 702);
            var l_右 = Get_乱数配列(3_000, 703);
            var l_骨格 = this.Get_アセンブリ(C_ファイル名_backbone_fasta, 41, l_左, l_右);
            var l_他 = this.Get_アセンブリ(C_ファイル名_other_fasta, 119, l_左 + l_欠け + l_右);
            var l_出力 = Path.Combine(this._作業ディレクトリ, C_ファイル名_rescued_fasta);
            var (l_本数, l_延長) = AssemblyComplementer.Get_補完(l_骨格, [l_骨格, l_他], static _ => true, C_アンカーk長, l_出力);
            var l_結果 = FastaReader.Get_全エントリ(l_出力);
            Assert.Equal(1, l_本数);
            Assert.Equal(l_欠け.Length + (2 * (C_アンカーk長 - 1)), l_延長);
            Assert.Equal(3, l_結果.Count);
            Assert.Equal(l_左, l_結果[0].A_配列);
            Assert.Equal(l_右, l_結果[1].A_配列);
            Assert.StartsWith(Consts.補完配列の接頭辞 + C_項目_119, l_結果[2].A_ID);
            Assert.Contains(l_欠け, l_結果[2].A_配列);
        }

        /// <summary>
        /// 最小長に届かない区間は足さず、出力も書かないことを確かめる
        /// </summary>
        [Fact]
        public void Get_補完_短い区間は足さない()
        {
            var l_左 = Get_乱数配列(3_000, 711);
            var l_欠け = Get_乱数配列(AssemblyComplementer.C_補う最小長 - (2 * C_アンカーk長), 712);
            var l_右 = Get_乱数配列(3_000, 713);
            var l_骨格 = this.Get_アセンブリ(C_ファイル名_backbone_fasta, 41, l_左, l_右);
            var l_他 = this.Get_アセンブリ(C_ファイル名_other_fasta, 119, l_左 + l_欠け + l_右);
            var l_出力 = Path.Combine(this._作業ディレクトリ, C_ファイル名_rescued_fasta);
            var (l_本数, _) = AssemblyComplementer.Get_補完(l_骨格, [l_骨格, l_他], static _ => true, C_アンカーk長, l_出力);
            Assert.Equal(0, l_本数);
            Assert.False(File.Exists(l_出力));
        }

        /// <summary>
        /// リードで信頼できない k-mer ばかりの区間は足さないことを確かめる
        /// </summary>
        [Fact]
        public void Get_補完_信頼できない配列は足さない()
        {
            var l_左 = Get_乱数配列(3_000, 721);
            var l_欠け = Get_乱数配列(1_500, 722);
            var l_右 = Get_乱数配列(3_000, 723);
            var l_骨格 = this.Get_アセンブリ(C_ファイル名_backbone_fasta, 41, l_左, l_右);
            var l_他 = this.Get_アセンブリ(C_ファイル名_other_fasta, 119, l_左 + l_欠け + l_右);
            var l_出力 = Path.Combine(this._作業ディレクトリ, C_ファイル名_rescued_fasta);
            var (l_本数, _) = AssemblyComplementer.Get_補完(l_骨格, [l_骨格, l_他], static _ => false, C_アンカーk長, l_出力);
            Assert.Equal(0, l_本数);
        }

        /// <summary>
        /// 複数の k が同じ配列を持っていても 1 回だけ、大きい k のものを足すことを確かめる
        /// </summary>
        [Fact]
        public void Get_補完_同じ配列は大きいkから1回だけ足す()
        {
            var l_左 = Get_乱数配列(3_000, 731);
            var l_欠け = Get_乱数配列(1_500, 732);
            var l_右 = Get_乱数配列(3_000, 733);
            var l_骨格 = this.Get_アセンブリ(C_ファイル名_backbone_fasta, 41, l_左, l_右);
            var l_他1 = this.Get_アセンブリ(C_ファイル名_other1_fasta, 83, l_左 + l_欠け + l_右);
            var l_他2 = this.Get_アセンブリ(C_ファイル名_other2_fasta, 119, Util.V_逆相補(l_左 + l_欠け + l_右));
            var l_出力 = Path.Combine(this._作業ディレクトリ, C_ファイル名_rescued_fasta);
            var (l_本数, _) = AssemblyComplementer.Get_補完(l_骨格, [l_他1, l_骨格, l_他2], static _ => true, C_アンカーk長, l_出力);
            var l_結果 = FastaReader.Get_全エントリ(l_出力);
            Assert.Equal(1, l_本数);
            Assert.StartsWith(Consts.補完配列の接頭辞 + C_項目_119, l_結果[^1].A_ID);
        }

        /// <summary>
        /// 骨格と同じ配列しか持たない候補からは何も足さないことを確かめる
        /// </summary>
        [Fact]
        public void Get_補完_骨格に既にある配列は足さない()
        {
            var l_配列 = Get_乱数配列(6_000, 741);
            var l_骨格 = this.Get_アセンブリ(C_ファイル名_backbone_fasta, 41, l_配列);
            var l_他 = this.Get_アセンブリ(C_ファイル名_other_fasta, 119, Util.V_逆相補(l_配列));
            var l_出力 = Path.Combine(this._作業ディレクトリ, C_ファイル名_rescued_fasta);
            var (l_本数, _) = AssemblyComplementer.Get_補完(l_骨格, [l_骨格, l_他], static _ => true, C_アンカーk長, l_出力);
            Assert.Equal(0, l_本数);
        }

        /// <summary>
        /// 骨格にある短い反復で割れた欠けは、つないで 1 本として足すことを確かめる
        /// </summary>
        [Fact]
        public void Get_補完_短い既知区間を挟んだ欠けはつないで足す()
        {
            var l_左 = Get_乱数配列(3_000, 761);
            var l_反復 = Get_乱数配列(300, 762);
            var l_欠け1 = Get_乱数配列(400, 763);
            var l_欠け2 = Get_乱数配列(400, 764);
            var l_右 = Get_乱数配列(3_000, 765);
            var l_骨格 = this.Get_アセンブリ(C_ファイル名_backbone_fasta, 41, l_左 + l_反復, l_右);
            var l_他 = this.Get_アセンブリ(C_ファイル名_other_fasta, 119, l_左 + l_欠け1 + l_反復 + l_欠け2 + l_右);
            var l_出力 = Path.Combine(this._作業ディレクトリ, C_ファイル名_rescued_fasta);
            var (l_本数, _) = AssemblyComplementer.Get_補完(l_骨格, [l_骨格, l_他], static _ => true, C_アンカーk長, l_出力);
            var l_結果 = FastaReader.Get_全エントリ(l_出力);
            Assert.Equal(1, l_本数);
            Assert.Contains(l_欠け1 + l_反復 + l_欠け2, l_結果[^1].A_配列);
        }

        /// <summary>
        /// 既知の区間が大半を占めるなら、まばらな欠けをつないで足さないことを確かめる
        /// </summary>
        [Fact]
        public void Get_つないだ区間_既知が大半なら返さない()
        {
            var l_配列 = Get_乱数配列(2_000, 771);
            var l_区間 = AssemblyComplementer.Get_つないだ区間(l_配列, [(0, 37), (400, 37), (800, 37), (1_200, 37)]).ToList();
            Assert.Empty(l_区間);
        }

        /// <summary>
        /// N を挟んだ欠けはつながないことを確かめる
        /// </summary>
        [Fact]
        public void Get_つないだ区間_Nを挟めばつながない()
        {
            var l_配列 = Get_乱数配列(600, 781) + C_塩基配列_NNNNN + Get_乱数配列(600, 782);
            var l_区間 = AssemblyComplementer.Get_つないだ区間(l_配列, [(0, 600), (605, 600)]).ToList();
            Assert.Equal([(0, 600), (605, 600)], l_区間);
        }

        /// <summary>
        /// 欠けた区間は N で途切れることを確かめる
        /// </summary>
        [Fact]
        public void Get_欠けた区間_Nで途切れる()
        {
            var l_配列 = Get_乱数配列(100, 751) + C_塩基配列_NNNNN + Get_乱数配列(200, 752);
            var l_区間 = AssemblyComplementer.Get_欠けた区間(l_配列, new HashSet<UInt128>(), C_アンカーk長).ToList();
            Assert.Equal([(0, 100), (105, 200)], l_区間);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 再現できる乱数配列を作る
        /// </summary>
        /// <param name="p_長さ">配列長</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>配列</returns>
        private static string Get_乱数配列(int p_長さ, int p_シード)
        {
            var l_乱数 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数.Next(4)]));
        }

        /// <summary>
        /// 配列を FASTA に書き出し、アセンブリの実行結果として返す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_配列群">書き出す配列</param>
        /// <returns>アセンブリの実行結果</returns>
        private アセンブリ実行結果 Get_アセンブリ(string p_ファイル名, int p_k長, params string[] p_配列群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_ファイル名);
            using (var l_書き込み = new FastaWriter(l_パス))
            {
                var l_連番 = 1;
                foreach (var l_配列 in p_配列群)
                {
                    l_書き込み.V_書き込み($"NODE{l_連番++}", l_配列);
                }
            }

            return new アセンブリ実行結果(p_k長, l_パス, l_パス, null, 2UL, 20.0D);
        }

        #endregion
    }
}
