using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアエンドから推定される「インサートサイズ」が、リードに挟まれた内側の未読区間ではなく、真のフラグメント長 (左リードの 5'端から右リードの 3'端まで) の単位になっていることを検証する
    /// </summary>
    public class InsertSizeEstimationTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki insertsize tests
        /// </summary>
        private const string C_項目_tsumiki_insertsize_tests = "tsumiki_insertsize_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 unitigs fasta
        /// </summary>
        private const string C_ファイル名_unitigs_fasta = "unitigs.fasta";

        /// <summary>
        /// ファイル名 r1 fq
        /// </summary>
        private const string C_ファイル名_r1_fq = "r1.fq";

        /// <summary>
        /// ファイル名 r2 fq
        /// </summary>
        private const string C_ファイル名_r2_fq = "r2.fq";

        /// <summary>
        /// 項目 pair 1
        /// </summary>
        private const string C_項目_pair_1 = "pair/1";

        /// <summary>
        /// 項目 pair 2
        /// </summary>
        private const string C_項目_pair_2 = "pair/2";

        /// <summary>
        /// FASTQ 品質区切り
        /// </summary>
        private const string C_FASTQ品質区切り = "+";

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _一時ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public InsertSizeEstimationTests()
        {
            this._一時ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_insertsize_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._一時ディレクトリ))
            {
                Directory.Delete(this._一時ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// V_SameUnitigSamples_MeasureFullFragmentLength_NotTheInnerDistance
        /// </summary>
        [Fact]
        public void V_SameUnitigSamples_MeasureFullFragmentLength_NotTheInnerDistance()
        {
            const int l_k長 = 21;
            const int l_unitig長 = 600;
            const int l_リード長 = 50;
            const int l_断片開始位置 = 100;
            const int l_真の断片長 = 350;
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_k長,
                A_スレッド数 = 1
            };
            var l_unitig配列 = V_生成_乱数配列(l_unitig長, p_シード: 12_345);
            var l_unitigパス = Path.Combine(this._一時ディレクトリ, C_ファイル名_unitigs_fasta);
            File.WriteAllText(l_unitigパス, $">1\n{l_unitig配列}\n");
            var l_先行リード = l_unitig配列.Substring(l_断片開始位置, l_リード長);
            var l_後続リード = Util.V_逆相補(l_unitig配列.Substring(l_断片開始位置 + l_真の断片長 - l_リード長, l_リード長));
            var l_先行パス = Path.Combine(this._一時ディレクトリ, C_ファイル名_r1_fq);
            var l_後続パス = Path.Combine(this._一時ディレクトリ, C_ファイル名_r2_fq);
            var l_ペア群 = Enumerable.Range(0, 5).ToList();
            V_書き込み_Fastq(l_先行パス, l_ペア群.Select(i => ($"pair{i}/1", l_先行リード)));
            V_書き込み_Fastq(l_後続パス, l_ペア群.Select(i => ($"pair{i}/2", l_後続リード)));
            var l_contig構築 = new ContigMaker(l_unitigパス);
            l_contig構築.V_マッピング_ペアリード(l_先行パス, l_後続パス);
            Assert.NotEmpty(l_contig構築.A_同一unitig標本);
            Assert.All(l_contig構築.A_同一unitig標本, l_標本 => Assert.Equal(l_真の断片長, l_標本));
        }

        /// <summary>
        /// フラグメント長を変えたときに推定値が同じだけ動くこと (定数ぶんのずれではなく、単位そのものが一致していること) を確認する
        /// </summary>
        /// <param name="p_真の断片長"></param>
        [Theory]
        [InlineData(200)]
        [InlineData(350)]
        [InlineData(500)]
        public void V_SameUnitigSamples_TrackTheActualFragmentLength(int p_真の断片長)
        {
            const int l_k長 = 21;
            const int l_unitig長 = 900;
            const int l_リード長 = 50;
            const int l_断片開始位置 = 120;
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = l_k長,
                A_スレッド数 = 1
            };
            var l_unitig配列 = V_生成_乱数配列(l_unitig長, p_シード: 777);
            var l_unitigパス = Path.Combine(this._一時ディレクトリ, $"unitigs_{p_真の断片長}.fasta");
            File.WriteAllText(l_unitigパス, $">1\n{l_unitig配列}\n");
            var l_先行リード = l_unitig配列.Substring(l_断片開始位置, l_リード長);
            var l_後続リード = Util.V_逆相補(l_unitig配列.Substring(l_断片開始位置 + p_真の断片長 - l_リード長, l_リード長));
            var l_先行パス = Path.Combine(this._一時ディレクトリ, $"r1_{p_真の断片長}.fq");
            var l_後続パス = Path.Combine(this._一時ディレクトリ, $"r2_{p_真の断片長}.fq");
            V_書き込み_Fastq(l_先行パス, [(C_項目_pair_1, l_先行リード)]);
            V_書き込み_Fastq(l_後続パス, [(C_項目_pair_2, l_後続リード)]);
            var l_contig構築 = new ContigMaker(l_unitigパス);
            l_contig構築.V_マッピング_ペアリード(l_先行パス, l_後続パス);
            Assert.NotEmpty(l_contig構築.A_同一unitig標本);
            Assert.All(l_contig構築.A_同一unitig標本, l_標本 => Assert.Equal(p_真の断片長, l_標本));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 決定的な擬似乱数で非反復的な塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">生成する配列長</param>
        /// <param name="p_シード">乱数シード</param>
        /// <returns></returns>
        private static string V_生成_乱数配列(int p_長さ, int p_シード)
        {
            var l_乱数 = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数.Next(4)]));
        }

        /// <summary>
        /// リードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_パス">書き出し先</param>
        /// <param name="p_リード一覧">書き出すリード</param>
        private static void V_書き込み_Fastq(string p_パス, IEnumerable<(string A_ID, string A_配列)> p_リード一覧)
        {
            using var l_ライター = new StreamWriter(p_パス);
            foreach (var (l_ID, l_配列)in p_リード一覧)
            {
                l_ライター.WriteLine($"@{l_ID}");
                l_ライター.WriteLine(l_配列);
                l_ライター.WriteLine(C_FASTQ品質区切り);
                l_ライター.WriteLine(new string('I', l_配列.Length));
            }
        }

        #endregion
    }
}
