using Tsumiki.Cores.Pipeline;
using Tsumiki.IO;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 統合結果から contig を書き出すときの N の扱い
    /// </summary>
    public class MergedContigSplitTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki merged contig
        /// </summary>
        private const string C_項目_tsumiki_merged_contig = "tsumiki_merged_contig_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 merged scaffolds fasta
        /// </summary>
        private const string C_ファイル名_merged_scaffolds_fasta = "merged_scaffolds.fasta";

        /// <summary>
        /// ファイル名 merged contigs fasta
        /// </summary>
        private const string C_ファイル名_merged_contigs_fasta = "merged_contigs.fasta";

        /// <summary>
        /// 書式 SCAFFOLD1 ACGTACGTNNNGGT AFFOLD2 TTTTAAAA
        /// </summary>
        private const string C_書式_SCAFFOLD1_ACGTACGTNNNGGT_AFFOLD2_TTTTAAAA = ">SCAFFOLD1\nACGTACGTNNNGGTTCCAA\n>SCAFFOLD2\nTTTTAAAA\n";

        /// <summary>
        /// 塩基配列 ACGTACGT
        /// </summary>
        private const string C_塩基配列_ACGTACGT = "ACGTACGT";

        /// <summary>
        /// 塩基配列 GGTTCCAA
        /// </summary>
        private const string C_塩基配列_GGTTCCAA = "GGTTCCAA";

        /// <summary>
        /// 塩基配列 TTTTAAAA
        /// </summary>
        private const string C_塩基配列_TTTTAAAA = "TTTTAAAA";

        /// <summary>
        /// 書式 SCAFFOLD1 NNACGTACGTN
        /// </summary>
        private const string C_書式_SCAFFOLD1_NNACGTACGTN = ">SCAFFOLD1\nNNACGTACGTN\n";

        #endregion

        #region 内部変数

        /// <summary>
        /// テスト専用の作業パス
        /// </summary>
        private readonly string _作業パス = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_merged_contig + Guid.NewGuid().ToString(C_GUID書式));

        #endregion

        #region コンストラクタ

        /// <summary>
        /// テスト用ディレクトリを用意する
        /// </summary>
        public MergedContigSplitTests()
        {
            _ = Directory.CreateDirectory(this._作業パス);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// テスト用ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            Directory.Delete(this._作業パス, recursive: true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// N の連続で分断し、N を含まない片だけを書き出す
        /// </summary>
        [Fact]
        public void V_書き出し_N分割_Nで分断して書き出す()
        {
            var l_入力 = Path.Combine(this._作業パス, C_ファイル名_merged_scaffolds_fasta);
            var l_出力 = Path.Combine(this._作業パス, C_ファイル名_merged_contigs_fasta);
            File.WriteAllText(l_入力, C_書式_SCAFFOLD1_ACGTACGTNNNGGT_AFFOLD2_TTTTAAAA);
            MultiKAssembler.V_書き出し_N分割(l_入力, l_出力);
            var l_配列群 = FastaReader.Get_全エントリ(l_出力).Select(x => x.A_配列).ToList();
            Assert.Equal([C_塩基配列_ACGTACGT, C_塩基配列_GGTTCCAA, C_塩基配列_TTTTAAAA], l_配列群);
        }

        /// <summary>
        /// 端の N は空の片を作らない
        /// </summary>
        [Fact]
        public void V_書き出し_N分割_端のNは空の片を作らない()
        {
            var l_入力 = Path.Combine(this._作業パス, C_ファイル名_merged_scaffolds_fasta);
            var l_出力 = Path.Combine(this._作業パス, C_ファイル名_merged_contigs_fasta);
            File.WriteAllText(l_入力, C_書式_SCAFFOLD1_NNACGTACGTN);
            MultiKAssembler.V_書き出し_N分割(l_入力, l_出力);
            var l_配列 = Assert.Single(FastaReader.Get_全エントリ(l_出力).Select(x => x.A_配列));
            Assert.Equal(C_塩基配列_ACGTACGT, l_配列);
        }

        #endregion
    }
}
