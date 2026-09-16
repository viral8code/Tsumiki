using Tsumiki.Cores.Pipeline;
using Tsumiki.IO;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 統合結果から contig を書き出すときの N の扱い
    /// </summary>
    /// <remarks>
    /// 統合は scaffold 同士を橋渡しするので統合結果には N が残る。contig は N を含まない連続配列として出す
    /// </remarks>
    public class MergedContigSplitTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// テスト専用の作業パス
        /// </summary>
        private readonly string _作業パス = Path.Combine(Path.GetTempPath(), "tsumiki_merged_contig_" + Guid.NewGuid().ToString("N"));

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
            var l_入力 = Path.Combine(this._作業パス, "merged_scaffolds.fasta");
            var l_出力 = Path.Combine(this._作業パス, "merged_contigs.fasta");
            File.WriteAllText(l_入力, ">SCAFFOLD1\nACGTACGTNNNGGTTCCAA\n>SCAFFOLD2\nTTTTAAAA\n");

            MultiKAssembler.V_書き出し_N分割(l_入力, l_出力);

            var l_配列群 = FastaReader.Get_全エントリ(l_出力).Select(x => x.A_配列).ToList();

            Assert.Equal(["ACGTACGT", "GGTTCCAA", "TTTTAAAA"], l_配列群);
        }

        /// <summary>
        /// 端の N は空の片を作らない
        /// </summary>
        [Fact]
        public void V_書き出し_N分割_端のNは空の片を作らない()
        {
            var l_入力 = Path.Combine(this._作業パス, "merged_scaffolds.fasta");
            var l_出力 = Path.Combine(this._作業パス, "merged_contigs.fasta");
            File.WriteAllText(l_入力, ">SCAFFOLD1\nNNACGTACGTN\n");

            MultiKAssembler.V_書き出し_N分割(l_入力, l_出力);

            var l_配列 = Assert.Single(FastaReader.Get_全エントリ(l_出力).Select(x => x.A_配列));

            Assert.Equal("ACGTACGT", l_配列);
        }

        #endregion
    }
}
