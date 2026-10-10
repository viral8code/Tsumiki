using Tsumiki.Cores.Evaluation;
using Tsumiki.IO;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終成果物の統計を Markdown の表に書き出す処理の検証
    /// </summary>
    public class AssemblyStatsTableTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki stats table tests
        /// </summary>
        private const string C_項目_tsumiki_stats_table_tests = "tsumiki_stats_table_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 assembly fasta
        /// </summary>
        private const string C_ファイル名_assembly_fasta = "assembly.fasta";

        /// <summary>
        /// 項目 assembly
        /// </summary>
        private const string C_項目_assembly = "assembly";

        /// <summary>
        /// 項目 scaffolds
        /// </summary>
        private const string C_項目_scaffolds = "scaffolds";

        /// <summary>
        /// ファイル名 none fasta
        /// </summary>
        private const string C_ファイル名_none_fasta = "none.fasta";

        /// <summary>
        /// 項目 assembly V 書き出し 統計表 ファイルごとに 3 行を書き存在しないファイルは飛ばす
        /// </summary>
        private const string C_項目_assembly_V_書き出し_統計表_ファイルごとに3行を書き存在しないファイルは飛ばす = "| assembly";

        /// <summary>
        /// 項目 assembly   all   2   1 310   1 010   300
        /// </summary>
        private const string C_項目_assembly___all___2___1_310___1_010___300 = "| assembly | all | 2 | 1,310 | 1,010 | 300 |";

        /// <summary>
        /// 項目 assembly      500bp   1   1 010
        /// </summary>
        private const string C_項目_assembly______500bp___1___1_010 = "| assembly | >= 500bp | 1 | 1,010 |";

        /// <summary>
        /// 項目 assembly   N split     500bp   1   600
        /// </summary>
        private const string C_項目_assembly___N_split_____500bp___1___600 = "| assembly | N-split, >= 500bp | 1 | 600 |";

        /// <summary>
        /// 項目 scaffolds V 書き出し 統計表 ファイルごとに 3 行を書き存在しないファイルは飛ばす
        /// </summary>
        private const string C_項目_scaffolds_V_書き出し_統計表_ファイルごとに3行を書き存在しないファイルは飛ばす = "| scaffolds";

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 一時ディレクトリを用意する
        /// </summary>
        public AssemblyStatsTableTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_stats_table_tests + Guid.NewGuid().ToString(C_GUID書式));
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
        /// ファイルごとに 3 通りの絞り込みの行を書き、存在しないファイルは飛ばす
        /// </summary>
        [Fact]
        public void V_書き出し_統計表_ファイルごとに3行を書き存在しないファイルは飛ばす()
        {
            var l_FASTAパス = Path.Combine(this._作業ディレクトリ, C_ファイル名_assembly_fasta);
            using (var l_書き込み = new FastaWriter(l_FASTAパス))
            {
                l_書き込み.V_書き込み(1, new string('A', 600) + new string('N', 10) + new string('C', 400));
                l_書き込み.V_書き込み(2, new string('G', 300));
            }

            var l_表 = AssemblyStatsReporter.Get_統計表([(C_項目_assembly, l_FASTAパス), (C_項目_scaffolds, Path.Combine(this._作業ディレクトリ, C_ファイル名_none_fasta))]);
            var l_行群 = l_表.Where(x => x.StartsWith(C_項目_assembly_V_書き出し_統計表_ファイルごとに3行を書き存在しないファイルは飛ばす, StringComparison.Ordinal)).ToList();
            Assert.Equal(3, l_行群.Count);
            Assert.StartsWith(C_項目_assembly___all___2___1_310___1_010___300, l_行群[0]);
            Assert.StartsWith(C_項目_assembly______500bp___1___1_010, l_行群[1]);
            Assert.StartsWith(C_項目_assembly___N_split_____500bp___1___600, l_行群[2]);
            Assert.DoesNotContain(l_表, x => x.StartsWith(C_項目_scaffolds_V_書き出し_統計表_ファイルごとに3行を書き存在しないファイルは飛ばす, StringComparison.Ordinal));
        }

        #endregion
    }
}
