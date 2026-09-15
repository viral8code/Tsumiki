using Tsumiki.Cores.Evaluation;
using Tsumiki.IO;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終成果物の統計を Markdown の表に書き出す処理の検証
    /// </summary>
    public class AssemblyStatsTableTests : IDisposable
    {
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
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_stats_table_tests_" + Guid.NewGuid().ToString("N"));
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
            var l_FASTAパス = Path.Combine(this._作業ディレクトリ, "assembly.fasta");
            using (var l_書き込み = new FastaWriter(l_FASTAパス))
            {
                l_書き込み.V_書き込み(1, new string('A', 600) + new string('N', 10) + new string('C', 400));
                l_書き込み.V_書き込み(2, new string('G', 300));
            }
            var l_表 = AssemblyStatsReporter.Get_統計表([("assembly", l_FASTAパス), ("scaffolds", Path.Combine(this._作業ディレクトリ, "none.fasta"))]);

            var l_行群 = l_表.Where(x => x.StartsWith("| assembly", StringComparison.Ordinal)).ToList();
            Assert.Equal(3, l_行群.Count);
            Assert.StartsWith("| assembly | all | 2 | 1,310 | 1,010 | 300 |", l_行群[0]);
            Assert.StartsWith("| assembly | >= 500bp | 1 | 1,010 |", l_行群[1]);
            Assert.StartsWith("| assembly | N-split, >= 500bp | 1 | 600 |", l_行群[2]);
            Assert.DoesNotContain(l_表, x => x.StartsWith("| scaffolds", StringComparison.Ordinal));
        }

        #endregion
    }
}
