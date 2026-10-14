using Tsumiki.Cores.Pipeline;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 作業ディレクトリの中間ファイルの削除 (gzip で置いた間引きリードを含む) の検証
    /// </summary>
    public class AssemblyWorkspaceTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki assembly workspace tests
        /// </summary>
        private const string C_項目_tsumiki_assembly_workspace_tests = "tsumiki_assembly_workspace_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 圧縮した間引きリード 1 側のファイル名
        /// </summary>
        private const string C_ファイル名_thinned_1_fq_gz = "thinned.1.fq.gz";

        /// <summary>
        /// 圧縮した間引きリード 2 側のチェックサムのファイル名
        /// </summary>
        private const string C_ファイル名_thinned_2_fq_gz_sha256 = "thinned.2.fq.gz.sha256";

        /// <summary>
        /// 中間ファイルではない FASTA のファイル名 (残るべき)
        /// </summary>
        private const string C_ファイル名_thinned_1_fasta = "thinned.1.fasta";

        #endregion

        #region 内部変数

        /// <summary>
        /// 作業用ディレクトリ
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の作業ディレクトリを用意する
        /// </summary>
        public AssemblyWorkspaceTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_assembly_workspace_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// gzip の間引きリードとそのチェックサムは削除され、中間ファイルでない FASTA は残る
        /// </summary>
        [Fact]
        public void V_削除_中間ファイル_gzipの間引きリードとsha256を消し_fastaは残す()
        {
            File.WriteAllText(Path.Combine(this._作業ディレクトリ, C_ファイル名_thinned_1_fq_gz), "gz");
            File.WriteAllText(Path.Combine(this._作業ディレクトリ, C_ファイル名_thinned_2_fq_gz_sha256), "sha");
            File.WriteAllText(Path.Combine(this._作業ディレクトリ, C_ファイル名_thinned_1_fasta), ">x");

            AssemblyWorkspace.V_削除_中間ファイル(this._作業ディレクトリ);

            var l_残り = Directory.EnumerateFiles(this._作業ディレクトリ).Select(Path.GetFileName).Order(StringComparer.Ordinal).ToList();
            Assert.Equal([C_ファイル名_thinned_1_fasta], l_残り);
        }

        /// <summary>
        /// 作業ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }

            GC.SuppressFinalize(this);
        }

        #endregion
    }
}
