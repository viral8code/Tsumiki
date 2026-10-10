using Tsumiki.Commons;
using Tsumiki.IO;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終成果物から、リード長より短い配列を落とす処理の検証
    /// </summary>
    public class ShortSequenceFilterTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki shortfilter tests
        /// </summary>
        private const string C_項目_tsumiki_shortfilter_tests = "tsumiki_shortfilter_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 項目 SCAFFOLD1
        /// </summary>
        private const string C_項目_SCAFFOLD1 = "SCAFFOLD1";

        /// <summary>
        /// 項目 SCAFFOLD2
        /// </summary>
        private const string C_項目_SCAFFOLD2 = "SCAFFOLD2";

        /// <summary>
        /// 項目 SCAFFOLD3
        /// </summary>
        private const string C_項目_SCAFFOLD3 = "SCAFFOLD3";

        /// <summary>
        /// 項目 SCAFFOLD4
        /// </summary>
        private const string C_項目_SCAFFOLD4 = "SCAFFOLD4";

        /// <summary>
        /// ファイル名 unitigs fasta
        /// </summary>
        private const string C_ファイル名_unitigs_fasta = "unitigs.fasta";

        /// <summary>
        /// ファイル名 contigs fasta
        /// </summary>
        private const string C_ファイル名_contigs_fasta = "contigs.fasta";

        /// <summary>
        /// ファイル名 corrected 1 fq
        /// </summary>
        private const string C_ファイル名_corrected_1_fq = "corrected.1.fq";

        /// <summary>
        /// ファイル名 preprocessed 1 fq
        /// </summary>
        private const string C_ファイル名_preprocessed_1_fq = "preprocessed.1.fq";

        /// <summary>
        /// ファイル名 polished fasta
        /// </summary>
        private const string C_ファイル名_polished_fasta = "polished.fasta";

        /// <summary>
        /// ファイル名 merged scaffolds fasta
        /// </summary>
        private const string C_ファイル名_merged_scaffolds_fasta = "merged_scaffolds.fasta";

        /// <summary>
        /// 項目 k93
        /// </summary>
        private const string C_項目_k93 = "k93";

        /// <summary>
        /// 項目 x
        /// </summary>
        private const string C_項目_x = "x";

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
        public ShortSequenceFilterTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_shortfilter_tests + Guid.NewGuid().ToString(C_GUID書式));
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

            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// リードより短い配列だけを落とすことを確かめる
        /// </summary>
        [Fact]
        public void V_除外_短い配列_リードより短い配列だけを落とす()
        {
            var l_パス = this.Get_書き出し((C_項目_SCAFFOLD1, new string('A', 300)), (C_項目_SCAFFOLD2, new string('C', 150)), (C_項目_SCAFFOLD3, new string('G', 149)), (C_項目_SCAFFOLD4, new string('T', 1)));
            Tsumiki.Cores.Pipeline.FinalAssemblyPipeline.V_除外_短い配列(l_パス, 150);
            var l_残り = FastaReader.Get_全エントリ(l_パス);
            Assert.Equal([C_項目_SCAFFOLD1, C_項目_SCAFFOLD2], l_残り.Select(x => x.A_ID));
            Assert.Equal(300, l_残り[0].A_配列.Length);
            Assert.Equal(150, l_残り[1].A_配列.Length);
        }

        /// <summary>
        /// 短い配列が無ければファイルはそのままであることを確かめる
        /// </summary>
        [Fact]
        public void V_除外_短い配列_短い配列が無ければファイルはそのまま()
        {
            var l_パス = this.Get_書き出し((C_項目_SCAFFOLD1, new string('A', 300)), (C_項目_SCAFFOLD2, new string('C', 200)));
            var l_元 = File.ReadAllBytes(l_パス);
            Tsumiki.Cores.Pipeline.FinalAssemblyPipeline.V_除外_短い配列(l_パス, 150);
            Assert.Equal(l_元, File.ReadAllBytes(l_パス));
        }

        /// <summary>
        /// リード長が不明なら何もしないことを確かめる
        /// </summary>
        [Fact]
        public void V_除外_短い配列_リード長が不明なら何もしない()
        {
            var l_パス = this.Get_書き出し((C_項目_SCAFFOLD1, new string('A', 10)));
            var l_元 = File.ReadAllBytes(l_パス);
            Tsumiki.Cores.Pipeline.FinalAssemblyPipeline.V_除外_短い配列(l_パス, null);
            Assert.Equal(l_元, File.ReadAllBytes(l_パス));
        }

        /// <summary>
        /// 最終成果物とログだけを残し、中間ファイルを削除することを確かめる
        /// </summary>
        [Fact]
        public void V_削除_中間ファイル_最終成果物とログだけ残す()
        {
            foreach (var l_名前 in new[]
            {
                C_ファイル名_unitigs_fasta,
                C_ファイル名_contigs_fasta,
                Consts.Scaffoldファイル名,
                Consts.GFAファイル名,
                Consts.レポートファイル名,
                Consts.曖昧箇所ファイル名,
                Consts.ログファイル名,
                C_ファイル名_corrected_1_fq,
                C_ファイル名_preprocessed_1_fq,
                C_ファイル名_polished_fasta,
                C_ファイル名_merged_scaffolds_fasta,
            }

            )
            {
                File.WriteAllText(Path.Combine(this._作業ディレクトリ, l_名前), l_名前);
            }

            _ = Directory.CreateDirectory(Path.Combine(this._作業ディレクトリ, C_項目_k93));
            File.WriteAllText(Path.Combine(this._作業ディレクトリ, C_項目_k93, C_ファイル名_unitigs_fasta), C_項目_x);
            Tsumiki.Cores.Pipeline.AssemblyWorkspace.V_削除_中間ファイル(this._作業ディレクトリ);
            var l_残り = Directory.EnumerateFiles(this._作業ディレクトリ).Select(Path.GetFileName).Order(StringComparer.Ordinal).ToList();
            Assert.Equal([
                Consts.ログファイル名,
                Consts.曖昧箇所ファイル名,
                Consts.GFAファイル名,
                Consts.レポートファイル名,
                C_ファイル名_contigs_fasta,
                Consts.Scaffoldファイル名,
                C_ファイル名_unitigs_fasta,
            ], l_残り);
            Assert.Empty(Directory.EnumerateDirectories(this._作業ディレクトリ));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_全件">書き出す名前と配列</param>
        /// <returns>書き出したパス</returns>
        private string Get_書き出し(params (string A_ID, string A_配列)[] p_全件)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, Consts.Scaffoldファイル名);
            using (var l_書き込み = new FastaWriter(l_パス))
            {
                foreach (var (l_ID, l_配列)in p_全件)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }

            return l_パス;
        }

        #endregion
    }
}
