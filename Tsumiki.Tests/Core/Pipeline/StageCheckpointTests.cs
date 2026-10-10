using Tsumiki.Cores.Pipeline;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Evaluation;
using Tsumiki.Commons;
using System.Text.Json;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 再開時の内容照合と実行境界を検証する
    /// </summary>
    public class StageCheckpointTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki checkpoint
        /// </summary>
        private const string C_項目_tsumiki_checkpoint = "tsumiki_checkpoint_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 raw fq
        /// </summary>
        private const string C_ファイル名_raw_fq = "raw.fq";

        /// <summary>
        /// ファイル名 source fasta
        /// </summary>
        private const string C_ファイル名_source_fasta = "source.fasta";

        /// <summary>
        /// 項目 output
        /// </summary>
        private const string C_項目_output = "output";

        /// <summary>
        /// 項目 complete
        /// </summary>
        private const string C_項目_complete = "complete";

        /// <summary>
        /// 項目 self check
        /// </summary>
        private const string C_項目_self_check = "self_check";

        /// <summary>
        /// 項目 missing kmers
        /// </summary>
        private const string C_項目_missing_kmers = "missing_kmers";

        /// <summary>
        /// 項目 checks
        /// </summary>
        private const string C_項目_checks = "checks";

        /// <summary>
        /// 項目 name
        /// </summary>
        private const string C_項目_name = "name";

        /// <summary>
        /// 項目 read support
        /// </summary>
        private const string C_項目_read_support = "read_support";

        /// <summary>
        /// 項目 result
        /// </summary>
        private const string C_項目_result = "result";

        /// <summary>
        /// 項目 fail
        /// </summary>
        private const string C_項目_fail = "fail";

        /// <summary>
        /// 項目 Weighted
        /// </summary>
        private const string C_項目_Weighted = "Weighted";

        /// <summary>
        /// 項目 assembly settings
        /// </summary>
        private const string C_項目_assembly_settings = "assembly_settings";

        /// <summary>
        /// 項目 copy number baseline requested
        /// </summary>
        private const string C_項目_copy_number_baseline_requested = "copy_number_baseline_requested";

        /// <summary>
        /// 項目 copy number baseline actual
        /// </summary>
        private const string C_項目_copy_number_baseline_actual = "copy_number_baseline_actual";

        /// <summary>
        /// 項目 trim low coverage ends
        /// </summary>
        private const string C_項目_trim_low_coverage_ends = "trim_low_coverage_ends";

        /// <summary>
        /// ファイル名 assembly provenance json
        /// </summary>
        private const string C_ファイル名_assembly_provenance_json = "assembly.provenance.json";

        /// <summary>
        /// ファイル名 assembly fasta
        /// </summary>
        private const string C_ファイル名_assembly_fasta = "assembly.fasta";

        /// <summary>
        /// 項目 assembly sha256
        /// </summary>
        private const string C_項目_assembly_sha256 = "assembly_sha256";

        /// <summary>
        /// ファイル名 corrected fq
        /// </summary>
        private const string C_ファイル名_corrected_fq = "corrected.fq";

        /// <summary>
        /// 塩基配列 AAAA
        /// </summary>
        private const string C_塩基配列_AAAA = "AAAA";

        /// <summary>
        /// 塩基配列 CCCC
        /// </summary>
        private const string C_塩基配列_CCCC = "CCCC";

        /// <summary>
        /// 塩基配列 AAAT
        /// </summary>
        private const string C_塩基配列_AAAT = "AAAT";

        /// <summary>
        /// 塩基配列 CCCT
        /// </summary>
        private const string C_塩基配列_CCCT = "CCCT";

        /// <summary>
        /// ファイル名 preprocessed 1 fq
        /// </summary>
        private const string C_ファイル名_preprocessed_1_fq = "preprocessed.1.fq";

        /// <summary>
        /// ファイル名 preprocessed 2 fq
        /// </summary>
        private const string C_ファイル名_preprocessed_2_fq = "preprocessed.2.fq";

        /// <summary>
        /// 塩基配列 TTTT
        /// </summary>
        private const string C_塩基配列_TTTT = "TTTT";

        /// <summary>
        /// 項目 upstream
        /// </summary>
        private const string C_項目_upstream = "upstream";

        /// <summary>
        /// ファイル名 sha256
        /// </summary>
        private const string C_ファイル名_sha256 = ".sha256";

        /// <summary>
        /// ファイル名 forward fq
        /// </summary>
        private const string C_ファイル名_forward_fq = "forward.fq";

        /// <summary>
        /// ファイル名 reverse fq
        /// </summary>
        private const string C_ファイル名_reverse_fq = "reverse.fq";

        /// <summary>
        /// 項目 input
        /// </summary>
        private const string C_項目_input = "input";

        /// <summary>
        /// 塩基配列 TTTA
        /// </summary>
        private const string C_塩基配列_TTTA = "TTTA";

        /// <summary>
        /// 項目 unknown
        /// </summary>
        private const string C_項目_unknown = "--unknown";

        /// <summary>
        /// 項目 invalid
        /// </summary>
        private const string C_項目_invalid = "invalid";

        /// <summary>
        /// ファイル名 notes txt
        /// </summary>
        private const string C_ファイル名_notes_txt = "notes.txt";

        /// <summary>
        /// 項目 user data
        /// </summary>
        private const string C_項目_user_data = "user-data";

        /// <summary>
        /// 項目 keep
        /// </summary>
        private const string C_項目_keep = "keep";

        #endregion

        #region 内部変数

        /// <summary>
        /// テスト専用の作業パス
        /// </summary>
        private readonly string _作業パス = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_checkpoint + Guid.NewGuid().ToString(C_GUID書式));

        #endregion

        #region コンストラクタ

        /// <summary>
        /// テスト用ディレクトリを用意する
        /// </summary>
        public StageCheckpointTests()
        {
            Directory.CreateDirectory(this._作業パス);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 最終レポートは過去の検査値を引き継がず原入力を再検査する
        /// </summary>
        [Fact]
        public void V_最終検査_古い合格を引き継がない()
        {
            var l_入力 = Path.Combine(this._作業パス, C_ファイル名_raw_fq);
            var l_配列 = Path.Combine(this._作業パス, C_ファイル名_source_fasta);
            var l_出力 = Path.Combine(this._作業パス, C_項目_output);
            Directory.CreateDirectory(l_出力);
            File.WriteAllText(l_入力, string.Concat(Enumerable.Range(0, 8).Select(i => $"@read{i}\n{new string('A', 250)}\n+\n{new string('I', 250)}\n")));
            File.WriteAllText(l_配列, $">contig\n{new string('C', 250)}\n");
            var l_原設定 = ConfigurationManager.A_実行時引数;
            var l_設定 = new Parameters
            {
                A_順リードのパス = l_入力,
                A_k長 = 31,
                A_スレッド数 = 1,
                A_コピー数基準の出所 = コピー数基準の出所.Weighted,
                A_Is低カバレッジ端トリミング = false
            };
            try
            {
                ConfigurationManager.A_実行時引数 = l_設定;
                var l_結果 = new アセンブリ実行結果(31, l_配列, l_配列, null, 2UL, 100D, A_整合性検査: new 整合性検査結果(1L, 1L, 1L, 0L, 0L, 0L), A_実際のコピー数基準: コピー数基準の出所.Weighted);
                FinalAssemblyPipeline.V_実行(l_結果, l_設定, l_出力, 100);
                using var l_レポート = JsonDocument.Parse(File.ReadAllText(Path.Combine(l_出力, Consts.レポートファイル名)));
                Assert.False(l_レポート.RootElement.GetProperty(C_項目_complete).GetBoolean());
                Assert.True(l_レポート.RootElement.GetProperty(C_項目_self_check).GetProperty(C_項目_missing_kmers).GetInt64() > 0L);
                Assert.Contains(l_レポート.RootElement.GetProperty(C_項目_checks).EnumerateArray(), x => x.GetProperty(C_項目_name).GetString() == C_項目_read_support && x.GetProperty(C_項目_result).GetString() == C_項目_fail);
                Assert.Equal(C_項目_Weighted, l_レポート.RootElement.GetProperty(C_項目_assembly_settings).GetProperty(C_項目_copy_number_baseline_requested).GetString());
                Assert.Equal(C_項目_Weighted, l_レポート.RootElement.GetProperty(C_項目_assembly_settings).GetProperty(C_項目_copy_number_baseline_actual).GetString());
                Assert.False(l_レポート.RootElement.GetProperty(C_項目_assembly_settings).GetProperty(C_項目_trim_low_coverage_ends).GetBoolean());
                using var l_出所 = JsonDocument.Parse(File.ReadAllText(Path.Combine(l_出力, C_ファイル名_assembly_provenance_json)));
                Assert.Equal(StageCheckpoint.Get_ハッシュ(Path.Combine(l_出力, C_ファイル名_assembly_fasta)), l_出所.RootElement.GetProperty(C_項目_assembly_sha256).GetString());
                Assert.Equal(C_項目_Weighted, l_出所.RootElement.GetProperty(C_項目_assembly_settings).GetProperty(C_項目_copy_number_baseline_requested).GetString());
                Assert.Equal(C_項目_Weighted, l_出所.RootElement.GetProperty(C_項目_assembly_settings).GetProperty(C_項目_copy_number_baseline_actual).GetString());
                Assert.False(l_出所.RootElement.GetProperty(C_項目_assembly_settings).GetProperty(C_項目_trim_low_coverage_ends).GetBoolean());
            }
            finally
            {
                ConfigurationManager.A_実行時引数 = l_原設定;
            }
        }

        /// <summary>
        /// テスト用ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            Directory.Delete(this._作業パス, recursive: true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// 同じ長さでも入力または出力が変われば再利用を拒否する
        /// </summary>
        [Fact]
        public void V_照合_内容変更を拒否()
        {
            var l_入力 = Path.Combine(this._作業パス, C_ファイル名_raw_fq);
            var l_出力 = Path.Combine(this._作業パス, C_ファイル名_corrected_fq);
            File.WriteAllText(l_入力, C_塩基配列_AAAA);
            File.WriteAllText(l_出力, C_塩基配列_CCCC);
            var l_設定 = new Parameters
            {
                A_順リードのパス = l_入力
            };
            var l_署名 = StageCheckpoint.Get_入力署名(l_設定);
            Assert.False(StageCheckpoint.Is再利用可能(l_署名, l_出力, null));
            StageCheckpoint.V_保存(l_署名, l_出力, null);
            Assert.True(StageCheckpoint.Is再利用可能(l_署名, l_出力, null));
            l_設定.A_Is再開 = true;
            Assert.Equal(l_署名, StageCheckpoint.Get_入力署名(l_設定));
            File.WriteAllText(l_入力, C_塩基配列_AAAT);
            Assert.False(StageCheckpoint.Is再利用可能(StageCheckpoint.Get_入力署名(l_設定), l_出力, null));
            File.WriteAllText(l_出力, C_塩基配列_CCCT);
            Assert.False(StageCheckpoint.Is再利用可能(l_署名, l_出力, null));
        }

        /// <summary>
        /// 前処理済みリードは実体だけ消し、再開に要る記録は残す
        /// </summary>
        [Fact]
        public void V_削除_前処理済みリードは記録を残して消す()
        {
            var l_前処理済み1 = Path.Combine(this._作業パス, C_ファイル名_preprocessed_1_fq);
            var l_前処理済み2 = Path.Combine(this._作業パス, C_ファイル名_preprocessed_2_fq);
            File.WriteAllText(l_前処理済み1, C_塩基配列_AAAA);
            File.WriteAllText(l_前処理済み2, C_塩基配列_TTTT);
            StageCheckpoint.V_保存(C_項目_upstream, l_前処理済み1, l_前処理済み2);
            ReadPreparationPipeline.V_削除_前処理済みリード(this._作業パス);
            Assert.False(File.Exists(l_前処理済み1));
            Assert.False(File.Exists(l_前処理済み2));
            Assert.True(File.Exists(l_前処理済み1 + C_ファイル名_sha256));
        }

        /// <summary>
        /// 上流工程の記録があれば、その出力の実体が無くても署名は変わらない
        /// </summary>
        [Fact]
        public void V_署名_上流の記録があれば入力の実体が無くても変わらない()
        {
            var l_上流1 = Path.Combine(this._作業パス, C_ファイル名_preprocessed_1_fq);
            var l_上流2 = Path.Combine(this._作業パス, C_ファイル名_preprocessed_2_fq);
            File.WriteAllText(l_上流1, C_塩基配列_AAAA);
            File.WriteAllText(l_上流2, C_塩基配列_TTTT);
            StageCheckpoint.V_保存(C_項目_upstream, l_上流1, l_上流2);
            var l_設定 = new Parameters
            {
                A_順リードのパス = l_上流1,
                A_逆リードのパス = l_上流2
            };
            var l_署名 = StageCheckpoint.Get_入力署名(l_設定);
            File.Delete(l_上流1);
            File.Delete(l_上流2);
            Assert.Equal(l_署名, StageCheckpoint.Get_入力署名(l_設定));
        }

        /// <summary>
        /// 上流工程の出力が作り直されて中身が変われば署名も変わる
        /// </summary>
        [Fact]
        public void V_署名_上流の出力が変われば変わる()
        {
            var l_上流 = Path.Combine(this._作業パス, C_ファイル名_preprocessed_1_fq);
            File.WriteAllText(l_上流, C_塩基配列_AAAA);
            StageCheckpoint.V_保存(C_項目_upstream, l_上流, null);
            var l_設定 = new Parameters
            {
                A_順リードのパス = l_上流
            };
            var l_署名 = StageCheckpoint.Get_入力署名(l_設定);
            File.WriteAllText(l_上流, C_塩基配列_AAAT);
            StageCheckpoint.V_保存(C_項目_upstream, l_上流, null);
            Assert.NotEqual(l_署名, StageCheckpoint.Get_入力署名(l_設定));
        }

        /// <summary>
        /// 片側の出力の欠落や破損を検出する
        /// </summary>
        [Fact]
        public void V_照合_対出力の破損を拒否()
        {
            var l_順鎖 = Path.Combine(this._作業パス, C_ファイル名_forward_fq);
            var l_逆鎖 = Path.Combine(this._作業パス, C_ファイル名_reverse_fq);
            File.WriteAllText(l_順鎖, C_塩基配列_AAAA);
            File.WriteAllText(l_逆鎖, C_塩基配列_TTTT);
            StageCheckpoint.V_保存(C_項目_input, l_順鎖, l_逆鎖);
            File.WriteAllText(l_逆鎖, C_塩基配列_TTTA);
            Assert.False(StageCheckpoint.Is再利用可能(C_項目_input, l_順鎖, l_逆鎖));
            File.Delete(l_逆鎖);
            Assert.False(StageCheckpoint.Is再利用可能(C_項目_input, l_順鎖, l_逆鎖));
        }

        /// <summary>
        /// 設定の複製は元の入力と明示指定を保つ
        /// </summary>
        [Fact]
        public void V_複製_原設定を保存()
        {
            var l_元 = new Parameters();
            l_元.Set_k長一覧([31, 63]);
            var l_複製 = l_元.Get_複製();
            l_複製.Set_k長一覧([21, 41]);
            Assert.Equal([31, 63], l_元.A_k長一覧);
            Assert.True(l_複製.A_Isk長明示指定);
        }

        /// <summary>
        /// 不正な引数を成功として扱わない
        /// </summary>
        [Fact]
        public void V_引数_エラー終了()
        {
            Assert.Throws<ArgumentException>(() => ArgumentsReader.Get_実行時引数([C_項目_unknown]));
            Assert.Equal(1, AssemblyApplication.Get_終了コード([Consts.引数キー.スレッド数, C_項目_invalid]));
            Assert.Equal(0, AssemblyApplication.Get_終了コード([Consts.引数キー.ヘルプ]));
        }

        /// <summary>
        /// 既存ディレクトリで失敗しても削除指定を実行しない
        /// </summary>
        [Fact]
        public void V_実行失敗_既存の入力を保護()
        {
            var l_入力 = Path.Combine(this._作業パス, C_ファイル名_raw_fq);
            File.WriteAllText(l_入力, $"@read\n{new string('A', 100)}\n+\n{new string('I', 100)}\n");
            var l_元 = File.ReadAllBytes(l_入力);
            Assert.Equal(1, AssemblyApplication.Get_終了コード([Consts.引数キー.リード1のパス, l_入力, Consts.引数キー.一時ディレクトリ, this._作業パス, Consts.引数キー.一時ディレクトリ削除, Consts.引数キー.再開]));
            Assert.Equal(l_元, File.ReadAllBytes(l_入力));
        }

        /// <summary>
        /// 中間成果物以外のファイルは削除しない
        /// </summary>
        [Fact]
        public void V_削除_未知の資料を保持()
        {
            var l_資料 = Path.Combine(this._作業パス, C_ファイル名_notes_txt);
            var l_データ = Path.Combine(this._作業パス, C_項目_user_data);
            File.WriteAllText(l_資料, C_項目_keep);
            Directory.CreateDirectory(l_データ);
            AssemblyWorkspace.V_削除_中間ファイル(this._作業パス);
            Assert.True(File.Exists(l_資料));
            Assert.True(Directory.Exists(l_データ));
        }

        #endregion
    }
}
