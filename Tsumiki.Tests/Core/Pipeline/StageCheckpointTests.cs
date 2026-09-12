using Tsumiki.Cores.Pipeline;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Evaluation;
using Tsumiki.Commons;
using System.Text.Json;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 再開時の内容照合と実行境界を検証する
    /// </summary>
    public class StageCheckpointTests : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// テスト専用の作業パス
        /// </summary>
        private readonly string _作業パス = Path.Combine(Path.GetTempPath(), "tsumiki_checkpoint_" + Guid.NewGuid().ToString("N"));

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
            var l_入力 = Path.Combine(this._作業パス, "raw.fq");
            var l_配列 = Path.Combine(this._作業パス, "source.fasta");
            var l_出力 = Path.Combine(this._作業パス, "output");
            Directory.CreateDirectory(l_出力);
            File.WriteAllText(l_入力, string.Concat(Enumerable.Range(0, 8).Select(i => $"@read{i}\n{new string('A', 250)}\n+\n{new string('I', 250)}\n")));
            File.WriteAllText(l_配列, $">contig\n{new string('C', 250)}\n");
            var l_原設定 = ConfigurationManager.A_実行時引数;
            var l_設定 = new Parameters { A_リード1のパス = l_入力, A_k長 = 31, A_スレッド数 = 1 };
            try
            {
                ConfigurationManager.A_実行時引数 = l_設定;
                var l_結果 = new アセンブリ実行結果(31, l_配列, l_配列, null, 2UL, 100D, A_整合性検査: new 整合性検査結果(1L, 1L, 1L, 0L, 0L, 0L));
                FinalAssemblyPipeline.V_実行(l_結果, l_設定, l_出力, 100);
                using var l_レポート = JsonDocument.Parse(File.ReadAllText(Path.Combine(l_出力, "assembly.report.json")));
                Assert.False(l_レポート.RootElement.GetProperty("complete").GetBoolean());
                Assert.True(l_レポート.RootElement.GetProperty("self_check").GetProperty("missing_kmers").GetInt64() > 0L);
                Assert.Contains(l_レポート.RootElement.GetProperty("checks").EnumerateArray(), x => x.GetProperty("name").GetString() == "read_support" && x.GetProperty("result").GetString() == "fail");
                using var l_出所 = JsonDocument.Parse(File.ReadAllText(Path.Combine(l_出力, "assembly.provenance.json")));
                Assert.Equal(StageCheckpoint.Get_ハッシュ(Path.Combine(l_出力, "assembly.fasta")), l_出所.RootElement.GetProperty("assembly_sha256").GetString());
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
            var l_入力 = Path.Combine(this._作業パス, "raw.fq");
            var l_出力 = Path.Combine(this._作業パス, "corrected.fq");
            File.WriteAllText(l_入力, "AAAA");
            File.WriteAllText(l_出力, "CCCC");
            var l_設定 = new Parameters { A_リード1のパス = l_入力 };
            var l_署名 = StageCheckpoint.Get_入力署名(l_設定);
            Assert.False(StageCheckpoint.Get_再利用可能(l_署名, l_出力, null));
            StageCheckpoint.V_保存(l_署名, l_出力, null);
            Assert.True(StageCheckpoint.Get_再利用可能(l_署名, l_出力, null));
            l_設定.A_再開するか = true;
            Assert.Equal(l_署名, StageCheckpoint.Get_入力署名(l_設定));
            File.WriteAllText(l_入力, "AAAT");
            Assert.False(StageCheckpoint.Get_再利用可能(StageCheckpoint.Get_入力署名(l_設定), l_出力, null));
            File.WriteAllText(l_出力, "CCCT");
            Assert.False(StageCheckpoint.Get_再利用可能(l_署名, l_出力, null));
        }

        /// <summary>
        /// 片側の出力の欠落や破損を検出する
        /// </summary>
        [Fact]
        public void V_照合_対出力の破損を拒否()
        {
            var l_順鎖 = Path.Combine(this._作業パス, "forward.fq");
            var l_逆鎖 = Path.Combine(this._作業パス, "reverse.fq");
            File.WriteAllText(l_順鎖, "AAAA");
            File.WriteAllText(l_逆鎖, "TTTT");
            StageCheckpoint.V_保存("input", l_順鎖, l_逆鎖);
            File.WriteAllText(l_逆鎖, "TTTA");
            Assert.False(StageCheckpoint.Get_再利用可能("input", l_順鎖, l_逆鎖));
            File.Delete(l_逆鎖);
            Assert.False(StageCheckpoint.Get_再利用可能("input", l_順鎖, l_逆鎖));
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
            Assert.True(l_複製.A_k長が明示指定されたか);
        }

        /// <summary>
        /// 不正な引数を成功として扱わない
        /// </summary>
        [Fact]
        public void V_引数_エラー終了()
        {
            Assert.Throws<ArgumentException>(() => ArgumentsReader.Get_実行時引数(["--unknown"]));
            Assert.Equal(1, AssemblyApplication.Get_終了コード(["-th", "invalid"]));
            Assert.Equal(0, AssemblyApplication.Get_終了コード(["-h"]));
        }

        /// <summary>
        /// 既存ディレクトリで失敗しても削除指定を実行しない
        /// </summary>
        [Fact]
        public void V_実行失敗_既存の入力を保護()
        {
            var l_入力 = Path.Combine(this._作業パス, "raw.fq");
            File.WriteAllText(l_入力, $"@read\n{new string('A', 100)}\n+\n{new string('I', 100)}\n");
            var l_元 = File.ReadAllBytes(l_入力);
            Assert.Equal(1, AssemblyApplication.Get_終了コード(["-1", l_入力, "-t", this._作業パス, "-rt", "-rs"]));
            Assert.Equal(l_元, File.ReadAllBytes(l_入力));
        }

        /// <summary>
        /// 中間成果物以外のファイルは削除しない
        /// </summary>
        [Fact]
        public void V_削除_未知の資料を保持()
        {
            var l_資料 = Path.Combine(this._作業パス, "notes.txt");
            var l_データ = Path.Combine(this._作業パス, "user-data");
            File.WriteAllText(l_資料, "keep");
            Directory.CreateDirectory(l_データ);
            AssemblyWorkspace.V_削除_中間ファイル(this._作業パス);
            Assert.True(File.Exists(l_資料));
            Assert.True(Directory.Exists(l_データ));
        }

        #endregion
    }
}
