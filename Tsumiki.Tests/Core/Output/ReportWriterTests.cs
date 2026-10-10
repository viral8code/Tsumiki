using Tsumiki.Commons;
using System.Text.Json;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Output;
using Tsumiki.Core;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Polishing;
using Tsumiki.Models.Reporting;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// レポートの書き出しを固定する
    /// </summary>
    public class ReportWriterTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki report
        /// </summary>
        private const string C_項目_tsumiki_report = "tsumiki_report_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 項目 scaffold1 circular
        /// </summary>
        private const string C_項目_scaffold1_circular = "scaffold1_circular";

        /// <summary>
        /// 項目 complete
        /// </summary>
        private const string C_項目_complete = "complete";

        /// <summary>
        /// 項目 Q5
        /// </summary>
        private const string C_項目_Q5 = "Q5";

        /// <summary>
        /// 項目 quality level
        /// </summary>
        private const string C_項目_quality_level = "quality_level";

        /// <summary>
        /// 項目 reason codes
        /// </summary>
        private const string C_項目_reason_codes = "reason_codes";

        /// <summary>
        /// 項目 k
        /// </summary>
        private const string C_項目_k = "k";

        /// <summary>
        /// 項目 unresolved gaps
        /// </summary>
        private const string C_項目_unresolved_gaps = "unresolved_gaps";

        /// <summary>
        /// 項目 circular replicons
        /// </summary>
        private const string C_項目_circular_replicons = "circular_replicons";

        /// <summary>
        /// 項目 polish
        /// </summary>
        private const string C_項目_polish = "polish";

        /// <summary>
        /// 項目 corrected bases
        /// </summary>
        private const string C_項目_corrected_bases = "corrected_bases";

        /// <summary>
        /// 項目 circular closure
        /// </summary>
        private const string C_項目_circular_closure = "circular_closure";

        /// <summary>
        /// 項目 spanning reads
        /// </summary>
        private const string C_項目_spanning_reads = "spanning_reads";

        /// <summary>
        /// 項目 n split contig stats
        /// </summary>
        private const string C_項目_n_split_contig_stats = "n_split_contig_stats";

        /// <summary>
        /// 項目 minimum sequence length
        /// </summary>
        private const string C_項目_minimum_sequence_length = "minimum_sequence_length";

        /// <summary>
        /// 項目 scaffold stats
        /// </summary>
        private const string C_項目_scaffold_stats = "scaffold_stats";

        /// <summary>
        /// 項目 sequences
        /// </summary>
        private const string C_項目_sequences = "sequences";

        /// <summary>
        /// 項目 unresolved gap
        /// </summary>
        private const string C_項目_unresolved_gap = "unresolved-gap";

        /// <summary>
        /// 項目 depth unmeasured
        /// </summary>
        private const string C_項目_depth_unmeasured = "depth-unmeasured";

        /// <summary>
        /// 項目 closure unverified
        /// </summary>
        private const string C_項目_closure_unverified = "closure-unverified";

        /// <summary>
        /// 項目 seq with quotes
        /// </summary>
        private const string C_項目_seq_with_quotes = "seq\"with\\quotes";

        /// <summary>
        /// 項目 id
        /// </summary>
        private const string C_項目_id = "id";

        /// <summary>
        /// 項目 unitig7
        /// </summary>
        private const string C_項目_unitig7 = "unitig7+";

        /// <summary>
        /// 項目 scaffold1 100 140
        /// </summary>
        private const string C_項目_scaffold1_100_140 = "scaffold1:100-140";

        /// <summary>
        /// 項目 k type location
        /// </summary>
        private const string C_項目_k_type_location = "k\ttype\tlocation";

        /// <summary>
        /// 項目 insufficient margin
        /// </summary>
        private const string C_項目_insufficient_margin = "insufficient-margin";

        /// <summary>
        /// 項目 unreachable
        /// </summary>
        private const string C_項目_unreachable = "unreachable";

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public ReportWriterTests()
        {
            this._一時ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_report + Guid.NewGuid().ToString(C_GUID書式));
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

            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// 完全長の判定をそのまま JSON へ載せることを確かめる
        /// </summary>
        [Fact]
        public void V_書き出し_レポート_完全長の判定をそのまま載せる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 0, p_整合性: new 整合性検査結果(1_000L, 1_000L, 1_000L, 10L, 0L, 0L), p_閉鎖検証: [new 環状閉鎖検証結果(C_項目_scaffold1_circular, 4_800_000, 30, 5)], p_ポリッシュ: new ポリッシュ統計(3, 5_000_000L, 1_000L, 0L, 12L, 0L, 5_000_000L, 90D), p_曖昧箇所: [], p_支持検査: Get_良好な支持());
            var l_JSON = this.Get_書き出したJSON(l_判定, new 整合性検査結果(1_000L, 1_000L, 1_000L, 10L, 0L, 0L), [new 環状閉鎖検証結果(C_項目_scaffold1_circular, 4_800_000, 30, 5)], new ポリッシュ統計(3, 5_000_000L, 1_000L, 0L, 12L, 0L, 5_000_000L, 90D));
            Assert.True(l_JSON.GetProperty(C_項目_complete).GetBoolean());
            Assert.Equal(C_項目_Q5, l_JSON.GetProperty(C_項目_quality_level).GetString());
            Assert.Equal(0, l_JSON.GetProperty(C_項目_reason_codes).GetArrayLength());
            Assert.Equal(63, l_JSON.GetProperty(C_項目_k).GetInt32());
            Assert.Equal(2, l_JSON.GetProperty(C_項目_unresolved_gaps).GetInt32());
            Assert.Equal(1, l_JSON.GetProperty(C_項目_circular_replicons).GetInt32());
            Assert.Equal(12, l_JSON.GetProperty(C_項目_polish).GetProperty(C_項目_corrected_bases).GetInt32());
            Assert.Equal(30, l_JSON.GetProperty(C_項目_circular_closure)[0].GetProperty(C_項目_spanning_reads).GetInt32());
            Assert.Equal(500, l_JSON.GetProperty(C_項目_n_split_contig_stats).GetProperty(C_項目_minimum_sequence_length).GetInt32());
            Assert.Equal(3, l_JSON.GetProperty(C_項目_scaffold_stats).GetProperty(C_項目_sequences).GetInt32());
        }

        /// <summary>
        /// 未達の理由をコードで載せ、測っていない項目は null で区別することを確かめる
        /// </summary>
        [Fact]
        public void V_書き出し_レポート_未達の理由をコードで載せる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(p_未解決ギャップ数: 2, p_整合性: new 整合性検査結果(1_000L, 1_000L, 1_000L, 10L, 0L, 0L), p_閉鎖検証: null, p_ポリッシュ: null, p_曖昧箇所: [], p_支持検査: Get_良好な支持());
            var l_JSON = this.Get_書き出したJSON(l_判定, new 整合性検査結果(1_000L, 1_000L, 1_000L, 10L, 0L, 0L));
            Assert.False(l_JSON.GetProperty(C_項目_complete).GetBoolean());
            var l_理由 = l_JSON.GetProperty(C_項目_reason_codes).EnumerateArray().Select(x => x.GetString()).ToList();
            Assert.Contains(C_項目_unresolved_gap, l_理由);
            Assert.Contains(C_項目_depth_unmeasured, l_理由);
            Assert.Contains(C_項目_closure_unverified, l_理由);
            Assert.Equal(JsonValueKind.Null, l_JSON.GetProperty(C_項目_polish).ValueKind);
            Assert.Equal(JsonValueKind.Null, l_JSON.GetProperty(C_項目_circular_closure).ValueKind);
        }

        /// <summary>
        /// 引用符を含む ID でも JSON が壊れないことを確かめる
        /// </summary>
        [Fact]
        public void V_書き出し_レポート_引用符を含むIDでも壊れない()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(0, null, null, null, [], null);
            var l_JSON = this.Get_書き出したJSON(l_判定, null, [new 環状閉鎖検証結果(C_項目_seq_with_quotes, 100, 1, 5)]);
            Assert.Equal(C_項目_seq_with_quotes, l_JSON.GetProperty(C_項目_circular_closure)[0].GetProperty(C_項目_id).GetString());
        }

        /// <summary>
        /// 見出しと各行を出力し、列数が見出しと一致することを確かめる
        /// </summary>
        [Fact]
        public void V_書き出し_曖昧箇所_見出しと各行を出す()
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, Consts.曖昧箇所ファイル名);
            ReportWriter.V_書き出し_曖昧箇所(l_パス, [
                new 曖昧箇所(63, 曖昧箇所の種別.僅差, C_項目_unitig7, 1.5D, 1.4D, 9L, 0.95D),
                new 曖昧箇所(63, 曖昧箇所の種別.到達不能, C_項目_scaffold1_100_140, 0D, 0D, 0L, 0D),
            ]);
            var l_行 = File.ReadAllLines(l_パス);
            Assert.Equal(3, l_行.Length);
            Assert.StartsWith(C_項目_k_type_location, l_行[0]);
            Assert.Contains(C_項目_insufficient_margin, l_行[1]);
            Assert.Contains(C_項目_unitig7, l_行[1]);
            Assert.Contains(C_項目_unreachable, l_行[2]);
            var l_列数 = l_行[0].Split('\t').Length;
            Assert.All(l_行, x => Assert.Equal(l_列数, x.Split('\t').Length));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// リードに裏付けの無い位置が一つも無い検査結果
        /// </summary>
        /// <returns>支持検査の結果</returns>
        private static 支持検査結果 Get_良好な支持()
        {
            return new 支持検査結果(A_r長: 31, A_調べた位置数: 100_000L, A_支持のない位置数: 0L, A_区間: []);
        }

        /// <summary>
        /// 検証に使うアセンブリ統計
        /// </summary>
        /// <returns>アセンブリ統計</returns>
        private static アセンブリ統計 Get_統計()
        {
            return new アセンブリ統計(3, 5_000_000L, 4_800_000, 900, 4_800_000, 1, 50.5D);
        }

        /// <summary>
        /// レポートを書き出し、その JSON を読み直して返す
        /// </summary>
        /// <param name="p_判定">完全長の判定結果</param>
        /// <param name="p_整合性">自己検査の結果</param>
        /// <param name="p_閉鎖検証">環状閉鎖の検証結果</param>
        /// <param name="p_ポリッシュ">ポリッシュの結果</param>
        /// <param name="p_曖昧箇所"></param>
        /// <returns>書き出した JSON</returns>
        private JsonElement Get_書き出したJSON(完全性判定結果 p_判定, 整合性検査結果? p_整合性 = null, IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証 = null, ポリッシュ統計? p_ポリッシュ = null, IReadOnlyList<曖昧箇所>? p_曖昧箇所 = null)
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, Consts.レポートファイル名);
            ReportWriter.V_書き出し_レポート(l_パス, 63, Get_統計(), 2, 1, p_判定, p_整合性, p_閉鎖検証, p_ポリッシュ, p_曖昧箇所 ?? []);
            return JsonDocument.Parse(File.ReadAllText(l_パス)).RootElement.Clone();
        }

        #endregion
    }
}
