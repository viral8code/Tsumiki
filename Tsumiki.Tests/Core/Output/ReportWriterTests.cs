using System.Text.Json;
using Tsumiki.Core.Evaluation;
using Tsumiki.Core.Output;
using Tsumiki.Core;
using Tsumiki.Model.Evaluation;
using Tsumiki.Model.Polishing;
using Tsumiki.Model.Reporting;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// レポートの書き出しを固定する。JSON は手で組み立てているので、
    /// 実際に構文として通ることと、判定の要点が載っていることを確かめる。
    /// </summary>
    public class ReportWriterTests : IDisposable
    {
        private readonly string _一時ディレクトリ;

        public ReportWriterTests()
        {
            this._一時ディレクトリ = Path.Combine(
                Path.GetTempPath(), "tsumiki_report_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._一時ディレクトリ);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._一時ディレクトリ))
            {
                Directory.Delete(this._一時ディレクトリ, recursive: true);
            }
            GC.SuppressFinalize(this);
        }

        /// <summary>リードに裏付けの無い位置が一つも無い検査結果。</summary>
        private static 支持検査結果 Get_良好な支持()
        {
            return new 支持検査結果(A_r長: 31, A_調べた位置数: 100000, A_支持のない位置数: 0, A_区間: []);
        }

        private static アセンブリ統計 Get_統計()
        {
            return new アセンブリ統計(3, 5_000_000, 4_800_000, 900, 4_800_000, 1, 50.5);
        }

        private JsonElement Get_書き出したJSON(
            完全性判定結果 p_判定,
            整合性検査結果? p_整合性 = null,
            IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証 = null,
            ポリッシュ統計? p_ポリッシュ = null,
            IReadOnlyList<曖昧箇所>? p_曖昧箇所 = null)
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, "assembly.report.json");
            ReportWriter.V_書き出し_レポート(
                l_パス, 63, Get_統計(), 2, 1, p_判定,
                p_整合性, p_閉鎖検証, p_ポリッシュ, p_曖昧箇所 ?? []);
            return JsonDocument.Parse(File.ReadAllText(l_パス)).RootElement.Clone();
        }

        [Fact]
        public void V_書き出し_レポート_完全長の判定をそのまま載せる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 0,
                p_整合性: new 整合性検査結果(1000, 1000, 1000, 10, 0, 0),
                p_閉鎖検証: [new 環状閉鎖検証結果("scaffold1_circular", 4_800_000, 30, 5)],
                p_ポリッシュ: new ポリッシュ統計(3, 5_000_000, 1000, 0, 12, 0, 5_000_000, 90),
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            var l_JSON = this.Get_書き出したJSON(
                l_判定,
                new 整合性検査結果(1000, 1000, 1000, 10, 0, 0),
                [new 環状閉鎖検証結果("scaffold1_circular", 4_800_000, 30, 5)],
                new ポリッシュ統計(3, 5_000_000, 1000, 0, 12, 0, 5_000_000, 90));

            Assert.True(l_JSON.GetProperty("complete").GetBoolean());
            Assert.Equal("Q5", l_JSON.GetProperty("quality_level").GetString());
            Assert.Equal(0, l_JSON.GetProperty("reason_codes").GetArrayLength());
            Assert.Equal(63, l_JSON.GetProperty("k").GetInt32());
            Assert.Equal(2, l_JSON.GetProperty("unresolved_gaps").GetInt32());
            Assert.Equal(1, l_JSON.GetProperty("circular_replicons").GetInt32());
            Assert.Equal(12, l_JSON.GetProperty("polish").GetProperty("corrected_bases").GetInt32());
            Assert.Equal(
                30,
                l_JSON.GetProperty("circular_closure")[0].GetProperty("spanning_reads").GetInt32());
        }

        [Fact]
        public void V_書き出し_レポート_未達の理由をコードで載せる()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(
                p_未解決ギャップ数: 2,
                p_整合性: new 整合性検査結果(1000, 1000, 1000, 10, 0, 0),
                p_閉鎖検証: null,
                p_ポリッシュ: null,
                p_曖昧箇所: [],
                p_支持検査: Get_良好な支持());

            var l_JSON = this.Get_書き出したJSON(
                l_判定, new 整合性検査結果(1000, 1000, 1000, 10, 0, 0));

            Assert.False(l_JSON.GetProperty("complete").GetBoolean());

            var l_理由 = l_JSON.GetProperty("reason_codes")
                .EnumerateArray().Select(x => x.GetString()).ToList();
            Assert.Contains("unresolved-gap", l_理由);
            Assert.Contains("depth-unmeasured", l_理由);
            Assert.Contains("closure-unverified", l_理由);

            // 測っていない項目は null として区別できる形で出す。
            Assert.Equal(JsonValueKind.Null, l_JSON.GetProperty("polish").ValueKind);
            Assert.Equal(JsonValueKind.Null, l_JSON.GetProperty("circular_closure").ValueKind);
        }

        [Fact]
        public void V_書き出し_レポート_引用符を含むIDでも壊れない()
        {
            var l_判定 = CompletenessValidator.Get_判定結果(0, null, null, null, [], null);
            var l_JSON = this.Get_書き出したJSON(
                l_判定, null, [new 環状閉鎖検証結果("seq\"with\\quotes", 100, 1, 5)]);

            Assert.Equal(
                "seq\"with\\quotes",
                l_JSON.GetProperty("circular_closure")[0].GetProperty("id").GetString());
        }

        [Fact]
        public void V_書き出し_曖昧箇所_見出しと各行を出す()
        {
            var l_パス = Path.Combine(this._一時ディレクトリ, "assembly.ambiguous.tsv");
            ReportWriter.V_書き出し_曖昧箇所(
                l_パス,
                [
                    new 曖昧箇所(63, 曖昧箇所の種別.僅差, "unitig7+", 1.5, 1.4, 9, 0.95),
                    new 曖昧箇所(63, 曖昧箇所の種別.到達不能, "scaffold1:100-140", 0, 0, 0, 0),
                ]);

            var l_行 = File.ReadAllLines(l_パス);
            Assert.Equal(3, l_行.Length);
            Assert.StartsWith("k\ttype\tlocation", l_行[0]);
            Assert.Contains("insufficient-margin", l_行[1]);
            Assert.Contains("unitig7+", l_行[1]);
            Assert.Contains("unreachable", l_行[2]);

            // 列数は見出しと一致していなければ表として読めない。
            var l_列数 = l_行[0].Split('\t').Length;
            Assert.All(l_行, x => Assert.Equal(l_列数, x.Split('\t').Length));
        }
    }
}
