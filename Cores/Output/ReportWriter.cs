using System.Globalization;
using System.Text;
using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Polishing;
using Tsumiki.Models.Reporting;

namespace Tsumiki.Cores.Output
{
    /// <summary>
    /// 完全長の判定と、その根拠になった数値をファイルへ書き出す
    /// </summary>
    internal static class ReportWriter
    {
        #region 公開メソッド

        /// <summary>
        /// 完全性レポートを JSON で書き出す
        /// </summary>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_k長">採用した k の長さ</param>
        /// <param name="p_統計">アセンブリ統計</param>
        /// <param name="p_未解決ギャップ数">埋まらなかったギャップの数</param>
        /// <param name="p_環状本数">環状に閉じた配列の本数</param>
        /// <param name="p_判定">完全長の判定結果</param>
        /// <param name="p_整合性">自己検査の結果</param>
        /// <param name="p_閉鎖検証">環状閉鎖の検証結果</param>
        /// <param name="p_ポリッシュ">ポリッシュの結果</param>
        /// <param name="p_曖昧箇所">決めきれなかった箇所</param>
        public static void V_書き出し_レポート(string p_出力パス, int p_k長, アセンブリ統計 p_統計, int p_未解決ギャップ数, int p_環状本数, 完全性判定結果 p_判定, 整合性検査結果? p_整合性, IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証, ポリッシュ統計? p_ポリッシュ, IReadOnlyList<曖昧箇所> p_曖昧箇所)
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine("{");
            V_追加(l_文, "  \"tsumiki_version\": {0},", Get_文字列(Consts.バージョン));
            V_追加(l_文, "  \"complete\": {0},", p_判定.A_Is完全長 ? "true" : "false");
            V_追加(l_文, "  \"quality_level\": {0},", Get_文字列("Q" + (int)p_判定.A_品質保証レベル));
            V_追加(l_文, "  \"reason_codes\": [{0}],", string.Join(", ", p_判定.A_未達理由.Select(x => Get_文字列(CompletenessValidator.Get_理由コード(x)))));
            V_追加(l_文, "  \"k\": {0},", p_k長);
            V_追加(l_文, "  \"sequences\": {0},", p_統計.A_配列数);
            V_追加(l_文, "  \"total_length\": {0},", p_統計.A_総延長);
            V_追加(l_文, "  \"largest\": {0},", p_統計.A_最大長);
            V_追加(l_文, "  \"n50\": {0},", p_統計.A_N50);
            V_追加(l_文, "  \"l50\": {0},", p_統計.A_L50);
            V_追加(l_文, "  \"gc_percent\": {0},", Get_数値(p_統計.A_GC率));
            V_追加(l_文, "  \"circular_replicons\": {0},", p_環状本数);
            V_追加(l_文, "  \"unresolved_gaps\": {0},", p_未解決ギャップ数);

            _ = l_文.AppendLine("  \"checks\": [");
            for (var i = 0; i < p_判定.A_検査項目.Count; i++)
            {
                var l_項目 = p_判定.A_検査項目[i];
                var l_末尾 = i == p_判定.A_検査項目.Count - 1 ? string.Empty : ",";
                V_追加(l_文, "    {{\"name\": {0}, \"result\": {1}, \"detail\": {2}}}{3}", Get_文字列(l_項目.A_キー), Get_文字列(CompletenessValidator.Get_判定コード(l_項目.A_判定)), Get_文字列(l_項目.A_内訳), l_末尾);
            }
            _ = l_文.AppendLine("  ],");

            V_追加_自己検査(l_文, p_整合性);
            V_追加_ポリッシュ(l_文, p_ポリッシュ);
            V_追加_閉鎖検証(l_文, p_閉鎖検証);
            V_追加(l_文, "  \"ambiguous_junctions\": {0}", p_曖昧箇所.Count);

            _ = l_文.AppendLine("}");

            File.WriteAllText(p_出力パス, l_文.ToString());
        }

        /// <summary>
        /// リードに裏付けの無い区間を TSV で書き出す
        /// </summary>
        /// <remarks>
        /// 位置は 1 始まり・両端を含む
        /// </remarks>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_区間">支持のない区間</param>
        /// <param name="p_r長">支持を問うた r-mer の長さ</param>
        public static void V_書き出し_未支持箇所(string p_出力パス, IReadOnlyList<支持のない区間> p_区間, int p_r長)
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine(string.Join('\t', "sequence", "start", "end", "length", "r"));
            foreach (var l_区間 in p_区間)
            {
                _ = l_文.AppendLine(string.Join('\t', l_区間.A_配列ID, l_区間.A_開始, l_区間.A_終了, l_区間.A_長さ, p_r長));
            }
            File.WriteAllText(p_出力パス, l_文.ToString());
        }

        /// <summary>
        /// 決めきれなかった箇所を TSV で書き出す
        /// </summary>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_曖昧箇所">決めきれなかった箇所</param>
        public static void V_書き出し_曖昧箇所(string p_出力パス, IReadOnlyList<曖昧箇所> p_曖昧箇所)
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine(string.Join('\t', "k", "type", "location", "top_support", "second_support", "margin", "raw_support", "confidence"));
            foreach (var l_箇所 in p_曖昧箇所)
            {
                _ = l_文.AppendLine(string.Join('\t', l_箇所.A_k長, Get_種別コード(l_箇所.A_種別), l_箇所.A_場所, Get_数値(l_箇所.A_首位の支持), Get_数値(l_箇所.A_次点の支持), Get_数値(l_箇所.A_余裕), l_箇所.A_首位の生支持数, Get_数値(l_箇所.A_確信度)));
            }
            File.WriteAllText(p_出力パス, l_文.ToString());
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// TSV に出す固定の種別名
        /// </summary>
        /// <param name="p_種別"></param>
        /// <remarks>
        /// 訳さない
        /// </remarks>
        /// <returns></returns>
        private static string Get_種別コード(曖昧箇所の種別 p_種別)
        {
            return p_種別 switch
            {
                曖昧箇所の種別.支持なし => "no-support",
                曖昧箇所の種別.僅差 => "insufficient-margin",
                曖昧箇所の種別.経路が一意でない => "multiple-paths",
                曖昧箇所の種別.到達不能 => "unreachable",
                _ => "inside-repeat",
            };
        }

        /// <summary>
        /// 自己検査の結果をレポートへ足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_整合性">自己検査の結果</param>
        private static void V_追加_自己検査(StringBuilder p_文, 整合性検査結果? p_整合性)
        {
            if (p_整合性 is not { } l_整合性)
            {
                _ = p_文.AppendLine("  \"self_check\": null,");
                return;
            }
            _ = p_文.AppendLine("  \"self_check\": {");
            V_追加(p_文, "    \"trusted_kmers\": {0},", l_整合性.A_信頼kmer数);
            V_追加(p_文, "    \"missing_kmers\": {0},", l_整合性.A_取りこぼし数);
            V_追加(p_文, "    \"missing_percent\": {0},", Get_数値(l_整合性.A_取りこぼし率));
            V_追加(p_文, "    \"excess_percent\": {0}", Get_数値(l_整合性.A_出しすぎ率));
            _ = p_文.AppendLine("  },");
        }

        /// <summary>
        /// ポリッシュの結果をレポートへ足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_ポリッシュ">ポリッシュの結果</param>
        private static void V_追加_ポリッシュ(StringBuilder p_文, ポリッシュ統計? p_ポリッシュ)
        {
            if (p_ポリッシュ is not { } l_ポリッシュ)
            {
                _ = p_文.AppendLine("  \"polish\": null,");
                return;
            }
            _ = p_文.AppendLine("  \"polish\": {");
            V_追加(p_文, "    \"mapped_reads\": {0},", l_ポリッシュ.A_マップされたリード数);
            V_追加(p_文, "    \"rejected_reads\": {0},", l_ポリッシュ.A_棄却されたリード数);
            V_追加(p_文, "    \"corrected_bases\": {0},", l_ポリッシュ.A_訂正した塩基数);
            V_追加(p_文, "    \"median_depth\": {0},", Get_数値(l_ポリッシュ.A_深度の中央値));
            V_追加(p_文, "    \"low_depth_fraction\": {0}", Get_数値(l_ポリッシュ.A_深度不足率));
            _ = p_文.AppendLine("  },");
        }

        /// <summary>
        /// 環状閉鎖の検証結果をレポートへ足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_閉鎖検証">環状閉鎖の検証結果</param>
        private static void V_追加_閉鎖検証(StringBuilder p_文, IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証)
        {
            if (p_閉鎖検証 is null)
            {
                _ = p_文.AppendLine("  \"circular_closure\": null,");
                return;
            }
            _ = p_文.AppendLine("  \"circular_closure\": [");
            for (var i = 0; i < p_閉鎖検証.Count; i++)
            {
                var l_検証 = p_閉鎖検証[i];
                var l_末尾 = i == p_閉鎖検証.Count - 1 ? string.Empty : ",";
                V_追加(p_文, "    {{\"id\": {0}, \"length\": {1}, \"spanning_reads\": {2}, \"required\": {3}, \"supported\": {4}}}{5}", Get_文字列(l_検証.A_配列ID), l_検証.A_長さ, l_検証.A_跨いだリード数, l_検証.A_必要本数, l_検証.A_Has支持 ? "true" : "false", l_末尾);
            }
            _ = p_文.AppendLine("  ],");
        }

        /// <summary>
        /// レポートへ 1 行足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_書式">行の書式</param>
        /// <param name="p_引数">書式へ埋める値</param>
        private static void V_追加(StringBuilder p_文, string p_書式, params object?[] p_引数)
        {
            _ = p_文.AppendLine(string.Format(CultureInfo.InvariantCulture, p_書式, p_引数));
        }

        /// <summary>
        /// JSON へ書ける形の数値にして返す
        /// </summary>
        /// <param name="p_値">元の値</param>
        /// <returns>JSON へ書ける文字列</returns>
        private static string Get_数値(double p_値)
        {
            return p_値.ToString("0.####", CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// JSON の文字列リテラル
        /// </summary>
        /// <param name="p_値"></param>
        /// <remarks>
        /// ID には引用符も含まれうる
        /// </remarks>
        /// <returns></returns>
        private static string Get_文字列(string? p_値)
        {
            if (p_値 is null)
            {
                return "null";
            }
            var l_文 = new StringBuilder("\"");
            foreach (var l_文字 in p_値)
            {
                _ = l_文字 switch
                {
                    '"' => l_文.Append("\\\""),
                    '\\' => l_文.Append("\\\\"),
                    '\n' => l_文.Append("\\n"),
                    '\r' => l_文.Append("\\r"),
                    '\t' => l_文.Append("\\t"),
                    _ => l_文字 < ' ' ? l_文.Append(CultureInfo.InvariantCulture, $"\\u{(int)l_文字:x4}") : l_文.Append(l_文字),
                };
            }
            return l_文.Append('"').ToString();
        }

        #endregion
    }
}
