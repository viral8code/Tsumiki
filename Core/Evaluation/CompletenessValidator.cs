using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model.Evaluation;
using Tsumiki.Model.Foundation;
using Tsumiki.Model.Polishing;
using Tsumiki.Model.Reporting;

namespace Tsumiki.Core.Evaluation
{
    /// <summary>
    /// 最終成果物が完全長を名乗れるかを判定する
    /// </summary>
    /// <remarks>
    /// 完全長は「最長の配列がゲノムサイズに近い」ことではない<br/>
    /// 必要な検査を
    /// すべて通ったことを指し、材料が足りない項目は不合格ではなく判定不能として
    /// 区別する<br/>
    /// 情報が足りないところを推測で埋めて完全長を名乗らせないための
    /// 仕組みであり、判定できないことが分かる状態のほうが下流にとって安全
    /// </remarks>
    internal static class CompletenessValidator
    {
        /// <summary>
        /// 信頼できる k-mer の取りこぼしとして許す割合 (%)
        /// </summary>
        private const double 取りこぼしの許容率 = 5.0D;

        /// <summary>
        /// コピー数の推定を超えて出している延べ数として許す割合 (%)
        /// </summary>
        private const double 出しすぎの許容率 = 1.0D;

        /// <summary>
        /// 深度が落ち込んだ位置として許す割合
        /// </summary>
        private const double 深度不足の許容率 = 0.01D;

        /// <summary>
        /// リードに裏付けの無い位置として許す数
        /// </summary>
        /// <remarks>
        /// 割合ではなく数で見るのは、この検査が「そう繋いだ読みが一つも無い」
        /// という白黒のはっきりした事実を数えているため<br/>
        /// 総延長で薄めると
        /// 数箇所の捏造が見えなくなる
        /// </remarks>
        private const int 支持のない位置の許容数 = 0;

        /// <summary>
        /// 集めた材料から完全長かどうかを判定する
        /// </summary>
        /// <remarks>
        /// p_閉鎖検証 が null なら閉じ目を調べていない、p_ポリッシュ が null なら
        /// 深度を測っていないことを意味し、いずれも判定不能として扱う
        /// </remarks>
        public static 完全性判定結果 Get_判定結果(
            int p_未解決ギャップ数,
            整合性検査結果? p_整合性,
            IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証,
            ポリッシュ統計? p_ポリッシュ,
            IReadOnlyList<曖昧箇所> p_曖昧箇所,
            支持検査結果? p_支持検査)
        {
            var l_僅差の数 = p_曖昧箇所.Count(x => x.A_種別 == 曖昧箇所の種別.僅差);

            List<検査項目> l_項目 = [];
            List<未達理由> l_理由 = [];

            var l_取りこぼし = p_整合性 is { } l_整合1
                ? Get_判定(l_整合1.A_取りこぼし率 <= 取りこぼしの許容率, 未達理由.取りこぼしが多い, l_理由)
                : Get_判定不能(未達理由.自己検査を行えなかった, l_理由);
            l_項目.Add(new 検査項目(
                "graph_coverage", メッセージID.検査項目_グラフ被覆, l_取りこぼし,
                p_整合性 is { } l_整合2
                    ? $"{l_整合2.A_取りこぼし率:F2}% <= {取りこぼしの許容率:F2}%"
                    : string.Empty));

            var l_出しすぎ = p_整合性 is { } l_整合3
                ? Get_判定(l_整合3.A_出しすぎ率 <= 出しすぎの許容率, 未達理由.出しすぎている, l_理由)
                : 検査判定.判定不能;
            l_項目.Add(new 検査項目(
                "copy_consistency", メッセージID.検査項目_コピー数整合, l_出しすぎ,
                p_整合性 is { } l_整合4
                    ? $"{l_整合4.A_出しすぎ率:F2}% <= {出しすぎの許容率:F2}%"
                    : string.Empty));

            var l_深度 = p_ポリッシュ is { } l_ポリッシュ1
                ? Get_判定(l_ポリッシュ1.A_深度不足率 <= 深度不足の許容率, 未達理由.深度が不連続, l_理由)
                : Get_判定不能(未達理由.深度を測っていない, l_理由);
            l_項目.Add(new 検査項目(
                "coverage_continuity", メッセージID.検査項目_深度の連続性, l_深度,
                p_ポリッシュ is { } l_ポリッシュ2
                    ? $"{l_ポリッシュ2.A_深度不足率 * 100:F2}% <= {深度不足の許容率 * 100:F2}%"
                    : string.Empty));

            var l_ギャップ = Get_判定(p_未解決ギャップ数 == 0, 未達理由.未解決のギャップが残る, l_理由);
            l_項目.Add(new 検査項目(
                "unsupported_join", メッセージID.検査項目_未解決のギャップ, l_ギャップ,
                p_未解決ギャップ数.ToString()));

            var l_接合点 = Get_判定(p_曖昧箇所.Count == 0, 未達理由.決めきれない分岐が残る, l_理由);
            l_項目.Add(new 検査項目(
                "junction_support", メッセージID.検査項目_接合点の支持, l_接合点,
                p_曖昧箇所.Count.ToString()));

            // 僅差で捨てた箇所が残っているということは、同程度の証拠を持つ
            // 別の経路が残っているということ
            var l_代替経路 = Get_判定(l_僅差の数 == 0, 未達理由.決めきれない分岐が残る, l_理由);
            l_項目.Add(new 検査項目(
                "no_alternative_path", メッセージID.検査項目_競合経路, l_代替経路,
                l_僅差の数.ToString()));

            var l_支持 = p_支持検査 is { } l_支持1
                ? Get_判定(
                    l_支持1.A_支持のない位置数 <= 支持のない位置の許容数,
                    未達理由.リードに裏付けの無い箇所がある, l_理由)
                : Get_判定不能(未達理由.リードの支持を調べていない, l_理由);
            l_項目.Add(new 検査項目(
                "read_support", メッセージID.検査項目_リードの支持, l_支持,
                p_支持検査 is { } l_支持2
                    ? $"{l_支持2.A_支持のない位置数} <= {支持のない位置の許容数} ({l_支持2.A_区間.Count} stretch(es))"
                    : string.Empty));

            var (l_閉鎖, l_閉鎖の内訳) = Get_閉鎖の判定(p_閉鎖検証, l_理由);
            l_項目.Add(new 検査項目(
                "circular_closure", メッセージID.検査項目_環状閉鎖, l_閉鎖, l_閉鎖の内訳));

            var l_レベル = Get_品質保証レベル(
                l_取りこぼし, l_出しすぎ, l_深度, l_ギャップ, l_接合点, l_代替経路, l_閉鎖);

            return new 完全性判定結果(
                A_完全長か: l_レベル == 品質保証レベル.完全長,
                A_品質保証レベル: l_レベル,
                A_検査項目: l_項目,
                A_未達理由: [.. l_理由.Distinct()]);
        }

        /// <summary>
        /// 完全長の判定結果をログへ出力する
        /// </summary>
        /// <param name="p_判定">完全長の判定結果</param>
        public static void V_出力_判定結果(完全性判定結果 p_判定)
        {
            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.完全性の見出し);
            foreach (var l_項目 in p_判定.A_検査項目)
            {
                Logger.V_出力(
                    メッセージID.完全性の検査行,
                    Messages.Get_文言(l_項目.A_見出し),
                    Messages.Get_文言(Get_判定の見出し(l_項目.A_判定)),
                    l_項目.A_内訳);
            }
            Logger.V_出力(
                p_判定.A_完全長か ? メッセージID.完全長と判定 : メッセージID.完全長に届かず,
                (int)p_判定.A_品質保証レベル);
            if (p_判定.A_未達理由.Count > 0)
            {
                Logger.V_出力(
                    メッセージID.完全性の未達理由,
                    string.Join(", ", p_判定.A_未達理由.Select(Get_理由コード)));
            }
        }

        /// <summary>
        /// レポートに出す固定の理由コード
        /// </summary>
        /// <remarks>
        /// 訳さない
        /// </remarks>
        public static string Get_理由コード(未達理由 p_理由)
        {
            return p_理由 switch
            {
                未達理由.取りこぼしが多い => "kmer-missing",
                未達理由.出しすぎている => "kmer-excess",
                未達理由.自己検査を行えなかった => "self-check-unavailable",
                未達理由.リードに裏付けの無い箇所がある => "unsupported-sequence",
                未達理由.リードの支持を調べていない => "read-support-unchecked",
                未達理由.未解決のギャップが残る => "unresolved-gap",
                未達理由.決めきれない分岐が残る => "ambiguous-junction",
                未達理由.深度が不連続 => "depth-discontinuity",
                未達理由.深度を測っていない => "depth-unmeasured",
                未達理由.環状に閉じていない => "not-circular",
                未達理由.閉じ目がリードで裏付けられない => "closure-unsupported",
                _ => "closure-unverified",
            };
        }

        /// <summary>
        /// レポートに出す固定の判定名
        /// </summary>
        /// <remarks>
        /// 訳さない
        /// </remarks>
        public static string Get_判定コード(検査判定 p_判定)
        {
            return p_判定 switch
            {
                検査判定.合格 => "pass",
                検査判定.不合格 => "fail",
                _ => "unknown",
            };
        }

        /// <summary>
        /// 検査判定に対応する見出しの文言を返す
        /// </summary>
        /// <param name="p_判定">検査判定</param>
        /// <returns>見出しの文言</returns>
        private static メッセージID Get_判定の見出し(検査判定 p_判定)
        {
            return p_判定 switch
            {
                検査判定.合格 => メッセージID.検査判定_合格,
                検査判定.不合格 => メッセージID.検査判定_不合格,
                _ => メッセージID.検査判定_判定不能,
            };
        }

        /// <summary>
        /// 合否から検査判定を作り、不合格なら理由を書き留める
        /// </summary>
        /// <param name="p_合格か">合格したか</param>
        /// <param name="p_理由">不合格のときの理由</param>
        /// <param name="p_理由一覧">書き留める先</param>
        /// <returns>検査判定</returns>
        private static 検査判定 Get_判定(bool p_合格か, 未達理由 p_理由, List<未達理由> p_理由一覧)
        {
            if (p_合格か)
            {
                return 検査判定.合格;
            }
            p_理由一覧.Add(p_理由);
            return 検査判定.不合格;
        }

        /// <summary>
        /// 判定不能を返し、その理由を書き留める
        /// </summary>
        /// <param name="p_理由">判定できない理由</param>
        /// <param name="p_理由一覧">書き留める先</param>
        /// <returns>判定不能</returns>
        private static 検査判定 Get_判定不能(未達理由 p_理由, List<未達理由> p_理由一覧)
        {
            p_理由一覧.Add(p_理由);
            return 検査判定.判定不能;
        }

        private static (検査判定 A_判定, string A_内訳) Get_閉鎖の判定(
            IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証, List<未達理由> p_理由一覧)
        {
            if (p_閉鎖検証 is null)
            {
                return (Get_判定不能(未達理由.環状閉鎖を検証していない, p_理由一覧), string.Empty);
            }
            if (p_閉鎖検証.Count == 0)
            {
                p_理由一覧.Add(未達理由.環状に閉じていない);
                return (検査判定.不合格, "0");
            }
            var l_支持数 = p_閉鎖検証.Count(x => x.A_支持されたか);
            var l_内訳 = $"{l_支持数}/{p_閉鎖検証.Count}";
            if (l_支持数 == p_閉鎖検証.Count)
            {
                return (検査判定.合格, l_内訳);
            }
            p_理由一覧.Add(未達理由.閉じ目がリードで裏付けられない);
            return (検査判定.不合格, l_内訳);
        }

        /// <summary>
        /// 検査の合否を段階に畳む
        /// </summary>
        /// <remarks>
        /// 下の段が通っていない限り上の段は名乗れない
        /// </remarks>
        private static 品質保証レベル Get_品質保証レベル(
            検査判定 p_取りこぼし, 検査判定 p_出しすぎ, 検査判定 p_深度,
            検査判定 p_ギャップ, 検査判定 p_接合点, 検査判定 p_代替経路, 検査判定 p_閉鎖)
        {
            return p_取りこぼし != 検査判定.合格 || p_出しすぎ != 検査判定.合格
                ? 品質保証レベル.出力のみ
                : p_深度 != 検査判定.合格
                ? 品質保証レベル.グラフ整合
                : p_ギャップ != 検査判定.合格
                ? 品質保証レベル.マッピング整合
                : p_接合点 != 検査判定.合格
                ? 品質保証レベル.ペア整合
                : p_代替経路 != 検査判定.合格 || p_閉鎖 != 検査判定.合格
                ? 品質保証レベル.接合点が支持済み
                : 品質保証レベル.完全長;
        }

        /// <summary>
        /// 環状として出力された配列の本数
        /// </summary>
        public static int Get_環状本数(string p_FASTAパス)
        {
            var l_数 = 0;
            using var l_読み込み = new FastaReader(p_FASTAパス);
            while (l_読み込み.Get_続きがあるか())
            {
                if (l_読み込み.Get_次の配列().A_ID.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase))
                {
                    l_数++;
                }
            }
            return l_数;
        }

        /// <summary>
        /// 埋まらずに残った N の連続区間の数
        /// </summary>
        /// <remarks>
        /// 埋められなかったギャップは
        /// 「そこを繋いだ根拠が無い」ことをそのまま表している
        /// </remarks>
        public static int Get_未解決ギャップ数(string p_FASTAパス)
        {
            var l_数 = 0;
            using var l_読み込み = new FastaReader(p_FASTAパス);
            while (l_読み込み.Get_続きがあるか())
            {
                var l_配列 = l_読み込み.Get_次の配列().A_配列;
                var l_直前がNか = false;
                foreach (var l_文字 in l_配列)
                {
                    var l_Nか = l_文字 is 'N' or 'n';
                    if (l_Nか && !l_直前がNか)
                    {
                        l_数++;
                    }
                    l_直前がNか = l_Nか;
                }
            }
            return l_数;
        }
    }
}
