using Tsumiki.Models.Foundation;

namespace Tsumiki.Models.Reporting
{
    /// <summary>
    /// 品質保証の段階
    /// </summary>
    /// <remarks>
    /// 完全長を名乗れるのは最上位だけとし、それ以外は「どこまでは言えるのか」を段階で示す
    /// </remarks>
    internal enum 品質保証レベル
    {
        /// <summary>
        /// 配列を出力できた
        /// </summary>
        出力のみ = 0,

        /// <summary>
        /// k-mer の取りこぼしと出しすぎが許容範囲に収まる
        /// </summary>
        グラフ整合 = 1,

        /// <summary>
        /// 元リードを貼り直しても深度が途切れない
        /// </summary>
        マッピング整合 = 2,

        /// <summary>
        /// 未解決のギャップが残っていない
        /// </summary>
        ペア整合 = 3,

        /// <summary>
        /// 決めきれずに打ち切った分岐が残っていない
        /// </summary>
        接合点が支持済み = 4,

        /// <summary>
        /// 閉じ目まで元リードで裏付けられている
        /// </summary>
        完全長 = 5,
    }

    /// <summary>
    /// 1 つの検査項目の結果
    /// </summary>
    /// <remarks>
    /// 材料が無い場合は不合格と区別する
    /// </remarks>
    internal enum 検査判定
    {
        /// <summary>
        /// 検査項目を満たした
        /// </summary>
        合格,

        /// <summary>
        /// 検査項目を満たさなかった
        /// </summary>
        不合格,

        /// <summary>
        /// 判定に必要な材料が揃っておらず、合否を言えない
        /// </summary>
        判定不能,
    }

    /// <summary>
    /// 完全長に届かなかった理由
    /// </summary>
    /// <remarks>
    /// レポートの reason_codes になる
    /// </remarks>
    internal enum 未達理由
    {
        /// <summary>
        /// グラフから期待される k-mer の取りこぼしが許容範囲を超えている
        /// </summary>
        取りこぼしが多い,

        /// <summary>
        /// グラフに含まれるべきでない k-mer の出現が許容範囲を超えている
        /// </summary>
        出しすぎている,

        /// <summary>
        /// グラフ整合の自己検査自体を実行できなかった
        /// </summary>
        自己検査を行えなかった,

        /// <summary>
        /// scaffold のギャップが埋まらずに残っている
        /// </summary>
        未解決のギャップが残る,

        /// <summary>
        /// 伸長中に一意へ絞り込めなかった分岐が残っている
        /// </summary>
        決めきれない分岐が残る,

        /// <summary>
        /// リードの再マッピング深度が接合点の前後で途切れている
        /// </summary>
        深度が不連続,

        /// <summary>
        /// 深度の連続性チェック自体を行っていない
        /// </summary>
        深度を測っていない,

        /// <summary>
        /// 配列が環状ゲノムとして閉じていない
        /// </summary>
        環状に閉じていない,

        /// <summary>
        /// 環状の閉じ目がリードで裏付けられていない
        /// </summary>
        閉じ目がリードで裏付けられない,

        /// <summary>
        /// 環状閉鎖の検証自体を行っていない
        /// </summary>
        環状閉鎖を検証していない,

        /// <summary>
        /// 元リードで裏付けられない箇所が配列中に残っている
        /// </summary>
        リードに裏付けの無い箇所がある,

        /// <summary>
        /// リードによる支持の有無を調べていない
        /// </summary>
        リードの支持を調べていない,
    }

    /// <summary>
    /// 検査 1 項目
    /// </summary>
    /// <remarks>
    /// A_キー はレポートに出す固定の英語キー、A_見出し はログに出す訳語
    /// </remarks>
    /// <param name="A_キー"></param>
    /// <param name="A_見出し"></param>
    /// <param name="A_判定"></param>
    /// <param name="A_内訳"></param>
    internal readonly record struct 検査項目(string A_キー, メッセージID A_見出し, 検査判定 A_判定, string A_内訳);

    /// <summary>
    /// 完全長かどうかの判定と、そう判定した根拠一式
    /// </summary>
    /// <remarks>
    /// 完全長は「長い配列が出た」ことではなく、必要な検査を全て通ったことを指す<br/>
    /// 情報が足りない箇所を推測で埋めて完全長を名乗らせないための型
    /// </remarks>
    /// <param name="A_Is完全長"></param>
    /// <param name="A_品質保証レベル"></param>
    /// <param name="A_検査項目"></param>
    /// <param name="A_未達理由"></param>
    internal sealed record 完全性判定結果(bool A_Is完全長, 品質保証レベル A_品質保証レベル, IReadOnlyList<検査項目> A_検査項目, IReadOnlyList<未達理由> A_未達理由);
}
