namespace Tsumiki.Models.Reporting
{
    /// <summary>
    /// 完全長かどうかの判定と、そう判定した根拠一式
    /// </summary>
    /// <param name="A_Is完全長"></param>
    /// <param name="A_品質保証レベル"></param>
    /// <param name="A_検査項目"></param>
    /// <param name="A_未達理由"></param>
    internal sealed record 完全性判定結果(bool A_Is完全長, 品質保証レベル A_品質保証レベル, IReadOnlyList<検査項目> A_検査項目, IReadOnlyList<未達理由> A_未達理由);
}
