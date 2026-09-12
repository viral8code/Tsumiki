namespace Tsumiki.Models.Correction
{
    /// <summary>
    /// 1 リードのエラー訂正結果 (訂正後の塩基列と、訂正した塩基数)
    /// </summary>
    /// <param name="A_塩基列"></param>
    /// <param name="A_訂正数"></param>
    internal readonly record struct 訂正結果(byte[] A_塩基列, int A_訂正数);
}
