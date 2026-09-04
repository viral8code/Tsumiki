namespace Tsumiki.Model
{
    /// <summary>1リードのエラー訂正結果(訂正後の塩基列と、訂正した塩基数)。</summary>
    internal readonly record struct 訂正結果(byte[] A_塩基列, int A_訂正数);
}
