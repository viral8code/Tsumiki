namespace Tsumiki.Models.Correction
{
    /// <summary>
    /// 1 ファイル分のエラー訂正の集計
    /// </summary>
    /// <param name="A_総リード数"></param>
    /// <param name="A_訂正されたリード数"></param>
    /// <param name="A_総訂正塩基数"></param>
    internal readonly record struct ファイル訂正統計(int A_総リード数, int A_訂正されたリード数, int A_総訂正塩基数);
}
