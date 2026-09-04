namespace Tsumiki.Model
{
    /// <summary>1ファイル分のエラー訂正の集計。</summary>
    internal readonly record struct ファイル訂正統計(int A_総リード数, int A_訂正されたリード数, int A_総訂正塩基数);
}
