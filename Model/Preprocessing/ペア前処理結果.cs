namespace Tsumiki.Model.Preprocessing
{
    /// <summary>
    /// 1 ペア分の前処理結果
    /// </summary>
    /// <remarks>
    /// 副作用のない純粋関数の戻り値
    /// </remarks>
    internal readonly record struct ペア前処理結果(
        string A_配列1, string A_クオリティ1, string A_配列2, string A_クオリティ2,
        bool A_アダプタを検出したか, int A_訂正塩基数);
}
