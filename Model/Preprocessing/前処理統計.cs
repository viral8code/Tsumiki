namespace Tsumiki.Model.Preprocessing
{
    /// <summary>
    /// ペア前処理 (オーバーラップ解析による アダプタ除去・相互訂正) 1 ファイル分の集計
    /// </summary>
    internal readonly record struct 前処理統計(int A_総ペア数, int A_アダプタ検出ペア数, int A_訂正塩基数);
}
