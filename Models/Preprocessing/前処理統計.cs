namespace Tsumiki.Models.Preprocessing
{
    /// <summary>
    /// ペア前処理 (オーバーラップ解析による アダプタ除去・相互訂正) 1 ファイル分の集計
    /// </summary>
    /// <param name="A_総ペア数"></param>
    /// <param name="A_アダプタ検出ペア数"></param>
    /// <param name="A_訂正塩基数"></param>
    /// <param name="A_総塩基数"></param>
    /// <param name="A_品質トリム塩基数"></param>
    /// <param name="A_品質トリム閾値"></param>
    /// <param name="A_polyGを削ったリード数">3' 末端の poly-G を削ったリードの数</param>
    /// <param name="A_polyGで削った塩基数">3' 末端の poly-G で削った塩基の数</param>
    internal readonly record struct 前処理統計(int A_総ペア数, int A_アダプタ検出ペア数, int A_訂正塩基数, long A_総塩基数 = 0L, long A_品質トリム塩基数 = 0L, int A_品質トリム閾値 = 0, int A_polyGを削ったリード数 = 0, long A_polyGで削った塩基数 = 0L);
}
