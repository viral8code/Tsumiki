namespace Tsumiki.Models.Preprocessing
{
    /// <summary>
    /// 1 ペア分の前処理結果
    /// </summary>
    /// <param name="A_順リード配列"></param>
    /// <param name="A_順リード品質"></param>
    /// <param name="A_逆リード配列"></param>
    /// <param name="A_逆リード品質"></param>
    /// <param name="A_Hasアダプタ検出"></param>
    /// <param name="A_訂正塩基数"></param>
    /// <param name="A_品質トリム塩基数">3' 末端の品質トリムで切った塩基数 (両側の合計)</param>
    internal readonly record struct ペア前処理結果(string A_順リード配列, string A_順リード品質, string A_逆リード配列, string A_逆リード品質, bool A_Hasアダプタ検出, int A_訂正塩基数, int A_品質トリム塩基数 = 0);
}
