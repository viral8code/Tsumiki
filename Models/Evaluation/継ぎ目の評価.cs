namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// 継ぎ目の候補に付けた、繋ぐ相手を誤っている確率
    /// </summary>
    /// <param name="A_候補"></param>
    /// <param name="A_組">較正に使った組</param>
    /// <param name="A_誤りの確率"></param>
    internal readonly record struct 継ぎ目の評価(継ぎ目候補 A_候補, 継ぎ目の組 A_組, double A_誤りの確率);
}
