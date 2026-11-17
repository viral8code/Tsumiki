namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 道が 2 本以上の N の列で、短い順に並べた道をリードの窓で比べた結果 (ReadGraphGapFiller)
    /// </summary>
    /// <param name="A_道の数">並べた道の数</param>
    /// <param name="A_Is打ち切り">並べる途中で、進めた回数が上限を超えたか</param>
    /// <param name="A_畳む形の数">並べた道のうち、畳む形 (埋めた間が負) の数</param>
    /// <param name="A_条件を満たす数">選ぶ条件を満たした道の数</param>
    /// <param name="A_選んだ順">選んだ道の短い順の番号 (1 始まり。選ばなければ 0)</param>
    /// <param name="A_窓長">窓の長さ</param>
    /// <param name="A_相手の窓の最少">選んだ道に対して、ほかの道の形にだけある窓の数の最少 (選ばなければ 0)</param>
    internal readonly record struct 窓で道を選んだ記録(int A_道の数, bool A_Is打ち切り, int A_畳む形の数, int A_条件を満たす数, int A_選んだ順, int A_窓長, int A_相手の窓の最少);
}
