namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// ロングリードの橋の間をどう埋めたか
    /// </summary>
    /// <param name="A_重なり">重なりで畳んだときの重なりの長さ、畳まなければ null</param>
    /// <param name="A_埋める配列">グラフの経路で埋めたときの間の配列、埋めなければ null</param>
    /// <param name="A_理由">繋がなかったときの理由 (no_overlap / ambiguous / unreachable / aborted / no_graph)、繋いだら null</param>
    /// <param name="A_Isロングリードで選んだ">グラフの道の分かれ道をロングリードの多数決で選んで埋めたか</param>
    internal readonly record struct ロングリードの繋ぎ方(int? A_重なり, string? A_埋める配列, string? A_理由, bool A_Isロングリードで選んだ = false);
}
