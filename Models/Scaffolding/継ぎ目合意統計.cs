namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 継ぎ目をリードの合意で埋める処理 (JunctionConsensusFiller) の結果の統計
    /// </summary>
    /// <param name="A_対象数">両側の片が十分に長い継ぎ目の数</param>
    /// <param name="A_埋めた数">リードの合意で埋めた継ぎ目の数</param>
    /// <param name="A_跨ぐリード無し数">両側の片に置けるリードが無かった継ぎ目の数</param>
    /// <param name="A_票が足りない数">最多の間の長さを示すリードが要る数に満たなかった継ぎ目の数</param>
    /// <param name="A_票が割れた数">最多の間の長さが次点の優勢比倍に届かなかった継ぎ目の数</param>
    /// <param name="A_合意が弱い数">列ごとの多数決で決まらない列があった継ぎ目の数 (片が短く余白を取れなかったものを含む)</param>
    /// <param name="A_縦列反復数">繋いだ配列が縦に並んだ反復の中なので埋めなかった継ぎ目の数</param>
    /// <param name="A_食い違い数">片側だけに置けるリードの多くが繋いだ形と食い違うので埋めなかった継ぎ目の数</param>
    /// <param name="A_片が収まる数">重なりが短いほうの片の長さ以上なので埋めなかった継ぎ目の数</param>
    internal readonly record struct 継ぎ目合意統計(int A_対象数, int A_埋めた数, int A_跨ぐリード無し数, int A_票が足りない数, int A_票が割れた数, int A_合意が弱い数, int A_縦列反復数, int A_食い違い数, int A_片が収まる数);
}
