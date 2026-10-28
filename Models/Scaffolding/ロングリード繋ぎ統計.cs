namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// ロングリードの橋の間を埋めて繋いだ結果の統計
    /// </summary>
    /// <param name="A_重なりで畳んだ数">重なりで畳んだ橋の数</param>
    /// <param name="A_経路で埋めた数">グラフの経路 (一意なもの) で埋めた橋の数</param>
    /// <param name="A_ロングリードで選んだ数">グラフの道の分かれ道をロングリードの多数決で選んで埋めた橋の数</param>
    /// <param name="A_決まらなかった数">どちらでも決まらず繋がなかった橋の数</param>
    /// <param name="A_前の本数">繋ぐ前の配列の本数</param>
    /// <param name="A_後の本数">繋いだ後の配列の本数</param>
    /// <param name="A_直した塩基数">埋めた配列のうち、ロングリードと短いリードで直した塩基の数</param>
    internal readonly record struct ロングリード繋ぎ統計(int A_重なりで畳んだ数, int A_経路で埋めた数, int A_ロングリードで選んだ数, int A_決まらなかった数, int A_前の本数, int A_後の本数, int A_直した塩基数);
}
