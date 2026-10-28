namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// ロングリードで橋の組を決めた結果の統計
    /// </summary>
    /// <param name="A_リード数">当て方に使ったリードの本数</param>
    /// <param name="A_当たり数">採った当たりの数</param>
    /// <param name="A_錨数">錨の数</param>
    /// <param name="A_多コピーの錨数">多コピーとして橋に使わなかった錨の数</param>
    /// <param name="A_候補数">証拠のある橋の候補の数 (向きを区別しない組)</param>
    /// <param name="A_橋数">決めた橋の数</param>
    internal readonly record struct ロングリード橋渡し統計(int A_リード数, int A_当たり数, int A_錨数, int A_多コピーの錨数, int A_候補数, int A_橋数);
}
