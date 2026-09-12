namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// カバレッジからのコピー数推定の結果
    /// </summary>
    /// <param name="A_単一コピー基準値"></param>
    /// <param name="A_カバレッジ"></param>
    /// <param name="A_コピー数"></param>
    internal readonly record struct コピー数推定結果(double A_単一コピー基準値, IReadOnlyDictionary<int, double> A_カバレッジ, IReadOnlyDictionary<int, int> A_コピー数);
}
