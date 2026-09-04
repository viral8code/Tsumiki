namespace Tsumiki.Model
{
    /// <summary>
    /// カバレッジからのコピー数推定の結果。
    /// </summary>
    internal readonly record struct コピー数推定結果(
        double A_単一コピー基準値,
        IReadOnlyDictionary<int, double> A_カバレッジ,
        IReadOnlyDictionary<int, int> A_コピー数);
}
