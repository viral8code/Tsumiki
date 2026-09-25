namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// カバレッジからのコピー数推定の結果
    /// </summary>
    /// <param name="A_単一コピー基準値"></param>
    /// <param name="A_カバレッジ"></param>
    /// <param name="A_コピー数"></param>
    /// <param name="A_分散診断">単一コピー集団の過分散診断、標本が少なすぎる等で求められない場合は null</param>
    /// <param name="A_コピー数区間">
    /// 観測された分散を踏まえた妥当なコピー数の範囲 (unitig ID -> 区間)<br/>
    /// A_分散診断 が求められない場合は null (根拠のない区間を作らない)
    /// </param>
    internal readonly record struct コピー数推定結果(double A_単一コピー基準値, コピー数基準の出所 A_基準の出所, IReadOnlyDictionary<int, double> A_カバレッジ, IReadOnlyDictionary<int, int> A_コピー数, カバレッジ分散診断? A_分散診断 = null, IReadOnlyDictionary<int, コピー数区間>? A_コピー数区間 = null);
}
