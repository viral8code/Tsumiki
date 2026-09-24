namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// 点推定ではなく、観測されたカバレッジの分散を踏まえた妥当なコピー数の範囲
    /// </summary>
    /// <param name="A_下限"></param>
    /// <param name="A_上限"></param>
    internal readonly record struct コピー数区間(int A_下限, int A_上限);
}
