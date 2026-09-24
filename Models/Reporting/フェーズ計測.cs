namespace Tsumiki.Models.Reporting
{
    /// <summary>
    /// 1 つの工程区間で計測した資源使用量
    /// </summary>
    /// <param name="A_工程"></param>
    /// <param name="A_経過秒"></param>
    /// <param name="A_CPU秒"></param>
    /// <param name="A_確保MB"></param>
    /// <param name="A_ワーキングセットMB"></param>
    /// <param name="A_ピークワーキングセットMB"></param>
    /// <param name="A_世代2回収回数"></param>
    internal readonly record struct フェーズ計測(string A_工程, double A_経過秒, double A_CPU秒, double A_確保MB, double A_ワーキングセットMB, double A_ピークワーキングセットMB, int A_世代2回収回数);
}
