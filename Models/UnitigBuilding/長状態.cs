namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// 128 bit 2 語で表した k-mer の順鎖・逆鎖パック値
    /// </summary>
    /// <param name="A_順上"></param>
    /// <param name="A_順下"></param>
    /// <param name="A_逆上"></param>
    /// <param name="A_逆下"></param>
    internal readonly record struct 長状態(UInt128 A_順上, UInt128 A_順下, UInt128 A_逆上, UInt128 A_逆下);
}
