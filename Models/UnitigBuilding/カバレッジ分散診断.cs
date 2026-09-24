namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// 単一コピー集団のカバレッジが、ポアソン分布が想定する以上にばらついていないかの診断
    /// </summary>
    /// <param name="A_平均"></param>
    /// <param name="A_分散"></param>
    /// <param name="A_分散指数"></param>
    /// <param name="A_Is過分散"></param>
    internal readonly record struct カバレッジ分散診断(double A_平均, double A_分散, double A_分散指数, bool A_Is過分散);
}
