namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// 単一コピー集団のカバレッジが、ポアソン分布が想定する以上にばらついていないかの診断
    /// </summary>
    /// <param name="A_平均"></param>
    /// <param name="A_分散"></param>
    /// <param name="A_分散指数"></param>
    /// <param name="A_Is過分散"></param>
    /// <remarks>
    /// 分散指数 (分散 ÷ 平均) はポアソン分布であれば 1 に近づく<br/>
    /// これを大きく超える場合、真のコピー数比較には単純な比 (暗黙のポアソン仮定) より、
    /// 負の二項分布のような過分散を許すモデルの方が適切な可能性がある、という診断に留める<br/>
    /// この診断だけでコピー数推定のロジックを変えることはしない (低信頼を理由に anchor 除外や枝刈りをしない、という方針を守るため)
    /// </remarks>
    internal readonly record struct カバレッジ分散診断(double A_平均, double A_分散, double A_分散指数, bool A_Is過分散);
}
