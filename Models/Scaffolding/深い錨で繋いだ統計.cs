namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 最終配列の N の列を、深い錨とほかの k の contig の配列で繋いだ結果の統計 (OtherKDeepGapBridger)<br/>
    /// 数はどれも N の列の数
    /// </summary>
    /// <param name="A_対象数">N の列の数</param>
    /// <param name="A_繋いだ数">繋いだ数 (畳んだ数を含む)</param>
    /// <param name="A_畳んだ数">繋いだうち、錨の間の配列が N の列のままの両脇の長さより短く、両脇の片の重なりを畳んだ数</param>
    /// <param name="A_決まらない数">ほかの k の錨の間の配列で、リードに続けて出てくるものがただ 1 つに決まらない数</param>
    /// <param name="A_リードで確かめられない数">左の錨から右の錨までがリード長より長いか、リードに無い数</param>
    internal readonly record struct 深い錨で繋いだ統計(int A_対象数, int A_繋いだ数, int A_畳んだ数, int A_決まらない数, int A_リードで確かめられない数);
}
