namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 最終配列の N の列を、ほかの k の contig の配列で埋めた結果の統計 (OtherKGapFiller)<br/>
    /// 数はどれも N の列の数
    /// </summary>
    /// <param name="A_対象数">N の列の数</param>
    /// <param name="A_埋めた数">ほかの k の配列で埋めた数</param>
    /// <param name="A_錨が取れない数">両脇の錨が取れないか、アセンブリの中で 1 回だけでない数</param>
    /// <param name="A_反復の端の数">両脇の片 (隣の N の列か配列の端まで) のどちらかが反復を調べる長さより短いか、両脇が反復の端である数</param>
    /// <param name="A_ほかのkに無い数">両脇の錨を同じ向きで順に含むほかの k の contig が無い数</param>
    /// <param name="A_リードに無い数">埋める配列の候補はあったが、リードに続けて出てくる候補がただ 1 つに決まらなかった数</param>
    internal readonly record struct ほかのkで埋めた統計(int A_対象数, int A_埋めた数, int A_錨が取れない数, int A_反復の端の数, int A_ほかのkに無い数, int A_リードに無い数);
}
