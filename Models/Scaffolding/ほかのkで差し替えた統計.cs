namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 最終配列と錨の間の長さが食い違う区間を、ほかの k の contig の形に差し替えた結果の統計 (OtherKRepeatReplacer)<br/>
    /// 数はどれも区間の数
    /// </summary>
    /// <param name="A_対象数">ほかの k の contig と錨の間の長さが食い違う区間の数</param>
    /// <param name="A_差し替えた数">ほかの k の形に差し替えた数</param>
    /// <param name="A_片が短い数">両脇の片 (隣の N の列か配列の端まで) のどちらかが、リードで確かめる窓より短い数</param>
    /// <param name="A_骨格がリードにある数">最終配列の形がリードに続けて出てくるので、そのまま残した数</param>
    /// <param name="A_決まらない数">リードに続けて出てくるほかの k の形が、ただ 1 つに決まらない数</param>
    /// <param name="A_錨を確かめられない数">ほかの k の形はただ 1 つに決まったが、左の錨から、最終配列とほかの k の区間の短いほうに錨の長さを足した長さの配列が、リード長より長いかリードに無い数</param>
    internal readonly record struct ほかのkで差し替えた統計(int A_対象数, int A_差し替えた数, int A_片が短い数, int A_骨格がリードにある数, int A_決まらない数, int A_錨を確かめられない数);
}
