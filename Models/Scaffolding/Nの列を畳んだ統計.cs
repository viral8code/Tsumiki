namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 最終配列の N の列のうち、両脇の片が完全に重なる所を畳んだ結果の統計 (OverlapGapFolder)<br/>
    /// 数はどれも N の列の数
    /// </summary>
    /// <param name="A_対象数">N の列の数</param>
    /// <param name="A_畳んだ数">重なりで畳んだ数</param>
    /// <param name="A_片が短い数">両脇の片 (隣の N の列か配列の端まで) のどちらかが、反復を調べる長さより短い数</param>
    /// <param name="A_重なりが無い数">左の片の末尾と右の片の先頭が完全に一致する重なりが無い数</param>
    /// <param name="A_反復の端の数">一致する重なりはあるが、どれも重なりを除いて両脇が反復の端である数</param>
    /// <param name="A_リードに無い数">反復の端でない重なりはあるが、畳んだ形がどれもリードに続けて出てこない数</param>
    /// <param name="A_決まらない数">畳んだ形がリードに続けて出てくる重なりの長さが 2 つ以上ある数</param>
    internal readonly record struct Nの列を畳んだ統計(int A_対象数, int A_畳んだ数, int A_片が短い数, int A_重なりが無い数, int A_反復の端の数, int A_リードに無い数, int A_決まらない数);
}
