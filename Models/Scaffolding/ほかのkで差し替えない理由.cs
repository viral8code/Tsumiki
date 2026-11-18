namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 最終配列とほかの k の contig で錨の間の長さが食い違う区間を、ほかの k の形に差し替えなかった理由 (OtherKRepeatReplacer)
    /// </summary>
    internal enum ほかのkで差し替えない理由
    {
        /// <summary>
        /// 両脇の片 (隣の N の列か配列の端まで) のどちらかが、リードで確かめる窓より短い
        /// </summary>
        片が短い,

        /// <summary>
        /// 最終配列の形がリードに続けて出てきて、窓で比べてもほかの k の形を選ばなかった
        /// </summary>
        骨格がリードにある,

        /// <summary>
        /// リードに続けて出てくるほかの k の形が、ただ 1 つに決まらない
        /// </summary>
        決まらない,

        /// <summary>
        /// ほかの k の形はただ 1 つに決まったが、錨の間の長さをリードで確かめられない
        /// </summary>
        錨を確かめられない,
    }
}
