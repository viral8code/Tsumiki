namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 最終配列の N の列をリードの道で埋めなかった理由 (ReadGraphGapFiller)
    /// </summary>
    internal enum リードの道で埋めない理由
    {
        /// <summary>
        /// 両脇の片のどちらかが、反復を調べる長さ (と k 長 + 端から削る長さの上限) より短い
        /// </summary>
        片が短い,

        /// <summary>
        /// 片の端から削る長さの上限まで内へずらしても、端の k-mer がリードに 2 本以上ない
        /// </summary>
        端がリードに無い,

        /// <summary>
        /// 始まりの k-mer から終わりの k-mer への道が無い
        /// </summary>
        道が無い,

        /// <summary>
        /// 道が 2 本以上あるか、進めた回数が上限を超えた
        /// </summary>
        道が決まらない,

        /// <summary>
        /// 道は 1 本で両脇の片が重なる形 (埋めた間が負) だが、畳んだ形の深い錨から深い錨までがリード長より長いか、リードに無いか、錨が無い
        /// </summary>
        畳む形を確かめられない,

        /// <summary>
        /// 道の始まりの錨から終わりの錨までがリード長より長いか、リードに無い
        /// </summary>
        錨から錨までを確かめられない,
    }
}
