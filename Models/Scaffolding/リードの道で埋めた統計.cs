namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 最終配列の N の列を、リードに 2 本以上ある k-mer の 1 本の道で埋めた結果の統計 (ReadGraphGapFiller)<br/>
    /// 数は k 長群を除いて、どれも N の列の数。埋めなかった理由は、最後に試した k での理由
    /// </summary>
    /// <param name="A_k長群">道を辿る k-mer の長さを試す順に並べたもの</param>
    /// <param name="A_対象数">N の列の数</param>
    /// <param name="A_埋めた数">埋めた数 (畳んだ数を含む)</param>
    /// <param name="A_畳んだ数">両脇の片が重なる形 (埋めた間が負) で畳んだ数</param>
    /// <param name="A_小さいkで埋めた数">最初の k では埋まらず、小さい k で埋めた数</param>
    /// <param name="A_窓で選んだ数">道が 2 本以上あり、短い順に並べた道をリードの窓で比べて選んだ道で埋めた数</param>
    /// <param name="A_片が短い数">両脇の片のどちらかが、反復を調べる長さ (と k 長 + 端から削る長さの上限) より短い数</param>
    /// <param name="A_端がリードに無い数">片の端から削る長さの上限まで内へずらしても、端の k-mer がリードに 2 本以上ない数</param>
    /// <param name="A_道が無い数">始まりの k-mer から終わりの k-mer への道が無い数</param>
    /// <param name="A_道が決まらない数">道が 2 本以上あるか、進めた回数が上限を超えて、窓で比べても選べなかった数</param>
    /// <param name="A_畳む形を確かめられない数">畳む形で、深い錨から深い錨までがリード長より長いか、リードに無いか、錨が無い数</param>
    /// <param name="A_錨から錨までを確かめられない数">埋める形で、道の始まりの錨から終わりの錨までがリード長より長いか、リードに無い数</param>
    internal readonly record struct リードの道で埋めた統計(int[] A_k長群, int A_対象数, int A_埋めた数, int A_畳んだ数, int A_小さいkで埋めた数, int A_窓で選んだ数, int A_片が短い数, int A_端がリードに無い数, int A_道が無い数, int A_道が決まらない数, int A_畳む形を確かめられない数, int A_錨から錨までを確かめられない数);
}
