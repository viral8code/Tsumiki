namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// ロングリードで決めた、互いに一意な端の組 (頂点 u の右端から頂点 x の左端へ)
    /// </summary>
    /// <param name="A_u">左の頂点 (配列 i の順向き 2i、逆相補 2i+1)</param>
    /// <param name="A_x">右の頂点</param>
    /// <param name="A_票">首位の票 (跨いだリードの本数)</param>
    /// <param name="A_次点の票">u から見た次点の票</param>
    /// <param name="A_間の長さ">間の長さの中央値 (bp、負は重なり)</param>
    internal readonly record struct ロングリードの橋(int A_u, int A_x, int A_票, int A_次点の票, int A_間の長さ);
}
