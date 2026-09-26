namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// グラフ整理の 1 回の反復で、1 本の unitig について下した判定
    /// </summary>
    /// <param name="A_塩基列">unitig の塩基 ID 列</param>
    /// <param name="A_平均カバレッジ">unitig の平均カバレッジ</param>
    /// <param name="A_Is除去">unitig 全体を tip として除くか</param>
    /// <param name="A_先頭から">先頭から削る k-mer 数</param>
    /// <param name="A_末尾から">末尾から削る k-mer 数</param>
    /// <param name="A_合流側kmer群">削る k-mer のうち、合流点に接する端の k-mer (合流点へ向かう向き)</param>
    internal record 整理判定(byte[] A_塩基列, double A_平均カバレッジ, bool A_Is除去, int A_先頭から, int A_末尾から, byte[][] A_合流側kmer群);
}
