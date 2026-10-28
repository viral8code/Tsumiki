namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// ロングリードが 1 つの錨に当たった区間 (minimizer の共線的な鎖)
    /// </summary>
    /// <param name="A_錨">錨の番号</param>
    /// <param name="A_Is順向き">リードが錨の順向きに当たったか</param>
    /// <param name="A_リード始">リード上の始まり (含む)</param>
    /// <param name="A_リード終">リード上の終わり (含まない)</param>
    /// <param name="A_錨始">錨上の始まり (含む)</param>
    /// <param name="A_錨終">錨上の終わり (含まない)</param>
    /// <param name="A_種の数">鎖の minimizer の数</param>
    internal readonly record struct ロングリードの当たり(int A_錨, bool A_Is順向き, int A_リード始, int A_リード終, int A_錨始, int A_錨終, int A_種の数);
}
