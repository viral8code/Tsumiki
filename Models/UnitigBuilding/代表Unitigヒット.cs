namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// ContigMaker でリード 1 本を unitig の索引で走査し、代表として選んだ unitig
    /// </summary>
    /// <param name="p_unitigID"></param>
    /// <param name="p_最終一致終端位置"></param>
    /// <param name="p_unitig長"></param>
    internal readonly struct 代表Unitigヒット(int p_unitigID, int p_最終一致終端位置, int p_unitig長)
    {
        #region 定数

        /// <summary>
        /// ヒットなし
        /// </summary>
        public static readonly 代表Unitigヒット C_ヒットなし = new(0, 0, 0);

        #endregion

        #region 内部変数

        /// <summary>
        /// マップ先 unitig ID
        /// </summary>
        public readonly int A_unitigID = p_unitigID;

        /// <summary>
        /// unitig の向きに揃えた最後の一致 k-mer の終端位置 (末尾の添字 + 1)
        /// </summary>
        public readonly int A_最終一致終端位置 = p_最終一致終端位置;

        /// <summary>
        /// マップ先 unitig の全長 (向きに依存せず同じ)
        /// </summary>
        public readonly int A_unitig長 = p_unitig長;

        #endregion

        #region プロパティ

        /// <summary>
        /// unitig の末尾から、リードが最後にヒットした位置までの残り塩基数
        /// </summary>
        public int A_末尾までの残り長 => Math.Max(0, this.A_unitig長 - this.A_最終一致終端位置);

        #endregion
    }
}
