namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// FASTA の 1 エントリ (ヘッダ行の ID と塩基配列)
    /// </summary>
    /// <param name="p_ID"></param>
    /// <param name="p_配列"></param>
    internal readonly struct 配列エントリ(string p_ID, string p_配列)
    {
        #region 内部変数

        /// <summary>
        /// 配列 ID
        /// </summary>
        public readonly string A_ID = p_ID;

        /// <summary>
        /// 配列
        /// </summary>
        public readonly string A_配列 = p_配列;

        #endregion
    }
}
