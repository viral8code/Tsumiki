namespace Tsumiki.Model.Foundation
{
    /// <summary>
    /// FASTA の 1 エントリ (ヘッダ行の ID と塩基配列)
    /// </summary>
    internal readonly struct 配列エントリ(string p_ID, string p_配列)
    {
        public readonly string A_ID = p_ID;
        /// <summary>
        /// 配列
        /// </summary>
        public readonly string A_配列 = p_配列;
    }
}
