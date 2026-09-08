namespace Tsumiki.Model
{
    /// <summary>FASTA の1エントリ(ヘッダ行のIDと塩基配列)。</summary>
    internal readonly struct 配列エントリ(string p_ID, string p_配列)
    {
        public readonly string A_ID = p_ID;
        public readonly string A_配列 = p_配列;
    }
}
