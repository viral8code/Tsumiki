namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// 単一コピー深度基準の由来
    /// </summary>
    internal enum コピー数基準の出所
    {
        /// <summary>
        /// k-mer スペクトル混合モデル
        /// </summary>
        Spectrum,

        /// <summary>
        /// unitig カバレッジの長さ加重中央値
        /// </summary>
        Weighted,
    }
}
