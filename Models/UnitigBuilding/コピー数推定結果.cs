namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// カバレッジからのコピー数推定の結果
    /// </summary>
    /// <param name="A_単一コピー基準値"></param>
    /// <param name="A_カバレッジ"></param>
    /// <param name="A_コピー数"></param>
    internal readonly record struct コピー数推定結果(double A_単一コピー基準値, コピー数基準の出所 A_基準の出所, IReadOnlyDictionary<int, double> A_カバレッジ, IReadOnlyDictionary<int, int> A_コピー数);

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
