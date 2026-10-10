namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// scaffold 辺の候補
    /// </summary>
    /// <param name="A_行き先">接続先の頂点番号</param>
    /// <param name="A_支持数">距離が揃っているペアの本数</param>
    /// <param name="A_ギャップ長">支持ペアから推定した未知区間の長さ</param>
    /// <param name="A_期待に対する比">観測本数 / 期待本数</param>
    /// <param name="A_期待本数">理想本数モデルによる期待本数 (NaN は分からない、つまり理想本数モデルが無い)</param>
    internal readonly record struct Scaffold候補(int A_行き先, ulong A_支持数, int A_ギャップ長, double A_期待に対する比, double A_期待本数 = double.NaN);
}
