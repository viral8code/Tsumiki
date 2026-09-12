namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// scaffold 辺の候補
    /// </summary>
    /// <remarks>
    /// 観測本数だけでなく期待本数に対する比を持つ
    /// </remarks>
    /// <param name="A_行き先">接続先の頂点番号</param>
    /// <param name="A_支持数">距離が揃っているペアの本数</param>
    /// <param name="A_ギャップ長">支持ペアから推定した未知区間の長さ</param>
    /// <param name="A_期待に対する比">観測本数 / 期待本数</param>
    internal readonly record struct スキャフォールド候補(int A_行き先, ulong A_支持数, int A_ギャップ長, double A_期待に対する比);
}
