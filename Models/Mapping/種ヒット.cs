namespace Tsumiki.Models.Mapping
{
    /// <summary>
    /// 索引中の種の位置
    /// </summary>
    internal readonly record struct 種ヒット(int A_配列番号, int A_参照位置, bool A_逆鎖か);
}
