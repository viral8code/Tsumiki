namespace Tsumiki.Models.Mapping
{
    /// <summary>
    /// 索引中の種の位置
    /// </summary>
    /// <param name="A_配列番号"></param>
    /// <param name="A_参照位置"></param>
    /// <param name="A_Is逆鎖"></param>
    internal readonly record struct 種ヒット(int A_配列番号, int A_参照位置, bool A_Is逆鎖);
}
