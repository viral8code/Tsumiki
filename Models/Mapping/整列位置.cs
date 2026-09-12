namespace Tsumiki.Models.Mapping
{
    /// <summary>
    /// リードと参照で対応した 1 塩基の位置
    /// </summary>
    /// <param name="A_リード位置"></param>
    /// <param name="A_参照位置"></param>
    internal readonly record struct 整列位置(int A_リード位置, int A_参照位置);
}
