namespace Tsumiki.Models.Mapping
{
    /// <summary>
    /// リードと参照で対応した 1 塩基の位置
    /// </summary>
    internal readonly record struct 整列位置(int A_リード位置, int A_参照位置);
}
