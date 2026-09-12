namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 局所アセンブリで埋める 1 箇所のギャップ
    /// </summary>
    /// <param name="A_足場番号"></param>
    /// <param name="A_開始"></param>
    /// <param name="A_長さ"></param>
    /// <param name="A_左アンカー"></param>
    /// <param name="A_右アンカー"></param>
    internal readonly record struct 局所ギャップ(int A_足場番号, int A_開始, int A_長さ, string A_左アンカー, string A_右アンカー);
}
