namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 1 つの継ぎ目の判定に使ったリードの数 (詳細ログに出す)
    /// </summary>
    /// <param name="A_跨ぐ数">両側の片に置けたリードの数</param>
    /// <param name="A_間">票で決まった間の長さ、決まらなければ null</param>
    /// <param name="A_食い違う数">片側だけに置けたリードのうち、繋いだ形のもう一方の片と確かな塩基で食い違うものの数</param>
    /// <param name="A_合う数">片側だけに置けたリードのうち、繋いだ形のもう一方の片と合うものの数</param>
    internal readonly record struct 継ぎ目合意の内訳(int A_跨ぐ数, int? A_間, int A_食い違う数, int A_合う数);
}
