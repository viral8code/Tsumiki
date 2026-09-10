namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// 次の k へ引き継ぐ、前段のアセンブリの 1 配列
    /// </summary>
    /// <remarks>
    /// カバレッジを配列と一緒に持つのが要点<br/>
    /// 引き継いだ k-mer に名目値を与えると、
    /// コピー数推定・低カバレッジ端のトリミング・自己検査がまとめて壊れる<br/>
    /// A_カバレッジ[i] は、この配列の位置 i から始まる A_k長 塩基の k-mer の出現回数
    /// </remarks>
    internal record 引き継ぎ配列(string A_配列, int[] A_カバレッジ, int A_k長);
}
