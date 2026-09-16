namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// 次の k へ引き継ぐ、前段のアセンブリの 1 配列
    /// </summary>
    /// <remarks>
    /// カバレッジを配列と一緒に持つのが要点<br/>
    /// 引き継いだ k-mer に名目値を与えると、コピー数推定・低カバレッジ端のトリミング・自己検査がまとめて壊れる<br/>
    /// A_カバレッジ[i] は、この配列の位置 i から始まる A_k長 塩基の k-mer の出現回数
    /// </remarks>
    /// <param name="A_配列"></param>
    /// <param name="A_カバレッジ"></param>
    /// <param name="A_k長"></param>
    /// <param name="A_Is確定経路">
    /// この段の scaffold/contig 全体から取った、既に確定した経路由来なら true<br/>
    /// バブル敗者や合成リードのような、経路として確定していない配列は false のままにする
    /// </param>
    /// <param name="A_分岐の継ぎ目位置">
    /// 前段 k で分岐のある継ぎ目を通った辺 ((k+1)-mer) が現れる開始位置、無ければ null<br/>
    /// この辺を丸ごと含む次の k の k-mer は足さない
    /// </param>
    internal record 引き継ぎ配列(string A_配列, int[] A_カバレッジ, int A_k長, bool A_Is確定経路 = false, IReadOnlyList<int>? A_分岐の継ぎ目位置 = null);
}
