namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// 初期値を変えて当てはめた 1 回ぶんの結果
    /// </summary>
    /// <param name="A_λ"></param>
    /// <param name="A_誤り平均"></param>
    /// <param name="A_誤り混合比"></param>
    /// <param name="A_コピー数減衰率"></param>
    /// <param name="A_対数尤度"></param>
    /// <param name="A_反復回数"></param>
    internal readonly record struct ヒストグラム試行結果(double A_λ, double A_誤り平均, double A_誤り混合比, double A_コピー数減衰率, double A_対数尤度, int A_反復回数);
}
