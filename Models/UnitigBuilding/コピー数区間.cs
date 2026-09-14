namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// 点推定ではなく、観測されたカバレッジの分散を踏まえた妥当なコピー数の範囲
    /// </summary>
    /// <param name="A_下限"></param>
    /// <param name="A_上限"></param>
    /// <remarks>
    /// 点推定 (コピー数推定結果.A_コピー数) は必ずこの区間に含まれる<br/>
    /// 反復配列を安全に何回まで通ってよいかの予算のように、誤推定で真の経路を消してはいけない場面で
    /// 点推定の代わりにこちらを使う<br/>
    /// 逆に、分岐選択のように保守的であるべき場面 (足場の単一コピー判定等) では、これまでどおり点推定を使う
    /// </remarks>
    internal readonly record struct コピー数区間(int A_下限, int A_上限);
}
