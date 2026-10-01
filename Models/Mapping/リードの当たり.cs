namespace Tsumiki.Models.Mapping
{
    /// <summary>
    /// 1 本のリードの配置を、配列上の整列範囲と両端の切れ端にまとめたもの
    /// </summary>
    /// <param name="A_配列番号"></param>
    /// <param name="A_開始">整列の始まり (0 始まり)</param>
    /// <param name="A_終了">整列の終わり (含まない)</param>
    /// <param name="A_左の切れ端">配列の左側で整列から外れたリードの塩基数</param>
    /// <param name="A_右の切れ端">配列の右側で整列から外れたリードの塩基数</param>
    /// <param name="A_Is逆鎖"></param>
    /// <param name="A_信頼度"></param>
    internal readonly record struct リードの当たり(int A_配列番号, int A_開始, int A_終了, int A_左の切れ端, int A_右の切れ端, bool A_Is逆鎖, int A_信頼度);
}
