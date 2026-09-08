namespace Tsumiki.Model
{
    /// <summary>スキャフォールドのギャップ充填の集計。</summary>
    internal readonly record struct ギャップ充填統計(
        int A_総ギャップ数,
        int A_埋めたギャップ数,
        int A_埋めた塩基数,
        int A_一意に定まらなかった数,
        int A_到達できなかった数);
}
