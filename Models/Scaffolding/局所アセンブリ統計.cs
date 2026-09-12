namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// スキャフォールドの局所アセンブリ (-la) によるギャップ充填の集計
    /// </summary>
    /// <param name="A_対象ギャップ数"></param>
    /// <param name="A_埋めたギャップ数"></param>
    /// <param name="A_埋めた塩基数"></param>
    /// <param name="A_局所リードが集まらなかった数"></param>
    /// <param name="A_一意に定まらなかった数"></param>
    /// <param name="A_到達できなかった数"></param>
    internal readonly record struct 局所アセンブリ統計(int A_対象ギャップ数, int A_埋めたギャップ数, int A_埋めた塩基数, int A_局所リードが集まらなかった数, int A_一意に定まらなかった数, int A_到達できなかった数);
}
