namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// リードの合意で埋める対象の継ぎ目
    /// </summary>
    /// <param name="A_配列番号">継ぎ目がある配列の番号</param>
    /// <param name="A_開始">継ぎ目 (ギャップの文字) の始まりの位置</param>
    /// <param name="A_終わり">継ぎ目の直後の位置</param>
    /// <param name="A_左の片">継ぎ目の左の片の末尾</param>
    /// <param name="A_右の片">継ぎ目の右の片の先頭</param>
    internal sealed record 合意の継ぎ目(int A_配列番号, int A_開始, int A_終わり, string A_左の片, string A_右の片);
}
