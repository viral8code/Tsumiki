namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 継ぎ目の埋め方 (左の片の末尾と右の片の先頭を削り、間に合意の配列を入れる)
    /// </summary>
    /// <param name="A_左から削る長さ">左の片の末尾から削る長さ</param>
    /// <param name="A_埋める配列">間に入れる配列</param>
    /// <param name="A_右から削る長さ">右の片の先頭から削る長さ</param>
    internal readonly record struct 継ぎ目の埋め方(int A_左から削る長さ, string A_埋める配列, int A_右から削る長さ);
}
