namespace Tsumiki.Model.Evaluation
{
    /// <summary>
    /// アセンブリ上の、どのリードにも現れない r-mer が連なっている区間。
    /// 位置は 1 始まり・両端を含む。
    /// </summary>
    internal readonly record struct 支持のない区間(string A_配列ID, int A_開始, int A_終了)
    {
        public int A_長さ => this.A_終了 - this.A_開始 + 1;
    }
}
