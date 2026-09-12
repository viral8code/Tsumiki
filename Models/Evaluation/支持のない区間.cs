namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// アセンブリ上の、どのリードにも現れない r-mer が連なっている区間
    /// </summary>
    /// <remarks>
    /// 位置は 1 始まり・両端を含む
    /// </remarks>
    /// <param name="A_配列ID"></param>
    /// <param name="A_開始"></param>
    /// <param name="A_終了"></param>
    internal readonly record struct 支持のない区間(string A_配列ID, int A_開始, int A_終了)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 区間が覆う塩基の数
        /// </summary>
        public int A_長さ => this.A_終了 - this.A_開始 + 1;

        #endregion
    }
}
