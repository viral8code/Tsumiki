namespace Tsumiki.Model.Evaluation
{
    /// <summary>
    /// 環状に閉じた 1 本について、閉じ目が元リードから裏付けられるかを調べた結果
    /// </summary>
    /// <remarks>
    /// 
    /// グラフ上で閉じたことと、その閉じ目を実際にリードが読んでいることは
    /// 別の主張であり、後者が無い環状は線状の断片と区別できない
    /// </remarks>
    internal readonly record struct 環状閉鎖検証結果(string A_配列ID, int A_長さ, int A_跨いだリード数, int A_必要本数)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 閉じ目を直接読んだリードが必要本数に達しているか
        /// </summary>
        public bool A_支持されたか => this.A_跨いだリード数 >= this.A_必要本数;

        #endregion
    }
}
