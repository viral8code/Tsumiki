namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// 環状に閉じた 1 本について、閉じ目が元リードから裏付けられるかを調べた結果
    /// </summary>
    /// <param name="A_配列ID"></param>
    /// <param name="A_長さ"></param>
    /// <param name="A_跨いだリード数"></param>
    /// <param name="A_必要本数"></param>
    internal readonly record struct 環状閉鎖検証結果(string A_配列ID, int A_長さ, int A_跨いだリード数, int A_必要本数)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 閉じ目を直接読んだリードが必要本数に達しているか
        /// </summary>
        public bool A_Has支持 => this.A_跨いだリード数 >= this.A_必要本数;

        #endregion
    }
}
