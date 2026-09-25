namespace Tsumiki.Models.Reporting
{
    /// <summary>
    /// 決めきれなかった 1 箇所の記録
    /// </summary>
    /// <param name="A_k長"></param>
    /// <param name="A_種別"></param>
    /// <param name="A_場所"></param>
    /// <param name="A_首位の支持"></param>
    /// <param name="A_次点の支持"></param>
    /// <param name="A_首位の生支持数"></param>
    /// <param name="A_確信度"></param>
    /// <param name="A_安定ID"></param>
    internal readonly record struct 曖昧箇所(int A_k長, 曖昧箇所の種別 A_種別, string A_場所, double A_首位の支持, double A_次点の支持, long A_首位の生支持数, double A_確信度, string A_安定ID = "")
    {
        #region カスタムプロパティ

        /// <summary>
        /// 首位と次点の差
        /// </summary>
        public double A_余裕 => this.A_首位の支持 - this.A_次点の支持;

        #endregion
    }
}
