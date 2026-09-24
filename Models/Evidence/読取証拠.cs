namespace Tsumiki.Models.Evidence
{
    /// <summary>
    /// 局所アセンブリへ回収した 1 本のリードと、その由来
    /// </summary>
    /// <param name="A_配列"></param>
    /// <param name="A_pairID">ペアの共通ID、ペア情報が無い場合は空文字列</param>
    internal readonly record struct 読取証拠(string A_配列, string A_pairID)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 支持を数える際の独立性の単位
        /// </summary>
        public string A_独立性キー => string.IsNullOrEmpty(this.A_pairID) ? this.A_配列 : this.A_pairID;

        #endregion
    }
}
