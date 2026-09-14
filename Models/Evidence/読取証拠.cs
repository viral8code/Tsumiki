namespace Tsumiki.Models.Evidence
{
    /// <summary>
    /// 局所アセンブリへ回収した 1 本のリードと、その由来
    /// </summary>
    /// <param name="A_配列"></param>
    /// <param name="A_pairID">ペアの共通ID、ペア情報が無い場合は空文字列</param>
    /// <remarks>
    /// 由来 (pair ID) を保持するのは、mate1/mate2 のように同一分子から出た 2 本を
    /// 独立な 2 件の支持として重複加算しないため
    /// </remarks>
    internal readonly record struct 読取証拠(string A_配列, string A_pairID)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 支持を数える際の独立性の単位
        /// </summary>
        /// <remarks>
        /// pair ID があればそれを使い、同一分子由来の複数本をまとめて 1 件に数える<br/>
        /// pair 情報が無い (シングルエンド等) 場合は配列そのものを単位にし、
        /// 同一配列の重複 (PCR 重複の疑いを含む) を 1 件にまとめる
        /// </remarks>
        public string A_独立性キー => string.IsNullOrEmpty(this.A_pairID) ? this.A_配列 : this.A_pairID;

        #endregion
    }
}
