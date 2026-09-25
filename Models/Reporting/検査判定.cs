namespace Tsumiki.Models.Reporting
{
    /// <summary>
    /// 1 つの検査項目の結果
    /// </summary>
    internal enum 検査判定
    {
        /// <summary>
        /// 検査項目を満たした
        /// </summary>
        合格,

        /// <summary>
        /// 検査項目を満たさなかった
        /// </summary>
        不合格,

        /// <summary>
        /// 判定に必要な材料が揃っておらず、合否を言えない
        /// </summary>
        判定不能,
    }
}
