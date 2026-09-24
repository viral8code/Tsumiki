namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// 画面へ出す量の段階
    /// </summary>
    internal enum ログ水準
    {
        /// <summary>
        /// 結論だけ
        /// </summary>
        最小 = 0,

        /// <summary>
        /// 進行状況と結果
        /// </summary>
        標準 = 1,

        /// <summary>
        /// 内部の判断過程まで含めて全部
        /// </summary>
        詳細 = 2,
    }
}
