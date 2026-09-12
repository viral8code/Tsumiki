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
        /// <remarks>
        /// 警告・エラーと、完全性の判定・レポートの出力先
        /// </remarks>
        最小 = 0,

        /// <summary>
        /// 進行状況と結果
        /// </summary>
        /// <remarks>
        /// 内部の判断過程は出さない
        /// </remarks>
        標準 = 1,

        /// <summary>
        /// 内部の判断過程まで含めて全部
        /// </summary>
        詳細 = 2,
    }
}
