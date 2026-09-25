namespace Tsumiki.Models.Reporting
{
    /// <summary>
    /// 品質保証の段階
    /// </summary>
    internal enum 品質保証レベル
    {
        /// <summary>
        /// 配列を出力できた
        /// </summary>
        出力のみ = 0,

        /// <summary>
        /// k-mer の取りこぼしと出しすぎが許容範囲に収まる
        /// </summary>
        グラフ整合 = 1,

        /// <summary>
        /// 元リードを貼り直しても深度が途切れない
        /// </summary>
        マッピング整合 = 2,

        /// <summary>
        /// 未解決のギャップが残っていない
        /// </summary>
        ペア整合 = 3,

        /// <summary>
        /// 決めきれずに打ち切った分岐が残っていない
        /// </summary>
        接合点が支持済み = 4,

        /// <summary>
        /// 閉じ目まで元リードで裏付けられている
        /// </summary>
        完全長 = 5,
    }
}
