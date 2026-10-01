namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// 最終配列へペアを当て直して集めた、配列ごとの深さ・切れ端・読み・断片・外れ錨
    /// </summary>
    internal sealed class 継ぎ目の証拠
    {
        #region プロパティ

        /// <summary>
        /// 当てた配列群
        /// </summary>
        public required IReadOnlyList<string> A_配列群 { get; init; }

        /// <summary>
        /// 位置ごとの、当たったリードの数
        /// </summary>
        public required int[][] A_深さ { get; init; }

        /// <summary>
        /// 位置ごとの、一意に当たったリードの数
        /// </summary>
        public required int[][] A_一意な深さ { get; init; }

        /// <summary>
        /// 位置ごとの、そこから整列が始まり左側が切れていたリードの数 (長さは配列長 + 1)
        /// </summary>
        public required int[][] A_左の切れ端 { get; init; }

        /// <summary>
        /// 位置ごとの、そこで整列が終わり右側が切れていたリードの数 (長さは配列長 + 1)
        /// </summary>
        public required int[][] A_右の切れ端 { get; init; }

        /// <summary>
        /// 当たったリードの整列範囲 (始まりの順)
        /// </summary>
        public required (int A_開始, int A_終了)[][] A_読み { get; init; }

        /// <summary>
        /// 向き正しく組めたペアの断片 (始まりの順)
        /// </summary>
        public required (int A_開始, int A_終了)[][] A_断片 { get; init; }

        /// <summary>
        /// 相方が正しく組めていない一意なリードの位置と向き (右を向くなら +1、左を向くなら -1、位置の順)
        /// </summary>
        public required (int A_位置, int A_向き)[][] A_外れ錨 { get; init; }

        #endregion
    }
}
