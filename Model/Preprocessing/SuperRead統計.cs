namespace Tsumiki.Model.Preprocessing
{
    /// <summary>
    /// ペアを合成リードへ統合する処理 1 回分の集計
    /// </summary>
    internal readonly record struct SuperRead統計(int A_総ペア数, int A_統合数, int A_重なり結合数, int A_曖昧で捨てた数)
    {
        #region カスタムプロパティ

        /// <summary>
        /// グラフ上の経路で橋渡しした数
        /// </summary>
        public int A_橋渡し数 => this.A_統合数 - this.A_重なり結合数;

        #endregion
    }
}
