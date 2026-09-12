namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// FASTQ のクオリティ文字列を一定数サンプリングした結果
    /// </summary>
    /// <remarks>
    /// Phred オフセット (33 or 64) の推定に使う
    /// </remarks>
    /// <param name="A_最小ASCII"></param>
    /// <param name="A_最大ASCII"></param>
    /// <param name="A_標本リード数"></param>
    /// <param name="A_標本文字数"></param>
    internal readonly record struct Phred標本(int A_最小ASCII, int A_最大ASCII, int A_標本リード数, int A_標本文字数)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 標本全体を通して ASCII コードが一切変化しなかったか (実機のシーケンサ出力では通常あり得ない、人工的/ビニング済みのクオリティである可能性を示す)
        /// </summary>
        public bool A_Is一様 => this.A_標本文字数 > 0 && this.A_最小ASCII == this.A_最大ASCII;

        #endregion
    }
}
