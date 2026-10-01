namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// 継ぎ目の候補を、ペアの証拠がどれだけ効くかで分けた組
    /// </summary>
    internal enum 継ぎ目の組
    {
        /// <summary>
        /// 反復の間に N のギャップがある
        /// </summary>
        ギャップ,

        /// <summary>
        /// ギャップがなく、ペアで跨げると期待できる
        /// </summary>
        組で跨げる,

        /// <summary>
        /// ギャップがなく、ペアで跨げない (反復がインサートより長い、配列の端、深さが薄い)
        /// </summary>
        組で跨げない,
    }
}
