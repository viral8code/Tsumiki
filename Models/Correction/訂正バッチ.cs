using Tsumiki.Cores.Preprocessing;

namespace Tsumiki.Models.Correction
{
    /// <summary>
    /// まとめて訂正するリードの束
    /// </summary>
    internal sealed class 訂正バッチ
    {
        #region プロパティ

        /// <summary>
        /// 入っているリード数
        /// </summary>
        public int A_件数 { get; set; }

        /// <summary>
        /// リードの ID
        /// </summary>
        public string[] A_ID群 { get; } = new string[ErrorCorrector.C_訂正バッチサイズ];

        /// <summary>
        /// リードのクオリティ
        /// </summary>
        public string[] A_クオリティ群 { get; } = new string[ErrorCorrector.C_訂正バッチサイズ];

        /// <summary>
        /// リードの塩基列
        /// </summary>
        public byte[][] A_塩基列群 { get; } = new byte[ErrorCorrector.C_訂正バッチサイズ][];

        /// <summary>
        /// 訂正の結果
        /// </summary>
        public 訂正結果[] A_結果群 { get; } = new 訂正結果[ErrorCorrector.C_訂正バッチサイズ];

        #endregion
    }
}
