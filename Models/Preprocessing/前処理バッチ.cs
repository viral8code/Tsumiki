using Tsumiki.Cores.Preprocessing;

namespace Tsumiki.Models.Preprocessing
{
    /// <summary>
    /// まとめて前処理するペアの束
    /// </summary>
    internal sealed class 前処理バッチ
    {
        #region プロパティ

        /// <summary>
        /// 入っているペア数
        /// </summary>
        public int A_件数 { get; set; }

        /// <summary>
        /// リード 1 の ID
        /// </summary>
        public string[] A_順リードID群 { get; } = new string[Preprocessor.C_前処理バッチサイズ];

        /// <summary>
        /// リード 2 の ID
        /// </summary>
        public string[] A_逆リードID群 { get; } = new string[Preprocessor.C_前処理バッチサイズ];

        /// <summary>
        /// リード 1 の配列
        /// </summary>
        public string[] A_順リード配列群 { get; } = new string[Preprocessor.C_前処理バッチサイズ];

        /// <summary>
        /// リード 2 の配列
        /// </summary>
        public string[] A_逆リード配列群 { get; } = new string[Preprocessor.C_前処理バッチサイズ];

        /// <summary>
        /// リード 1 のクオリティ
        /// </summary>
        public string[] A_順リード品質群 { get; } = new string[Preprocessor.C_前処理バッチサイズ];

        /// <summary>
        /// リード 2 のクオリティ
        /// </summary>
        public string[] A_逆リード品質群 { get; } = new string[Preprocessor.C_前処理バッチサイズ];

        /// <summary>
        /// 前処理の結果
        /// </summary>
        public ペア前処理結果[] A_結果群 { get; } = new ペア前処理結果[Preprocessor.C_前処理バッチサイズ];

        #endregion
    }
}
