using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Commons
{
    /// <summary>
    /// 実行中どこからでも参照する設定の置き場
    /// </summary>
    internal static class ConfigurationManager
    {
        #region プロパティ

        /// <summary>
        /// 実行時引数
        /// </summary>
        public static Parameters A_実行時引数 { get; set; } = new();

        /// <summary>
        /// kmer インデックス
        /// </summary>
        public static TrustedKmerIndex A_kmerインデックス { get; set; } = null!;

        /// <summary>
        /// 直近の <see cref="KmerCutoffSelector.V_解決_kmerカットオフ"/> で k-mer スペクトルに適合した混合モデル
        /// </summary>
        public static 混合スペクトル解析結果? A_スペクトルモデル { get; set; }

        #endregion
    }
}
