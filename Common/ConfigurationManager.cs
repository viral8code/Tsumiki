using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Common
{
    internal static class ConfigurationManager
    {
        public static Parameters A_実行時引数 { get; set; } = new();

        public static TrustedKmerIndex A_kmerインデックス { get; set; } = null!;

        /// <summary>
        /// 直近の KmerCutoffSelector.V_解決_kmerカットオフ で k-mer スペクトルに
        /// 適合した混合モデル<br/>
        /// CopyNumberEstimator の単一コピー基準値や GraphSimplifier の tip 判定の
        /// 「無条件に信頼する下限」が、カットオフと同じモデルから導いた値を共有するために公開する<br/>
        /// 適合に失敗した場合や -kc が明示指定された場合は null のままで、各所は自前の推定にフォールバックする
        /// </summary>
        public static 混合スペクトル解析結果? A_スペクトルモデル { get; set; }
    }
}
