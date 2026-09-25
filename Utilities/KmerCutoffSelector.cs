using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// k-mer スペクトルから k-mer カットオフ (-kc) を自動選択する
    /// </summary>
    internal static class KmerCutoffSelector
    {
        #region 定数

        /// <summary>
        /// 混合モデルのカットオフに対して、スペクトルの谷がこの倍率を超えたら谷を採る
        /// </summary>
        private const ulong C_谷を採る倍率 = 2UL;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 混合モデルのカットオフを、スペクトルの谷と明確に食い違うときだけ谷へ引き上げる
        /// </summary>
        /// <param name="p_モデルのカットオフ"></param>
        /// <param name="p_谷">スペクトルの谷、求められなければ null</param>
        /// <returns>採用するカットオフ</returns>
        internal static ulong Get_谷で補正したカットオフ(ulong p_モデルのカットオフ, ulong? p_谷)
        {
            return p_谷 is { } l_谷 && l_谷 > p_モデルのカットオフ * C_谷を採る倍率 ? l_谷 : p_モデルのカットオフ;
        }

        /// <summary>
        /// -kc が未指定の場合に限り、スペクトルから求めた値を適用する
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_kmerインデックス"></param>
        public static void V_解決_kmerカットオフ(Parameters p_引数, TrustedKmerIndex p_kmerインデックス)
        {
            ConfigurationManager.A_スペクトルモデル = null;

            if (p_引数.A_Iskmerカットオフ明示指定)
            {
                return;
            }

            var l_ヒストグラム = p_kmerインデックス.Get_出現回数ヒストグラム();

            if (KmerSpectrumMixtureModel.Get_解析結果(l_ヒストグラム) is { } l_混合モデル)
            {
                ConfigurationManager.A_スペクトルモデル = l_混合モデル;
                var l_谷 = KmerHistogram.Get_解析結果(l_ヒストグラム)?.A_谷;
                var l_カットオフ = Get_谷で補正したカットオフ(l_混合モデル.A_カットオフ, l_谷);
                if (l_カットオフ != p_引数.A_kmerカットオフ)
                {
                    p_引数.Set_推定kmerカットオフ(l_カットオフ);
                }

                var l_信頼下限 = l_混合モデル.A_信頼下限 == ulong.MaxValue ? "unknown (non-monotonic posterior)" : l_混合モデル.A_信頼下限.ToString(System.Globalization.CultureInfo.InvariantCulture);
                Logger.V_出力(メッセージID.kmerカットオフ_混合モデル, l_混合モデル.A_カットオフ, l_混合モデル.A_単一コピー平均, l_信頼下限, l_混合モデル.A_過分散, l_混合モデル.A_反復回数);
                if (l_カットオフ != l_混合モデル.A_カットオフ)
                {
                    Logger.V_出力_そのまま(FormattableString.Invariant($"[Info] k-mer cutoff raised from {l_混合モデル.A_カットオフ} to the spectrum valley {l_カットオフ} (the mixture model cut far below the valley)"));
                }

                return;
            }

            if (KmerHistogram.Get_推奨カットオフ(l_ヒストグラム) is not { } l_推奨値)
            {
                Logger.V_出力(メッセージID.kmerカットオフ_谷が不明, p_引数.A_kmerカットオフ);
                return;
            }

            if (l_推奨値 == p_引数.A_kmerカットオフ)
            {
                return;
            }

            p_引数.Set_推定kmerカットオフ(l_推奨値);
            Logger.V_出力(メッセージID.kmerカットオフ_スペクトル, l_推奨値);
        }

        #endregion
    }
}
