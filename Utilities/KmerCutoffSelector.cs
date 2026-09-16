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
        /// <remarks>
        /// 混合モデルは単一コピーをポアソン分布で表すので、GC に偏ったライブラリのように同じ位置で繰り返す系統誤りが裾を引くと、谷よりずっと手前で誤りの範囲を打ち切る<br/>
        /// そのまま使うと誤りの k-mer が大量に残り、グラフが unitig 数の上限を超えて k ごと飛ぶ<br/>
        /// 谷とモデルがおおむね一致するデータでは、谷の取り方の揺れで結果を動かさないよう、明確に食い違うときだけ採る
        /// </remarks>
        private const ulong 谷を採る倍率 = 2UL;

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
            return p_谷 is { } l_谷 && l_谷 > p_モデルのカットオフ * 谷を採る倍率 ? l_谷 : p_モデルのカットオフ;
        }

        /// <summary>
        /// -kc が未指定の場合に限り、スペクトルから求めた値を適用する
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <remarks>
        /// ヒストグラムはカットオフ適用前に読む必要があるため統合ファイルをもう一度走査するが、明示指定時はこの走査自体を行わない
        /// </remarks>
        public static void V_解決_kmerカットオフ(Parameters p_引数, TrustedKmerIndex p_kmerインデックス)
        {
            // 前回 (別の k、あるいは ErrorCorrector 用の一時インデックス) の適合結果を
            // 持ち越さない
            // 適合に成功した場合のみ、この下で改めて設定し直す
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
                Logger.V_出力(メッセージID.kmerカットオフ_混合モデル, l_混合モデル.A_カットオフ, l_混合モデル.A_単一コピー平均, l_信頼下限, l_混合モデル.A_反復回数);
                if (l_カットオフ != l_混合モデル.A_カットオフ)
                {
                    Logger.V_出力_そのまま(FormattableString.Invariant($"[Info] k-mer cutoff raised from {l_混合モデル.A_カットオフ} to the spectrum valley {l_カットオフ} (the mixture model cut far below the valley)"));
                }
                return;
            }

            // フォールバック: 谷検出
            // 混合モデルの適合に失敗するのは、データがこの
            // 2 成分モデルにうまく当てはまらない (EM が収束しない、あるいは単一コピー
            // 成分と誤り成分を分離できない) 場合
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
