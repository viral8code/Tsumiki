using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer スペクトルから k-mer カットオフ(-kc)を自動選択する。
    /// 方針は「エラー由来が集合を支配しない範囲でできるだけ低く」。
    ///
    /// まず2成分混合モデル(<see cref="KmerSpectrumMixtureModel"/>)の適合を試みる。
    /// これは谷の目視判定に頼らず事後誤り確率から閾値を導くため、低カバレッジなど
    /// 谷が視認できないデータでも働く。適合に失敗した場合のみ、谷検出
    /// (<see cref="KmerHistogram.Get_推奨カットオフ"/>)にフォールバックする。
    ///
    /// 混合モデルが適合できた場合、その単一コピー平均・信頼下限を
    /// ConfigurationManager.A_スペクトルモデル に公開し、CopyNumberEstimator の
    /// 単一コピー基準値と GraphSimplifier の tip 判定が同じモデルを共有できるようにする。
    /// </summary>
    internal static class KmerCutoffSelector
    {
        /// <summary>
        /// -kc が未指定の場合に限り、スペクトルから求めた値を適用する。
        /// ヒストグラムはカットオフ適用前に読む必要があるため統合ファイルを
        /// もう一度走査するが、明示指定時はこの走査自体を行わない。
        /// </summary>
        public static void V_解決_kmerカットオフ(Parameters p_引数, TrustedKmerIndex p_kmerインデックス)
        {
            // 前回(別のk、あるいはErrorCorrector用の一時インデックス)の適合結果を
            // 持ち越さない。適合に成功した場合のみ、この下で改めて設定し直す。
            ConfigurationManager.A_スペクトルモデル = null;

            if (p_引数.A_kmerカットオフが明示指定されたか)
            {
                return;
            }

            var l_ヒストグラム = p_kmerインデックス.Get_出現回数ヒストグラム();

            if (KmerSpectrumMixtureModel.Get_解析結果(l_ヒストグラム) is { } l_混合モデル)
            {
                ConfigurationManager.A_スペクトルモデル = l_混合モデル;
                if (l_混合モデル.A_カットオフ != p_引数.A_kmerカットオフ)
                {
                    p_引数.Set_推定kmerカットオフ(l_混合モデル.A_カットオフ);
                }
                Logger.V_出力(
                    メッセージID.kmerカットオフ_混合モデル, l_混合モデル.A_カットオフ, l_混合モデル.A_単一コピー平均, l_混合モデル.A_信頼下限, l_混合モデル.A_反復回数);
                return;
            }

            // フォールバック: 谷検出。混合モデルの適合に失敗するのは、データがこの
            // 2成分モデルにうまく当てはまらない(EMが収束しない、あるいは単一コピー
            // 成分と誤り成分を分離できない)場合。
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
    }
}
