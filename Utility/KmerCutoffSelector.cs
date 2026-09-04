using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer スペクトルから k-mer カットオフ(-kc)を自動選択する。
    ///
    /// 選ぶ方針は「エラー由来の k-mer が集合を支配しない範囲で、できるだけ低く」。
    /// 谷で切ってはいけない理由と、この方針を裏付ける実測値は
    /// <see cref="KmerHistogram.Get_推奨カットオフ"/> に書いてある。
    ///
    /// 自動選択する意味は主にメモリにある。エラー由来の k-mer の絶対数は
    /// カバレッジとともに増え(実データの 100x では出現回数2だけで 177 万種類、
    /// 35x では 33 万種類)、固定値ではどちらかの帯で無駄を抱える。
    /// 品質はどちらの帯でも下限(2)と変わらないことを実測で確認している。
    /// </summary>
    internal static class KmerCutoffSelector
    {
        /// <summary>
        /// -kc が未指定の場合に限り、スペクトルから求めた値を適用する。
        ///
        /// ヒストグラムはカットオフを掛ける前に読む必要があるため、統合ファイルを
        /// もう一度走査する。-kc が明示指定されている場合はこの走査自体を行わない
        /// ので、ディスク I/O は従来どおりのまま。
        /// </summary>
        public static void V_解決_kmerカットオフ(Parameters p_引数, TrustedKmerIndex p_kmerインデックス)
        {
            if (p_引数.A_kmerカットオフが明示指定されたか)
            {
                return;
            }

            var l_ヒストグラム = p_kmerインデックス.Get_出現回数ヒストグラム();
            if (KmerHistogram.Get_推奨カットオフ(l_ヒストグラム) is not { } l_推奨値)
            {
                Console.WriteLine(
                    "[Info] Could not identify a clear histogram valley " +
                    $"(the spectrum may not be bimodal at this coverage); keeping -kc {p_引数.A_kmerカットオフ}.");
                return;
            }

            if (l_推奨値 == p_引数.A_kmerカットオフ)
            {
                return;
            }

            p_引数.Set_推定kmerカットオフ(l_推奨値);
            Console.WriteLine(
                $"[Info] k-mer cutoff auto-selected as {l_推奨値} from the k-mer spectrum. " +
                "Pass -kc explicitly to override.");
        }
    }
}
