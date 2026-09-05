using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer スペクトルから k-mer カットオフ(-kc)を自動選択する。
    /// 方針は「エラー由来が集合を支配しない範囲でできるだけ低く」で、
    /// その理由は <see cref="KmerHistogram.Get_推奨カットオフ"/> にある。
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
