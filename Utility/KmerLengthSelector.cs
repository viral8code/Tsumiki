using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// リード長から k 長を自動選択する。
    ///
    /// 自動選択する理由: 適正な k はリード長にほぼ比例するのに対し、既定値は
    /// 固定である。実測では 150bp のリードに対して k=31 と k=63 で
    /// unitig の N50 が 21,946 と 115,997、つまり5倍以上違った。
    /// k が短いとその長さを超える反復配列がすべて潰れてグラフが繋がらなくなる。
    /// 逆に k をリード長に近づけすぎると、1リードから取れる k-mer の本数
    /// (リード長 - k + 1)が減ってカバレッジが痩せる。
    /// </summary>
    internal static class KmerLengthSelector
    {
        /// <summary>
        /// -k 未指定時に、リード長から k を決めて適用する。
        /// 明示指定されている場合はユーザーの判断を尊重し、明らかに成立しない
        /// 場合(k がリード長以上)だけ警告する。
        /// </summary>
        public static void V_解決_k長(Parameters p_引数, int? p_リード長)
        {
            if (p_リード長 is not { } l_リード長)
            {
                Console.WriteLine(
                    $"[Info] Could not sample a read length; keeping -k {p_引数.A_k長}.");
                return;
            }

            if (p_引数.A_k長が明示指定されたか)
            {
                if (p_引数.A_k長 >= l_リード長)
                {
                    Console.WriteLine(
                        $"[Warning] -k {p_引数.A_k長} is not shorter than the observed read length ({l_リード長} bp); " +
                        "no k-mer can be extracted from a read. Lower -k.");
                }
                return;
            }

            if (Get_推奨k長(l_リード長) is not { } l_推奨値)
            {
                Console.WriteLine(
                    $"[Info] Observed read length ({l_リード長} bp) is too short to pick a k automatically; " +
                    $"keeping -k {p_引数.A_k長}.");
                return;
            }

            if (l_推奨値 == p_引数.A_k長)
            {
                return;
            }

            p_引数.Set_推定k長(l_推奨値);
            Console.WriteLine(
                $"[Info] k auto-selected as {l_推奨値} from the observed read length ({l_リード長} bp). " +
                "Pass -k explicitly to override.");
        }

        /// <summary>
        /// リード長に対する推奨 k 長。
        ///
        /// リード長の <see cref="Consts.自動k長のリード長比"/> 倍を目安に、
        /// 上限 <see cref="Consts.自動k長の上限"/> で頭を抑える。上限があるのは
        /// k が 64 を超えると 2bit パックが UInt128 に収まらなくなり、
        /// 高速経路から外れて実行時間もメモリも大きく悪化するため。
        ///
        /// 偶数の k は避ける。k が偶数だと k-mer 自身がその逆相補と一致しうる
        /// (回文)ため、正規形が縮退して隣接関係の判定が壊れる。
        /// </summary>
        public static int? Get_推奨k長(int p_リード長)
        {
            if (p_リード長 < Consts.自動k長に必要な最小リード長)
            {
                return null;
            }

            var l_候補 = (int)(p_リード長 * Consts.自動k長のリード長比);
            l_候補 = Math.Min(l_候補, Consts.自動k長の上限);
            if (l_候補 % 2 == 0)
            {
                l_候補 -= 1;
            }
            return l_候補;
        }
    }
}
