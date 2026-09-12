using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// リード長から k 長を自動選択する
    /// </summary>
    /// <remarks>
    /// k が短いとその長さを超える反復配列が潰れ、長すぎると 1 リードから取れる k-mer の本数 (リード長 - k + 1) が減ってカバレッジが痩せる
    /// </remarks>
    internal static class KmerLengthSelector
    {
        #region 公開メソッド

        /// <summary>
        /// -k 未指定時に、リード長から k を決めて適用する
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 明示指定されている場合はユーザーの判断を尊重し、明らかに成立しない場合 (k がリード長以上) だけ警告する
        /// </remarks>
        public static void V_解決_k長(Parameters p_引数, int? p_リード長)
        {
            if (p_リード長 is not { } l_リード長)
            {
                Logger.V_出力(メッセージID.k自動選択_リード長不明, p_引数.A_k長);
                return;
            }

            if (p_引数.A_k長が明示指定されたか)
            {
                if (p_引数.A_k長 >= l_リード長)
                {
                    Logger.V_出力(メッセージID.k自動選択_kが長すぎる, p_引数.A_k長, l_リード長);
                }
                return;
            }

            if (Get_推奨k長(l_リード長) is not { } l_推奨値)
            {
                Logger.V_出力(メッセージID.k自動選択_リード長が短い, l_リード長, p_引数.A_k長);
                return;
            }

            if (l_推奨値 == p_引数.A_k長)
            {
                return;
            }

            p_引数.Set_推定k長(l_推奨値);
            Logger.V_出力(メッセージID.k自動選択, l_推奨値, l_リード長);
        }

        /// <summary>
        /// リード長に対する推奨 k 長
        /// </summary>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 上限を設けるのは、k が 64 を超えると 2 bit パックが UInt128 に収まらず高速経路から外れるため<br/>
        /// 偶数を避けるのは、k-mer 自身がその逆相補と一致しうる (回文) と正規形が縮退して隣接判定が壊れるため
        /// </remarks>
        /// <returns></returns>
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

        #endregion
    }
}
