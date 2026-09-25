using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// リード長から k 長を自動選択する
    /// </summary>
    internal static class KmerLengthSelector
    {
        #region 定数

        /// <summary>
        /// 自動選択する k 長のリード長に対する比
        /// </summary>
        private const double C_自動k長のリード長比 = 0.6D;

        /// <summary>
        /// k 長の自動選択に必要な最小リード長
        /// </summary>
        private const int C_自動k長に必要な最小リード長 = 32;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// -k 未指定時に、リード長から k を決めて適用する
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_リード長"></param>
        public static void V_解決_k長(Parameters p_引数, int? p_リード長)
        {
            if (p_リード長 is not { } l_リード長)
            {
                Logger.V_出力(メッセージID.k自動選択_リード長不明, p_引数.A_k長);
                return;
            }

            if (p_引数.A_Isk長明示指定)
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
        /// <returns></returns>
        public static int? Get_推奨k長(int p_リード長)
        {
            if (p_リード長 < C_自動k長に必要な最小リード長)
            {
                return null;
            }

            var l_候補 = (int)(p_リード長 * C_自動k長のリード長比);
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
