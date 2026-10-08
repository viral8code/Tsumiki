namespace Tsumiki.Utilities
{
    /// <summary>
    /// 128 塩基までの k-mer の値 (上位と下位の組) を比べる比較器。既定のハッシュより速いハッシュを使う
    /// </summary>
    internal sealed class UInt128組比較器 : IEqualityComparer<(UInt128 A_上位, UInt128 A_下位)>
    {
        #region 公開フィールド

        /// <summary>
        /// 既定の比較器
        /// </summary>
        public static readonly UInt128組比較器 A_既定 = new();

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 上位と下位がともに等しいか
        /// </summary>
        /// <param name="p_x">比べる値</param>
        /// <param name="p_y">比べる値</param>
        /// <returns>両方の成分が等しければ true</returns>
        public bool Equals((UInt128 A_上位, UInt128 A_下位) p_x, (UInt128 A_上位, UInt128 A_下位) p_y)
        {
            return p_x.A_上位 == p_y.A_上位 && p_x.A_下位 == p_y.A_下位;
        }

        /// <summary>
        /// 上位と下位を混ぜ合わせたハッシュ値
        /// </summary>
        /// <param name="p_値">ハッシュを求める値</param>
        /// <returns>ハッシュ値</returns>
        public int GetHashCode((UInt128 A_上位, UInt128 A_下位) p_値)
        {
            var l_値 = (ulong)p_値.A_下位
                ^ ((ulong)(p_値.A_下位 >> 64) * 0x9E37_79B9_7F4A_7C15UL)
                ^ ((ulong)p_値.A_上位 * 0xC2B2_AE3D_27D4_EB4FUL)
                ^ ((ulong)(p_値.A_上位 >> 64) * 0x1656_67B1_9E37_79F9UL);
            l_値 *= 0xBF58_476D_1CE4_E5B9UL;
            l_値 ^= l_値 >> 31;
            return (int)l_値 ^ (int)(l_値 >> 32);
        }

        #endregion
    }

    /// <summary>
    /// リードへの配置候補 (配列番号・向き・対角線) を比べる比較器
    /// </summary>
    internal sealed class 配置候補比較器 : IEqualityComparer<(int A_配列番号, bool A_Is逆鎖, int A_対角線)>
    {
        #region 公開フィールド

        /// <summary>
        /// 既定の比較器
        /// </summary>
        public static readonly 配置候補比較器 A_既定 = new();

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 3 つの成分がすべて等しいか
        /// </summary>
        /// <param name="p_x">比べる候補</param>
        /// <param name="p_y">比べる候補</param>
        /// <returns>すべての成分が等しければ true</returns>
        public bool Equals((int A_配列番号, bool A_Is逆鎖, int A_対角線) p_x, (int A_配列番号, bool A_Is逆鎖, int A_対角線) p_y)
        {
            return p_x.A_配列番号 == p_y.A_配列番号 && p_x.A_Is逆鎖 == p_y.A_Is逆鎖 && p_x.A_対角線 == p_y.A_対角線;
        }

        /// <summary>
        /// 3 つの成分を混ぜ合わせたハッシュ値
        /// </summary>
        /// <param name="p_値">ハッシュを求める候補</param>
        /// <returns>ハッシュ値</returns>
        public int GetHashCode((int A_配列番号, bool A_Is逆鎖, int A_対角線) p_値)
        {
            var l_値 = ((ulong)(uint)p_値.A_配列番号 << 32) ^ (uint)p_値.A_対角線 ^ (p_値.A_Is逆鎖 ? 0x8000_0000_0000_0000UL : 0UL);
            l_値 *= 0x9E37_79B9_7F4A_7C15UL;
            return (int)(l_値 >> 32) ^ (int)l_値;
        }

        #endregion
    }
}
