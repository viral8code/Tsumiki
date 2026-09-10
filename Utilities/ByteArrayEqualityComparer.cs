namespace Tsumiki.Utilities
{
    internal class ByteArrayEqualityComparer : IEqualityComparer<byte[]>
    {
        /// <summary>
        /// 同じ中身のバイト列か
        /// </summary>
        /// <param name="p_x">比べるバイト列</param>
        /// <param name="p_y">比べるバイト列</param>
        /// <returns>同じなら true</returns>
        public bool Equals(byte[]? p_x, byte[]? p_y)
        {
            if (p_x == p_y)
            {
                return true;
            }

            if (p_x == null)
            {
                return p_y == null;
            }

            if (p_y == null)
            {
                return false;
            }

            if (p_x.Length != p_y.Length)
            {
                return false;
            }

            for (var i = 0; i < p_x.Length; i++)
            {
                if (p_x[i] != p_y[i])
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// バイト列のハッシュ値を返す
        /// </summary>
        /// <param name="p_obj">対象のバイト列</param>
        /// <returns>ハッシュ値</returns>
        public int GetHashCode(byte[] p_obj)
        {
            var l_ハッシュ = 17;
            foreach (var l_バイト in p_obj)
            {
                l_ハッシュ = (l_ハッシュ * 31) + l_バイト;
            }
            return l_ハッシュ;
        }
    }
}
