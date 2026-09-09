namespace Tsumiki.Utility
{
    internal class ByteArrayEqualityComparer : IEqualityComparer<byte[]>
    {
        /// <summary>
        /// 同じ中身のバイト列か
        /// </summary>
        /// <param name="x">比べるバイト列</param>
        /// <param name="y">比べるバイト列</param>
        /// <returns>同じなら true</returns>
        public bool Equals(byte[]? x, byte[]? y)
        {
            if (x == y)
            {
                return true;
            }

            if (x == null)
            {
                return y == null;
            }

            if (y == null)
            {
                return false;
            }

            if (x.Length != y.Length)
            {
                return false;
            }

            for (var i = 0; i < x.Length; i++)
            {
                if (x[i] != y[i])
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// バイト列のハッシュ値を返す
        /// </summary>
        /// <param name="obj">対象のバイト列</param>
        /// <returns>ハッシュ値</returns>
        public int GetHashCode(byte[] obj)
        {
            var l_ハッシュ = 17;
            foreach (var l_バイト in obj)
            {
                l_ハッシュ = (l_ハッシュ * 31) + l_バイト;
            }
            return l_ハッシュ;
        }
    }
}
