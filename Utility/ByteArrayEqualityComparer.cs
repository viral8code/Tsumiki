namespace Tsumiki.Utility
{
    internal class ByteArrayEqualityComparer : IEqualityComparer<byte[]>
    {
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
