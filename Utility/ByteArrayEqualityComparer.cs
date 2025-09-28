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
            var hash = 17;
            foreach (var b in obj)
            {
                hash = (hash * 31) + b;
            }
            return hash;
        }
    }
}
