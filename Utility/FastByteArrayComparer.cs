namespace Tsumiki.Utility
{
    internal class FastByteArrayComparer : IComparer<byte[]>
    {
        public int Compare(byte[]? x, byte[]? y)
        {
            if (x == y)
            {
                return 0;
            }

            if (x == null)
            {
                return -1;
            }

            if (y == null)
            {
                return 1;
            }

            var hx = GetHash(x);
            var hy = GetHash(y);

            if (hx != hy)
            {
                return hx.CompareTo(hy);
            }

            var len = Math.Min(x.Length, y.Length);
            for (var i = 0; i < len; i++)
            {
                var cmp = x[i].CompareTo(y[i]);
                if (cmp != 0)
                {
                    return cmp;
                }
            }
            return x.Length.CompareTo(y.Length);
        }

        private static ulong GetHash(byte[] data)
        {
            unchecked
            {
                var hash = 2166136261UL;
                foreach (var b in data)
                {
                    hash ^= b;
                    hash *= 16777619;
                }
                return hash;
            }
        }
    }
}
