namespace Tsumiki.Utility
{
    internal class ByteArrayComparer : IComparer<byte[]>
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

            var len = Math.Min(x.Length, y.Length);
            for (var i = 0; i < len; i++)
            {
                if (x[i] != y[i])
                {
                    return x[i].CompareTo(y[i]);
                }
            }

            return x.Length.CompareTo(y.Length);
        }
    }
}
