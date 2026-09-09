namespace Tsumiki.Utility
{
    internal class ByteArrayComparer : IComparer<byte[]>
    {
        /// <summary>
        /// バイト列を辞書式順序で比べる
        /// </summary>
        /// <param name="x">比べるバイト列</param>
        /// <param name="y">比べるバイト列</param>
        /// <returns>x が小さければ負、大きければ正、等しければ 0</returns>
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

            var l_比較長 = Math.Min(x.Length, y.Length);
            for (var i = 0; i < l_比較長; i++)
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
