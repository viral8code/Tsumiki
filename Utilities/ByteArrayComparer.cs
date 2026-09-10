namespace Tsumiki.Utilities
{
    internal class ByteArrayComparer : IComparer<byte[]>
    {
        /// <summary>
        /// バイト列を辞書式順序で比べる
        /// </summary>
        /// <param name="p_x">比べるバイト列</param>
        /// <param name="p_y">比べるバイト列</param>
        /// <returns>p_x が小さければ負、大きければ正、等しければ 0</returns>
        public int Compare(byte[]? p_x, byte[]? p_y)
        {
            if (p_x == p_y)
            {
                return 0;
            }

            if (p_x == null)
            {
                return -1;
            }

            if (p_y == null)
            {
                return 1;
            }

            var l_比較長 = Math.Min(p_x.Length, p_y.Length);
            for (var i = 0; i < l_比較長; i++)
            {
                if (p_x[i] != p_y[i])
                {
                    return p_x[i].CompareTo(p_y[i]);
                }
            }

            return p_x.Length.CompareTo(p_y.Length);
        }
    }
}
