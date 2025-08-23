namespace Tsumiki.Utility
{
    internal class LongBitArray(ulong size)
    {
        private const int BitSize = 64;

        private readonly ulong[] _bitArray = new ulong[(size + BitSize - 1) / BitSize];

        public bool this[ulong index]
        {
            get
            {
                var subIndex = index / BitSize;
                var bitIndex = index % BitSize;
                return (this._bitArray[subIndex] & (1UL << (int)bitIndex)) > 0UL;
            }
            set
            {
                var subIndex = index / BitSize;
                var bitIndex = index % BitSize;
                if (value)
                {
                    this._bitArray[subIndex] |= 1UL << (int)bitIndex;
                }
                else
                {
                    this._bitArray[subIndex] &= ~(1UL << (int)bitIndex);
                }
            }
        }

        public void Clear()
        {
            Array.Fill(this._bitArray, 0UL);
        }
    }
}
