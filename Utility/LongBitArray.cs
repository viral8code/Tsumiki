namespace Tsumiki.Utility
{
    internal class LongBitArray
    {
        private const int BitSize = 64;

        private const ulong MaxArrayLength = 1 << 25;

        private readonly ulong[][] _bitArray;

        private readonly ulong _length;

        public LongBitArray(ulong size)
        {
            this._length = size;

            var bitsPerArray = MaxArrayLength * BitSize;

            var arrayCount = (size + bitsPerArray - 1) / bitsPerArray;
            this._bitArray = new ulong[arrayCount][];
            for (ulong i = 0; i < arrayCount; i++)
            {
                var start = i * bitsPerArray;
                var end = Math.Min(size, start + bitsPerArray);
                var bits = end - start;
                var arrayLength = (bits + BitSize - 1) / BitSize;
                this._bitArray[i] = new ulong[arrayLength];
            }
        }

        public bool this[ulong index]
        {
            get
            {
                var (chunk, subIndex, bitIndex) = this.GetPosition(index);
                return (this._bitArray[chunk][subIndex] & (1UL << bitIndex)) != 0UL;
            }
            set
            {
                var (chunk, subIndex, bitIndex) = this.GetPosition(index);
                if (value)
                {
                    this._bitArray[chunk][subIndex] |= 1UL << bitIndex;
                }
                else
                {
                    this._bitArray[chunk][subIndex] &= ~(1UL << bitIndex);
                }
            }
        }

        public void Clear()
        {
            foreach (var arr in this._bitArray)
            {
                Array.Fill(arr, 0UL);
            }
        }

        private (int chunk, int subIndex, int bitIndex) GetPosition(ulong index)
        {
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(index, this._length);

            var bitsPerChunk = MaxArrayLength * BitSize;
            var chunk = (int)(index / bitsPerChunk);
            var offset = index % bitsPerChunk;

            var subIndex = (int)(offset / BitSize);
            var bitIndex = (int)(offset % BitSize);

            return (chunk, subIndex, bitIndex);
        }
    }
}
