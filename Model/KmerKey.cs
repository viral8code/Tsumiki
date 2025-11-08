using Tsumiki.Common;

namespace Tsumiki.Model
{
    internal readonly struct KmerKey : IEquatable<KmerKey>
    {
        public readonly ulong[] Data;

        public KmerKey(ReadOnlySpan<char> kmer)
        {
            this.Data = new ulong[(kmer.Length + 31) >> 5];
            for (var i = 0; i < kmer.Length; i++)
            {
                var index = i >> 5;
                var shift = (31 ^ (i & 31)) << 1;
                var val = (ulong)Util.GetNucleotideIDs(kmer[i])[0] - 1;
                this.Data[index] = val << shift;
            }
        }

        private KmerKey(ulong[] Data)
        {
            this.Data = Data;
        }

        public KmerKey ReverseComprement()
        {
            int kmerLength = ConfigurationManager.Arguments.Kmer;
            var offSet = kmerLength & 31;
            var reversedKmer = new ulong[this.Data.Length];
            for (var i = 0; i < reversedKmer.Length; i++)
            {
                reversedKmer[i] = ~ReverseBit(this.Data[^(i + 1)]);
            }
            if (offSet > 0)
            {
                var shift = 31 - offSet;
                var temp = 0UL;
                for (var i = 1; i <= reversedKmer.Length; i++)
                {
                    var sub = reversedKmer[^i] >> offSet;
                    reversedKmer[^i] <<= shift;
                    reversedKmer[^i] |= temp;
                    temp = sub;
                }
            }
            return new KmerKey(reversedKmer);
        }

        private static ulong ReverseBit(ulong x)
        {
            x = ((x & 0x5555555555555555UL) << 1) | ((x >> 1) & 0x5555555555555555UL);
            x = ((x & 0x3333333333333333UL) << 2) | ((x >> 2) & 0x3333333333333333UL);
            x = ((x & 0x0F0F0F0F0F0F0F0FUL) << 4) | ((x >> 4) & 0x0F0F0F0F0F0F0F0FUL);
            x = ((x & 0x00FF00FF00FF00FFUL) << 8) | ((x >> 8) & 0x00FF00FF00FF00FFUL);
            x = ((x & 0x0000FFFF0000FFFFUL) << 16) | ((x >> 16) & 0x0000FFFF0000FFFFUL);
            x = (x << 32) | (x >> 32);
            return x;
        }

        public bool Equals(KmerKey other)
        {
            if (this.Data.Length != other.Data.Length)
            {
                return false;
            }

            for (var i = 0; i < this.Data.Length; i++)
            {
                if (this.Data[i] != other.Data[i])
                {
                    return false;
                }
            }
            return true;
        }

        public override bool Equals(object? obj)
        {
            return obj is KmerKey other && this.Equals(other);
        }

        public override int GetHashCode()
        {
            var hash = 1469598103934665603UL;
            foreach (var v in this.Data)
            {
                hash ^= v;
                hash *= 1099511628211UL;
            }
            return (int)(hash ^ (hash >> 32));
        }
    }

}
