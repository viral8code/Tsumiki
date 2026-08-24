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
                // GetNucleotideIDs は曖昧塩基対応のため List<int> を確保するが、
                // ContigMaker 側では badBase/revBadBase によって曖昧塩基を含む
                // 区間はそもそも KmerKey 化されない(呼ばれない)ため、
                // ここでは List 確保のない軽量な単一塩基変換で十分。
                // (曖昧塩基が来た場合は GetSimpleNucleotideID が InvalidBase(5) を
                //  返すが、そのようなケースは呼び出し元で事前に除外されている前提。)
                var val = (ulong)Util.GetSimpleNucleotideID(kmer[i]) - 1;
                // 32塩基ごとに同じ ulong 要素(2bit x 32 = 64bit)を共有するため、
                // 代入(=)ではなく OR(|=)で詰め込まないと、直前までに書き込んだ
                // 塩基の情報が上書きで消えてしまう。
                // (この不具合により、同じ ulong 要素に収まる k-mer 同士が
                //  実質「末尾の数文字だけで同一視される」形になっていた。)
                this.Data[index] |= val << shift;
            }
        }

        private KmerKey(ulong[] Data)
        {
            this.Data = Data;
        }

        public KmerKey ReverseComprement()
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
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