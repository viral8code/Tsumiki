using System.Runtime.InteropServices;
using Tsumiki.Common;

namespace Tsumiki.Utility
{
    internal class CountingBloomFilter(ulong bitSize) : IDisposable
    {
        private const string KmerFileName = "kmers";

        private static readonly List<int> ShiftValues = [1, 3, 4];

        private readonly LongBitArray _bitArray = new(bitSize);

        private readonly ulong _mod = bitSize;

        private CountingDB? _counter = new();

        public void Add(Span<byte[]> read)
        {
            this.Regist(read);
            this.Regist(Util.ReverseComprement(read));
        }

        public bool Contains(Span<byte> read)
        {
            var hashList = this.GetHashList(read);
            var flag = true;
            foreach (var hash in hashList)
            {
                flag &= this._bitArray[hash];
            }
            if (flag)
            {
                return true;
            }
            read = Util.ReverseComprement(read);
            hashList = this.GetHashList(read);
            flag = true;
            foreach (var hash in hashList)
            {
                flag &= this._bitArray[hash];
            }
            return flag;
        }

        public bool Contains(ulong read)
        {
            return this._bitArray[read];
        }

        public List<byte[]> Cutoff(ulong bounds)
        {
            var filePath = this._counter!.MergeAll();
            this._counter.Dispose();
            this._counter = null;
            var Length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            using (var reader = new BinaryReader(File.Open(filePath, FileMode.Open, FileAccess.Read)))
            {
                using var writer = new BinaryWriter(File.Open(KmerFileName, FileMode.Create, FileAccess.Write));
                ulong addedKmer = 0;
                ulong countKmer = 0;
                while (reader.BaseStream.Position < reader.BaseStream.Length)
                {
                    var read = reader.ReadBytes(Length);
                    var count = reader.ReadUInt64();
                    countKmer += 1;
                    if (count >= bounds)
                    {
                        addedKmer += 1;
                        List<byte> bytes = [];
                        foreach (var b in read)
                        {
                            bytes.AddRange(Util.ByteToNucleotideSequence(b));
                        }
                        var kmer = CollectionsMarshal.AsSpan(bytes)[..ConfigurationManager.Arguments.Kmer];
                        this.SetHash(kmer);
                        writer.Write(kmer);
                    }
                }
                Console.WriteLine("kmer count: " + countKmer);
                Console.WriteLine("good kmer: " + addedKmer);
            }
            File.Delete(filePath);
            Console.WriteLine("Search First k-mer");
            List<byte[]> kmers = [];
            using (var reader = new BinaryReader(File.Open(KmerFileName, FileMode.Open, FileAccess.Read)))
            {
                while (reader.BaseStream.Position < reader.BaseStream.Length)
                {
                    var read = reader.ReadBytes(ConfigurationManager.Arguments.Kmer);
                    if (this.IsFirstKmer(read))
                    {
                        kmers.Add(read);
                    }
                }
            }
            File.Delete(KmerFileName);
            return kmers;
        }

        private bool IsFirstKmer(Span<byte> kmer)
        {
            List<ulong> hashList = [];
            foreach (var shift in ShiftValues)
            {
                ulong hashValue = 0UL;

                foreach (var id in kmer[..^1])
                {
                    hashValue = (hashValue * (ulong)shift) + id;
                }

                hashList.Add(hashValue);
            }
            var exp = ConfigurationManager.Arguments.Kmer - 1;
            var count = 0;

            for (ulong i = 1; i <= 4; i++)
            {
                var isContains = true;
                for (var j = 0; j < hashList.Count; j++)
                {
                    var index = (hashList[j] + i * Util.Pow((ulong)ShiftValues[j], exp)) % this._mod;
                    isContains &= this._bitArray[index];
                }
                if (isContains)
                {
                    count++;
                }
            }
            return count != 1;
        }

        private void Regist(Span<byte[]> read)
        {
            this._counter?.Add(read);
        }

        private void SetHash(Span<byte> read)
        {
            var hashList = this.GetHashList(read);
            foreach (var hash in hashList)
            {
                this._bitArray[hash] = true;
            }
        }

        private List<ulong> GetHashList(Span<byte> read)
        {
            var hashList = new HashSet<ulong>();

            foreach (var shift in ShiftValues)
            {
                ulong hashValue = 0UL;

                foreach (var id in read)
                {
                    hashValue = (hashValue * (ulong)shift) + id;
                }

                _ = hashList.Add(hashValue);
            }

            return [.. hashList.Select(num => num % this._mod)];
        }

        public void Dispose()
        {
            this._counter?.Dispose();
        }
    }
}
