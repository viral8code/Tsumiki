using System.Text;
using Tsumiki.Common;

namespace Tsumiki.Utility
{
    internal class CountingBloomFilter(ulong bitSize) : IDisposable
    {
        private static readonly List<int> ShiftValues = [2, 3, 5];

        private readonly LongBitArray _bitArray = new(bitSize);

        private readonly ulong _mod = bitSize;

        private CountingDB _counter = new();

        public void Add(string read)
        {
            this.SetHash(read);
            this.SetHash(Util.ReverseComprement(read));
        }

        public bool Contains(string read)
        {
            var hashList = this.GetHashList(read);
            bool flag = true;
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

        public List<string> Cutoff(ulong bounds)
        {
            var filePath = this._counter.MergeAll();
            this._counter.Dispose();
            this._counter = null!;
            var Length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            this._bitArray.Clear();
            List<string> kmers = [];
            using (var reader = new BinaryReader(File.Open(filePath, FileMode.Open, FileAccess.Read)))
            {
                while (reader.BaseStream.Position < reader.BaseStream.Length)
                {
                    var read = reader.ReadBytes(Length);
                    var count = reader.ReadUInt64();
                    if (count > bounds)
                    {
                        StringBuilder sb = new();
                        foreach (var b in read)
                        {
                            _ = sb.Append(Util.ByteToNucleotideSequence(b));
                        }
                        var kmer = sb.ToString()[..ConfigurationManager.Arguments.Kmer];
                        this.Add(kmer);
                        kmers.Add(kmer);
                    }
                }
            }
            File.Delete(filePath);
            return kmers;
        }

        private void SetHash(string read)
        {
            var hashList = this.GetHashList(read);
            foreach (var hash in hashList)
            {
                this._bitArray[hash] = true;
            }
            this._counter?.Add(read);
        }

        private List<ulong> GetHashList(string read)
        {
            var hashList = new HashSet<ulong>();
            foreach (var shift in ShiftValues)
            {
                var hashValues = new List<ulong>
                {
                    0UL
                };

                foreach (var c in read)
                {
                    var ids = Util.GetNucleotideIDs(c);
                    var next = new List<ulong>();
                    foreach (var id in ids)
                    {
                        foreach (var val in hashValues)
                        {
                            var hash = (val << shift) | (uint)id;
                            next.Add(hash % this._mod);
                        }
                    }
                    hashValues = next;
                }

                foreach (var hash in hashValues)
                {
                    _ = hashList.Add(hash);
                }
            }
            return [.. hashList];
        }

        public void Dispose()
        {
            this._counter?.Dispose();
        }
    }
}
