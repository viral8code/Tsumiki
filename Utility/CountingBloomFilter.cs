using System.Text;
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

        public void Add(string read)
        {
            this.SetHash(read);
            this.SetHash(Util.ReverseComprement(read));
        }

        public bool Contains(string read)
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

        public List<string> Cutoff(ulong bounds)
        {
            var filePath = this._counter!.MergeAll();
            this._counter.Dispose();
            this._counter = null;
            var Length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            using (var reader = new BinaryReader(File.Open(filePath, FileMode.Open, FileAccess.Read)))
            {
                using var writer = new BinaryWriter(File.Open(KmerFileName, FileMode.Create, FileAccess.Write));
                while (reader.BaseStream.Position < reader.BaseStream.Length)
                {
                    var read = reader.ReadBytes(Length);
                    var count = reader.ReadUInt64();
                    if (count >= bounds)
                    {
                        StringBuilder sb = new();
                        foreach (var b in read)
                        {
                            _ = sb.Append(Util.ByteToNucleotideSequence(b));
                        }
                        var kmer = sb.ToString()[..ConfigurationManager.Arguments.Kmer];
                        if (kmer.Contains("N"))
                        {
                            Console.WriteLine(kmer);
                            Environment.Exit(0);
                        }
                        this.Add(kmer);
                        writer.Write(kmer);
                    }
                }
            }
            File.Delete(filePath);
            Console.WriteLine("Search First k-mer");
            List<string> kmers = [];
            using (var reader = new BinaryReader(File.Open(KmerFileName, FileMode.Open, FileAccess.Read)))
            {
                while (reader.BaseStream.Position < reader.BaseStream.Length)
                {
                    var read = reader.ReadString();
                    if (this.IsFirstKmer(read))
                    {
                        kmers.Add(read);
                    }
                }
            }
            File.Delete(KmerFileName);
            return kmers;
        }

        private bool IsFirstKmer(string kmer)
        {
            var count = 0;
            for (byte i = 1; i <= 4; i++)
            {
                var prevRead = Util.ByteToBase(i) + kmer[..^1];
                if (this.Contains(prevRead))
                {
                    count++;
                }
            }
            return count != 1;
        }

        private void SetHash(string read)
        {
            if (this._counter != null)
            {
                this._counter.Add(read);
            }
            else
            {
                var hashList = this.GetHashList(read);
                foreach (var hash in hashList)
                {
                    this._bitArray[hash] = true;
                }
            }
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
                    foreach (var id in ids.Select(v => (ulong)v))
                    {
                        foreach (var val in hashValues)
                        {
                            var hash = (val * (ulong)shift) + id;
                            next.Add(hash);
                        }
                    }
                    hashValues = next;
                }

                foreach (var hash in hashValues)
                {
                    _ = hashList.Add(hash);
                }
            }
            return [.. hashList.Select(num => num % this._mod)];
        }

        public void Dispose()
        {
            this._counter?.Dispose();
        }
    }
}
