using System.Text;
using Tsumiki.Common;

namespace Tsumiki.Utility
{
    internal class CountingBloomFilter(ulong bitSize) : IDisposable
    {
        private static readonly List<int> ShiftValues = [2, 3, 5];

        private readonly LargeBitArray _bitArray = new(bitSize);

        private readonly ulong _mod = bitSize;

        private readonly CountingDB _counter = new();

        public void Add(string read)
        {
            read = Util.CanonicalRead(read);
            this.SetHash(read);
        }

        public bool Contains(string read)
        {
            read = Util.CanonicalRead(read);
            var hashList = GetHashList(read);
            var hash = hashList.FirstOrDefault();
            return this._bitArray[hash];
        }

        public void Cutoff(ulong bounds)
        {
            var filePath = this._counter.MergeAll();
            var Length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            this._bitArray.Clear();
            using var reader = new BinaryReader(File.Open(filePath, FileMode.Open, FileAccess.Read));
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
                    this.Add(sb.ToString()[..ConfigurationManager.Arguments.Kmer]);
                }
            }
        }

        private void SetHash(string read)
        {
            var hashList = GetHashList(read);
            foreach (var hash in hashList)
            {
                this._bitArray[hash] = true;
            }
            this._counter.Add(read);
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
                            next.Add(hash % _mod);
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
            this._counter.Dispose();
        }

        private class LargeBitArray(ulong size)
        {
            private const int BitSize = 64;

            private readonly ulong[] _bitArray = new ulong[(size + BitSize - 1) / BitSize];

            public bool this[ulong index]
            {
                get
                {
                    var subIndex = index / BitSize;
                    var bitIndex = index % BitSize;
                    return (this._bitArray[subIndex] & (1UL << (int)bitIndex)) == 1;
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
                Array.Clear(this._bitArray);
            }
        }

        private class CountingDB : IDisposable
        {
            private const int MaxCount = 1_000_000;

            private readonly FastByteArrayComparer _comparator;

            private readonly string filePrefix;

            private int _fileCount;

            private int _count;

            private BinaryWriter _writer;

            public CountingDB()
            {
                this.filePrefix = Guid.NewGuid().ToString("N");
                this._comparator = new();

                this._fileCount = 1;
                this._count = 0;
                var fileName = $"{this.filePrefix}_{this._fileCount}";
                this._writer = new BinaryWriter(File.Open(fileName, FileMode.Create, FileAccess.Write));
            }

            private void CreateNewFile()
            {
                this._writer.Close();
                this._fileCount += 1;
                this._count = 0;
                var newFileName = $"{this.filePrefix}_{this._fileCount}";
                this._writer = new BinaryWriter(File.Open(newFileName, FileMode.Create, FileAccess.Write));
            }

            public void Add(string key)
            {
                List<byte[]> bytes = [[]];
                for (var i = 0; i < key.Length; i += 4)
                {
                    List<int> next = [0];
                    for (var j = 0; j < 4 && i + j < key.Length; j++)
                    {
                        List<int> subNext = [];
                        var ids = Util.GetNucleotideIDs(key[i + j]);
                        foreach (var id in ids)
                        {
                            foreach (var b in next)
                            {
                                subNext.Add((b << 2) | id);
                            }
                        }
                        next = subNext;
                    }

                    List<byte[]> nextBytes = [];
                    foreach (var bs in bytes)
                    {
                        foreach (var b in next)
                        {
                            nextBytes.Add([.. bs, (byte)b]);
                        }
                    }
                    bytes = nextBytes;
                }

                foreach (var bs in bytes)
                {
                    this.Add(bs);
                }
            }

            private void Add(byte[] values)
            {
                if (this._count == MaxCount)
                {
                    this.CreateNewFile();
                }
                this._count++;

                this._writer?.Write(values);
            }

            public string MergeAll()
            {
                this._writer.Close();
                this._writer = null!;
                var Length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
                var mergedFileList = new List<string>();
                for (var i = 1; i <= this._fileCount; i++)
                {
                    var fileName = $"{this.filePrefix}_{i}";
                    var dict = new SortedDictionary<byte[], ulong>(this._comparator);
                    using (var reader = new BinaryReader(File.Open(fileName, FileMode.Open, FileAccess.Read)))
                    {
                        while (reader.BaseStream.Position < reader.BaseStream.Length)
                        {
                            var read = reader.ReadBytes(Length);
                            if (dict.ContainsKey(read))
                            {
                                dict[read] += 1;
                            }
                            else
                            {
                                dict[read] = 1;
                            }
                        }
                    }
                    var mergedFileName = $"{this.filePrefix}_merged_{i}";
                    using (var writer = new BinaryWriter(File.Open(mergedFileName, FileMode.Create, FileAccess.Write)))
                    {
                        foreach (var kv in dict)
                        {
                            writer.Write(kv.Key);
                            writer.Write(kv.Value);
                        }
                    }
                    mergedFileList.Add(mergedFileName);
                }
                var index = this._fileCount + 1;
                while (mergedFileList.Count > 1)
                {
                    var file1 = mergedFileList[0];
                    var file2 = mergedFileList[1];
                    var mergedFileName = $"{this.filePrefix}_merged_{index++}";
                    mergedFileList.RemoveRange(0, 2);
                    using (var reader1 = new BinaryReader(File.Open(file1, FileMode.Open, FileAccess.Read)))
                    {
                        using var reader2 = new BinaryReader(File.Open(file2, FileMode.Open, FileAccess.Read));
                        using var writer = new BinaryWriter(File.Open(mergedFileName, FileMode.Create, FileAccess.Write));
                        var read1 = reader1.ReadBytes(Length);
                        var read2 = reader2.ReadBytes(Length);
                        while (read1 != null && read2 != null)
                        {
                            var result = this._comparator.Compare(read1, read2);
                            if (result == 0)
                            {
                                var sum = reader1.ReadUInt64() + reader2.ReadUInt64();
                                writer.Write(read1);
                                writer.Write(sum);
                                read1 = reader1.BaseStream.Position < reader1.BaseStream.Length ? reader1.ReadBytes(Length) : null;
                                read2 = reader2.BaseStream.Position < reader2.BaseStream.Length ? reader2.ReadBytes(Length) : null;
                            }
                            else if (result < 0)
                            {
                                var sum = reader1.ReadUInt64();
                                writer.Write(read1);
                                writer.Write(sum);
                                read1 = reader1.BaseStream.Position < reader1.BaseStream.Length ? reader1.ReadBytes(Length) : null;
                            }
                            else
                            {
                                var sum = reader2.ReadUInt64();
                                writer.Write(read2);
                                writer.Write(sum);
                                read2 = reader2.BaseStream.Position < reader2.BaseStream.Length ? reader2.ReadBytes(Length) : null;
                            }
                        }
                        while (read1 != null)
                        {
                            var sum = reader1.ReadUInt64();
                            writer.Write(read1);
                            writer.Write(sum);
                            read1 = reader1.BaseStream.Position < reader1.BaseStream.Length ? reader1.ReadBytes(Length) : null;
                        }
                        while (read2 != null)
                        {
                            var sum = reader2.ReadUInt64();
                            writer.Write(read2);
                            writer.Write(sum);
                            read2 = reader2.BaseStream.Position < reader2.BaseStream.Length ? reader2.ReadBytes(Length) : null;
                        }
                    }
                    mergedFileList.Add(mergedFileName);
                }
                return mergedFileList[0];
            }

            public void Dispose()
            {
                this._writer.Close();
            }

            private class FastByteArrayComparer : IComparer<byte[]>
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

                    var hx = GetHash(x);
                    var hy = GetHash(y);

                    if (hx != hy)
                    {
                        return hx.CompareTo(hy);
                    }

                    var len = Math.Min(x.Length, y.Length);
                    for (var i = 0; i < len; i++)
                    {
                        var cmp = x[i].CompareTo(y[i]);
                        if (cmp != 0)
                        {
                            return cmp;
                        }
                    }
                    return x.Length.CompareTo(y.Length);
                }

                private static int GetHash(byte[] data)
                {
                    unchecked
                    {
                        var hash = (int)2166136261;
                        foreach (var b in data)
                        {
                            hash ^= b;
                            hash *= 16777619;
                        }
                        return hash;
                    }
                }
            }
        }
    }
}
