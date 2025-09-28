using Tsumiki.Common;
using static Tsumiki.Common.Consts;

namespace Tsumiki.Utility
{
    internal class CountingDB : IDisposable
    {
        private const int MaxCount = 128 * 1024 * 1024;

        private readonly ByteArrayComparer _comparator;

        private readonly ByteArrayEqualityComparer _equalityComparator;

        private readonly string TempDirectory;

        private readonly string filePrefix;

        private int _fileCount;

        private BinaryWriter? _writer;

        public CountingDB(string tempDirectory)
        {
            this.filePrefix = Guid.NewGuid().ToString("N");
            this._comparator = new();
            this._equalityComparator = new();
            this.TempDirectory = tempDirectory;

            this._fileCount = 1;
            var fileName = Path.Combine(this.TempDirectory, $"{this.filePrefix}_{this._fileCount}");
            this._writer = new BinaryWriter(new BufferedStream(File.Open(fileName, FileMode.Create, FileAccess.Write)));
        }

        private void CreateNewFile()
        {
            this._writer?.Close();
            this._fileCount += 1;
            var newFileName = Path.Combine(this.TempDirectory, $"{this.filePrefix}_{this._fileCount}");
            this._writer = new BinaryWriter(new BufferedStream(File.Open(newFileName, FileMode.Create, FileAccess.Write)));
        }

        public void Add(Span<byte[]> key)
        {
            this.CreateByteArray(key, 0, new byte[(key.Length + 3) >> 2]);
        }

        private void CreateByteArray(Span<byte[]> key, int now, byte[] buffer)
        {
            if (now == key.Length)
            {
                this.Add(buffer);
                return;
            }
            var index = now >> 2;
            var shift = (3 - (now & 3)) << 1;
            foreach (var b in key[now])
            {
                var val = (byte)(b - 1 << shift);
                buffer[index] |= val;
                this.CreateByteArray(key, now + 1, buffer);
                buffer[index] &= (byte)~val;
            }
        }

        public void Add(Span<byte> key)
        {
            var arr = new byte[(key.Length + 3) / 4];
            var idx = 0;
            for (var i = 0; i < key.Length; i += 4)
            {
                var b = 0;
                for (var j = 0; j < 4; j++)
                {
                    var id = i + j < key.Length ? key[i + j] : NucleotideID.A;
                    b <<= 2;
                    b |= id - 1;
                }
                arr[idx++] = (byte)b;
            }

            this.Add(arr);
        }

        private void Add(byte[] values)
        {
            if (this._writer!.BaseStream.Position >= MaxCount)
            {
                this.CreateNewFile();
            }

            this._writer!.Write(values);
        }

        public string MergeAll()
        {
            this._writer!.Close();
            this._writer = null;
            var Length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            var mergedFileList = new List<string>();
            for (var i = 1; i <= this._fileCount; i++)
            {
                GC.Collect();
                var fileName = Path.Combine(this.TempDirectory, $"{this.filePrefix}_{i}");
                var dict = new Dictionary<byte[], ulong>(MaxCount / Length, this._equalityComparator);
                using (var reader = new BinaryReader(File.Open(fileName, FileMode.Open, FileAccess.Read)))
                {
                    while (Util.HasNext(reader))
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
                File.Delete(fileName);
                var mergedFileName = Path.Combine(this.TempDirectory, $"{this.filePrefix}_merged_{i}");
                using (var writer = new BinaryWriter(File.Open(mergedFileName, FileMode.Create, FileAccess.Write)))
                {
                    var arr = dict.ToArray();
                    Array.Sort(arr, (item1, item2) => this._comparator.Compare(item1.Key, item2.Key));
                    foreach (var kv in arr)
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
                var mergedFileName = Path.Combine(this.TempDirectory, $"{this.filePrefix}_merged_{index++}");
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
                            read1 = Util.HasNext(reader1) ? reader1.ReadBytes(Length) : null;
                            read2 = Util.HasNext(reader2) ? reader2.ReadBytes(Length) : null;
                        }
                        else if (result < 0)
                        {
                            var sum = reader1.ReadUInt64();
                            writer.Write(read1);
                            writer.Write(sum);
                            read1 = Util.HasNext(reader1) ? reader1.ReadBytes(Length) : null;
                        }
                        else
                        {
                            var sum = reader2.ReadUInt64();
                            writer.Write(read2);
                            writer.Write(sum);
                            read2 = Util.HasNext(reader2) ? reader2.ReadBytes(Length) : null;
                        }
                    }
                    while (read1 != null)
                    {
                        var sum = reader1.ReadUInt64();
                        writer.Write(read1);
                        writer.Write(sum);
                        read1 = Util.HasNext(reader1) ? reader1.ReadBytes(Length) : null;
                    }
                    while (read2 != null)
                    {
                        var sum = reader2.ReadUInt64();
                        writer.Write(read2);
                        writer.Write(sum);
                        read2 = Util.HasNext(reader2) ? reader2.ReadBytes(Length) : null;
                    }
                }
                File.Delete(file1);
                File.Delete(file2);
                mergedFileList.Add(mergedFileName);
            }
            return mergedFileList[0];
        }

        public void Dispose()
        {
            this._writer?.Close();
        }
    }
}
