using Tsumiki.Common;
using static Tsumiki.Common.Consts;

namespace Tsumiki.Utility
{
    internal class CountingDB : IDisposable
    {
        // フラッシュ前にメモリ上へ保持するエントリ数の上限を決める基準サイズ(バイト)。
        // 従来はここが「1ファイルあたりの生バイト数上限(128MB)」だったが、
        // 事前集約方式ではエントリ数(=ユニークなk-mer数)で管理する。
        private const int MaxCount = 256 * 1024 * 1024;

        // FileStream に渡すバッファサイズ。8バイト単位の細かい書き込みでも
        // システムコールが頻発しないよう大きめに確保する。
        private const int IoBufferSize = 16 * 1024 * 1024;

        private readonly ByteArrayComparer _comparator;

        private readonly ByteArrayEqualityComparer _equalityComparator;

        private readonly string TempDirectory;

        private readonly string filePrefix;

        private readonly int _length;

        private readonly int _flushThreshold;

        private int _fileCount;

        private Dictionary<byte[], ulong> _buffer;

        private readonly List<string> _flushedFiles = [];

        public CountingDB(string tempDirectory)
        {
            this.filePrefix = Guid.NewGuid().ToString("N");
            this._comparator = new();
            this._equalityComparator = new();
            this.TempDirectory = tempDirectory;
            this._length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            this._flushThreshold = Math.Max(1, MaxCount / Math.Max(1, this._length + sizeof(ulong)));
            this._buffer = new Dictionary<byte[], ulong>(this._flushThreshold, this._equalityComparator);
            this._fileCount = 0;
        }

        private static FileStream CreateWriteStream(string fileName)
        {
            return new FileStream(
                fileName,
                FileMode.Create,
                FileAccess.Write,
                FileShare.None,
                IoBufferSize,
                FileOptions.SequentialScan);
        }

        private static FileStream CreateReadStream(string fileName)
        {
            return new FileStream(
                fileName,
                FileMode.Open,
                FileAccess.Read,
                FileShare.Read,
                IoBufferSize,
                FileOptions.SequentialScan);
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
                var val = (byte)((b - 1) << shift);
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

        /// <summary>
        /// k-mer を1件登録する。従来はここで即ディスクに書き込んでいたが、
        /// メモリ上の Dictionary でカウントを集約することで、同一 k-mer の
        /// 再出現をディスク書き込みに変換しないようにする。
        /// 閾値に達したら整列済みの状態でディスクへフラッシュする。
        /// </summary>
        private void Add(byte[] values)
        {
            if (this._buffer.TryGetValue(values, out var count))
            {
                this._buffer[values] = count + 1;
            }
            else
            {
                this._buffer[values] = 1;
                if (this._buffer.Count >= this._flushThreshold)
                {
                    this.Flush();
                }
            }
        }

        /// <summary>
        /// メモリ上の集約済みカウントをキー順にソートしてディスクへ書き出す。
        /// フラッシュ後のファイルは常にソート済み・集約済みであるため、
        /// MergeAll 側では再集計(Dictionary への読み直し)が不要になる。
        /// </summary>
        private void Flush()
        {
            if (this._buffer.Count == 0)
            {
                return;
            }

            this._fileCount += 1;
            var fileName = Path.Combine(this.TempDirectory, $"{this.filePrefix}_{this._fileCount}");

            var arr = this._buffer.ToArray();
            Array.Sort(arr, (item1, item2) => this._comparator.Compare(item1.Key, item2.Key));

            using (var writer = new BinaryWriter(CreateWriteStream(fileName)))
            {
                foreach (var kv in arr)
                {
                    writer.Write(kv.Key);
                    writer.Write(kv.Value);
                }
            }

            this._flushedFiles.Add(fileName);
            this._buffer = new Dictionary<byte[], ulong>(this._flushThreshold, this._equalityComparator);
        }

        public string MergeAll()
        {
            // メモリ上に残っている未フラッシュ分を書き出す。
            this.Flush();

            var Length = this._length;
            var mergedFileList = new List<string>(this._flushedFiles);

            if (mergedFileList.Count == 0)
            {
                // 登録された k-mer が一件もない場合でも、空ファイルを返す。
                var emptyFileName = Path.Combine(this.TempDirectory, $"{this.filePrefix}_empty");
                using (CreateWriteStream(emptyFileName))
                {
                }
                return emptyFileName;
            }

            var index = this._fileCount + 1;
            while (mergedFileList.Count > 1)
            {
                GC.Collect();
                var file1 = mergedFileList[0];
                var file2 = mergedFileList[1];
                var mergedFileName = Path.Combine(this.TempDirectory, $"{this.filePrefix}_merged_{index++}");
                mergedFileList.RemoveRange(0, 2);
                using (var reader1 = new BinaryReader(CreateReadStream(file1)))
                {
                    using var reader2 = new BinaryReader(CreateReadStream(file2));
                    using var writer = new BinaryWriter(CreateWriteStream(mergedFileName));
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
            // 未フラッシュのデータは MergeAll 側で処理される想定だが、
            // MergeAll を呼ばずに破棄された場合に備えて残存ファイルを掃除する。
            foreach (var file in this._flushedFiles)
            {
                if (File.Exists(file))
                {
                    File.Delete(file);
                }
            }
        }
    }
}