using Tsumiki.Common;
using static Tsumiki.Common.Consts;

namespace Tsumiki.Utility
{
    internal class CountingDB : IDisposable
    {
        // メモリ上に保持する k-mer カウントの総量(全シャード合計)の目安。
        //
        // 以前はこの値を「1シャードあたり」の基準として使っていたため、
        // -th 16 で実行すると 16 倍に膨らみ、実データ(100x)で
        // ピーク 12.5GB を消費してノートPCでは実行が困難だった。
        // 全シャードで分け合う総量として扱い、既定値も実測に基づいて下げる。
        //
        // 1エントリあたりの実消費は、キーの byte[](オブジェクトヘッダ24B +
        // 中身)と Dictionary のエントリ構造体を合わせて概ね 80B 前後。
        private const long DefaultTotalBudgetBytes = 768L * 1024 * 1024;

        private const int EstimatedBytesPerEntry = 80;

        // FileStream に渡すバッファサイズ。8バイト単位の細かい書き込みでも
        // システムコールが頻発しないよう大きめに確保する。
        private const int IoBufferSize = 1 << 20; // 1MB

        private readonly ByteArrayComparer _comparator;

        private readonly ByteArrayEqualityComparer _equalityComparator;

        private readonly string TempDirectory;

        private readonly string filePrefix;

        private readonly int _length;

        private readonly int _flushThreshold;

        private int _fileCount;

        private Dictionary<byte[], ulong> _buffer;

        private readonly List<string> _flushedFiles = [];

        /// <summary>
        /// shardCount には、同時に生きている CountingDB の総数を渡す。
        /// メモリ予算を等分するために使う。
        /// </summary>
        public CountingDB(string tempDirectory, int shardCount = 1)
        {
            this.filePrefix = Guid.NewGuid().ToString("N");
            this._comparator = new();
            this._equalityComparator = new();
            this.TempDirectory = tempDirectory;
            this._length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            // 総予算をシャード数で分け合う。shardCount は呼び出し側が渡す。
            var perShardBytes = DefaultTotalBudgetBytes / Math.Max(1, shardCount);
            this._flushThreshold = (int)Math.Max(1024, Math.Min(int.MaxValue, perShardBytes / EstimatedBytesPerEntry));
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
        public void AddPacked(byte[] values)
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
                    // BinaryReader.ReadBytes は EOF でも null ではなく長さ0の配列を
                    // 返すため、初回読み取りを保護しないと空ファイルを
                    // 「まだ中身がある」と誤認し、続く ReadUInt64 で破綻する。
                    // k-mer をハッシュでシャードへ振り分けるようにして以降、
                    // 空のシャードが普通に発生するようになったため必須。
                    var read1 = Util.HasNext(reader1) ? reader1.ReadBytes(Length) : null;
                    var read2 = Util.HasNext(reader2) ? reader2.ReadBytes(Length) : null;
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

            // フラッシュ済みファイルが1件のみだった場合、マージが一度も走らず
            // その元ファイル(_flushedFiles に登録済み)がそのまま返される。
            // 登録したままだと、この直後に Dispose() が呼ばれた際
            // _flushedFiles を掃除する処理で削除されてしまい、
            // 呼び出し元に返したパスが消える(FileNotFoundException の原因)。
            // 呼び出し元へ所有権を渡すため、返す前に登録を外しておく。
            var finalFile = mergedFileList[0];
            _ = this._flushedFiles.Remove(finalFile);
            return finalFile;
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

        /// <summary>
        /// 複数の CountingDB インスタンス(スレッドごとに用意したもの)が
        /// それぞれ MergeAll() で出力したソート済み・集約済みファイルを
        /// さらにペアワイズマージして1つに統合する。
        /// 並列読み込み(スレッドごとに独立した CountingDB を使う設計)で、
        /// 最後にワーカーごとの結果を1本化するために使用する。
        /// </summary>
        public static string MergeExternalFiles(string tempDirectory, List<string> filePaths)
        {
            var comparator = new ByteArrayComparer();
            var length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            var mergedFileList = new List<string>(filePaths);
            var prefix = Guid.NewGuid().ToString("N");
            var index = 1;

            if (mergedFileList.Count == 0)
            {
                var emptyFileName = Path.Combine(tempDirectory, $"{prefix}_empty");
                using (CreateWriteStream(emptyFileName))
                {
                }
                return emptyFileName;
            }

            while (mergedFileList.Count > 1)
            {
                GC.Collect();
                var file1 = mergedFileList[0];
                var file2 = mergedFileList[1];
                var mergedFileName = Path.Combine(tempDirectory, $"{prefix}_workermerge_{index++}");
                mergedFileList.RemoveRange(0, 2);
                using (var reader1 = new BinaryReader(CreateReadStream(file1)))
                {
                    using var reader2 = new BinaryReader(CreateReadStream(file2));
                    using var writer = new BinaryWriter(CreateWriteStream(mergedFileName));
                    // 上と同じ理由で、初回読み取りを HasNext で保護する。
                    var read1 = Util.HasNext(reader1) ? reader1.ReadBytes(length) : null;
                    var read2 = Util.HasNext(reader2) ? reader2.ReadBytes(length) : null;
                    while (read1 != null && read2 != null)
                    {
                        var result = comparator.Compare(read1, read2);
                        if (result == 0)
                        {
                            var sum = reader1.ReadUInt64() + reader2.ReadUInt64();
                            writer.Write(read1);
                            writer.Write(sum);
                            read1 = Util.HasNext(reader1) ? reader1.ReadBytes(length) : null;
                            read2 = Util.HasNext(reader2) ? reader2.ReadBytes(length) : null;
                        }
                        else if (result < 0)
                        {
                            var sum = reader1.ReadUInt64();
                            writer.Write(read1);
                            writer.Write(sum);
                            read1 = Util.HasNext(reader1) ? reader1.ReadBytes(length) : null;
                        }
                        else
                        {
                            var sum = reader2.ReadUInt64();
                            writer.Write(read2);
                            writer.Write(sum);
                            read2 = Util.HasNext(reader2) ? reader2.ReadBytes(length) : null;
                        }
                    }
                    while (read1 != null)
                    {
                        var sum = reader1.ReadUInt64();
                        writer.Write(read1);
                        writer.Write(sum);
                        read1 = Util.HasNext(reader1) ? reader1.ReadBytes(length) : null;
                    }
                    while (read2 != null)
                    {
                        var sum = reader2.ReadUInt64();
                        writer.Write(read2);
                        writer.Write(sum);
                        read2 = Util.HasNext(reader2) ? reader2.ReadBytes(length) : null;
                    }
                }
                File.Delete(file1);
                File.Delete(file2);
                mergedFileList.Add(mergedFileName);
            }

            return mergedFileList[0];
        }
    }
}