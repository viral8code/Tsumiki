using System.Runtime.InteropServices;
using Tsumiki.Common;

namespace Tsumiki.Utility
{
    internal class CountingBloomFilter : IDisposable
    {
        private readonly LongBitArray _bitArray;

        private readonly ulong _mod;

        private readonly string _tempDirectory;

        // 並列読み込み用に、ワーカースレッドの数だけ CountingDB を用意する。
        // 各スレッドは自分専用のインスタンスにのみ書き込むため、
        // Add 呼び出し時のロックが不要になる。
        private CountingDB[]? _counters;

        public CountingBloomFilter(ulong bitSize, string tempDirectory)
        {
            this._bitArray = new(bitSize);
            this._mod = bitSize;
            this._tempDirectory = tempDirectory;
            var workerCount = Math.Max(1, ConfigurationManager.Arguments.ThreadCount);
            this._counters = new CountingDB[workerCount];
            for (var i = 0; i < workerCount; i++)
            {
                this._counters[i] = new CountingDB(tempDirectory);
            }
        }

        /// <summary>
        /// 指定したワーカー番号(0 始まり)専用の CountingDB に登録する。
        /// 並列読み込み時、各スレッドは自分の workerIndex を固定して呼び出すことで
        /// スレッドセーフに(ロックなしで)k-mer を登録できる。
        /// </summary>
        public void Add(Span<byte[]> read, int workerIndex)
        {
            this.Regist(read, workerIndex);
        }

        /// <summary>
        /// 指定したワーカー番号(0 始まり)専用の CountingDB に登録する。
        /// </summary>
        public void Add(Span<byte> read, int workerIndex)
        {
            this.Regist(read, workerIndex);
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
            // 各ワーカーの CountingDB をそれぞれ MergeAll し、
            // 出来上がった複数のソート済みファイルをさらに1本にマージする。
            var mergedFiles = new List<string>();
            foreach (var counter in this._counters!)
            {
                mergedFiles.Add(counter.MergeAll());
                counter.Dispose();
            }
            this._counters = null;

            var filePath = CountingDB.MergeExternalFiles(this._tempDirectory, mergedFiles);

            var Length = (ConfigurationManager.Arguments.Kmer + 3) / 4;
            var kmerPath = Path.Combine(this._tempDirectory, Consts.KmerFileName);
            using (var reader = new BinaryReader(File.Open(filePath, FileMode.Open, FileAccess.Read)))
            {
                using var writer = new BinaryWriter(File.Open(kmerPath, FileMode.Create, FileAccess.Write));
                ulong addedKmer = 0;
                ulong countKmer = 0;
                while (Util.HasNext(reader))
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
            using (var reader = new BinaryReader(File.Open(kmerPath, FileMode.Open, FileAccess.Read)))
            {
                while (Util.HasNext(reader))
                {
                    var read = reader.ReadBytes(ConfigurationManager.Arguments.Kmer);
                    if (this.IsFirstKmer(read))
                    {
                        kmers.Add(read);
                    }
                }
            }
            File.Delete(kmerPath);
            return kmers;
        }

        /// <summary>
        /// 与えられた k-mer が unitig の開始点(入次数が1でない = 0個または2個以上の
        /// prefix 拡張が存在する)かどうかを判定する。
        /// 入次数は「kmer の末尾 k-1 文字」の先頭に任意の1塩基を付加した k-mer が
        /// Bloom filter に存在するかで数える。逆相補鎖側からの接続も存在としてカウントする。
        /// </summary>
        private bool IsFirstKmer(Span<byte> kmer)
        {
            var count = this.CountInEdges(kmer);
            return count != 1;
        }

        /// <summary>
        /// kmer への入次数（前方に接続しうる異なる1塩基拡張の数）を数える。
        /// 順鎖・逆相補鎖の両方を試す点は Contains と同様。
        /// </summary>
        private int CountInEdges(Span<byte> kmer)
        {
            // kmer の末尾 k-1 文字（先頭の1文字を除いたもの）をベースに、
            // 先頭に1塩基 i を追加したときのハッシュを計算する。
            var suffixHashList = new List<ulong>();
            foreach (var shift in Consts.ShiftValues)
            {
                var hashValue = 0UL;
                foreach (var id in kmer[1..])
                {
                    hashValue = (hashValue * (ulong)shift) + id;
                }
                suffixHashList.Add(hashValue);
            }

            var exp = ConfigurationManager.Arguments.Kmer - 1;
            var count = 0;

            for (ulong i = 1; i <= 4; i++)
            {
                var isContains = true;
                for (var j = 0; j < suffixHashList.Count; j++)
                {
                    // 先頭に付加する1塩基 i は最上位桁(shift^exp)に乗る
                    var index = ((i * Util.Pow((ulong)Consts.ShiftValues[j], exp)) + suffixHashList[j]) % this._mod;
                    isContains &= this._bitArray[index];
                }

                if (!isContains)
                {
                    // 順鎖側で見つからない場合、逆相補鎖側での存在も確認する。
                    // 「先頭に塩基 i を追加した k-mer」の逆相補鎖は
                    // 「kmer の逆相補鎖の末尾に comp(i) を追加した k-mer」に等しい。
                    var candidate = new byte[kmer.Length];
                    kmer.CopyTo(candidate);
                    candidate[0] = (byte)i;
                    var revCandidate = Util.ReverseComprement(candidate);
                    isContains = this.ContainsExact(revCandidate);
                }

                if (isContains)
                {
                    count++;
                }
            }

            return count;
        }

        /// <summary>
        /// 与えられた k-mer が（逆相補は取らずに）そのままの向きで
        /// Bloom filter に存在するかを判定する。
        /// </summary>
        private bool ContainsExact(Span<byte> read)
        {
            var hashList = this.GetHashList(read);
            var flag = true;
            foreach (var hash in hashList)
            {
                flag &= this._bitArray[hash];
            }
            return flag;
        }

        private void Regist(Span<byte[]> read, int workerIndex)
        {
            this._counters?[workerIndex].Add(read);
        }

        private void Regist(Span<byte> read, int workerIndex)
        {
            this._counters?[workerIndex].Add(read);
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

            foreach (var shift in Consts.ShiftValues)
            {
                var hashValue = 0UL;

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
            if (this._counters != null)
            {
                foreach (var counter in this._counters)
                {
                    counter.Dispose();
                }
            }
        }
    }
}