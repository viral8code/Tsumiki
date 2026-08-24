using System.Collections.Concurrent;
using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            try
            {
                Run(args);
            }
            catch (Exception ex)
            {
                Logger.PrintError("Unhandled Tsumiki's method", ex);
            }
            finally
            {
                var tempDir = Path.Combine(Environment.CurrentDirectory, ConfigurationManager.Arguments.TempDirectory);
                if (Directory.Exists(tempDir))
                {
                    Directory.Delete(tempDir, true);
                }
            }
        }

        private static void Run(string[] args)
        {
            if (args.Length == 0)
            {
                Console.WriteLine(Consts.DetailsText);
                Environment.Exit(0);
            }

            var param = ArgumentsReader.ReadArguments(args);
            ConfigurationManager.Arguments = param;

            if (param.IsHelpMode)
            {
                Console.WriteLine(Consts.HelpText);
                Environment.Exit(0);
            }

            Console.WriteLine(param);

            Logger.PrintTimeStamp();

            var tempDir = Path.Combine(Environment.CurrentDirectory, param.TempDirectory);

            if (Path.Exists(tempDir))
            {
                Console.WriteLine($"{param.TempDirectory} already exists");
                Console.WriteLine("Please check path!");
                Environment.Exit(0);
            }

            _ = Directory.CreateDirectory(tempDir);

            Console.WriteLine("Start construction Bloom filter");
            var bloomFilter = new CountingBloomFilter(param.RowBitSize, tempDir);
            ConfigurationManager.BloomFilter = bloomFilter;

            if (string.IsNullOrWhiteSpace(param.ReadPath2))
            {
                Console.WriteLine("Loading File");
            }
            else
            {
                Console.WriteLine("Loading File1");
            }

            if (param.AllowAmbiguousBases)
            {
                LoadReadFileToBloomFilterWithAmbiguity(param.ReadPath1, bloomFilter);
            }
            else
            {
                LoadReadFileToBloomFilterIgnoreAmbiguity(param.ReadPath1, bloomFilter);
            }

            if (!string.IsNullOrWhiteSpace(param.ReadPath2))
            {
                Console.WriteLine("Loading File2");
                if (param.AllowAmbiguousBases)
                {
                    LoadReadFileToBloomFilterWithAmbiguity(param.ReadPath2, bloomFilter);
                }
                else
                {
                    LoadReadFileToBloomFilterIgnoreAmbiguity(param.ReadPath2, bloomFilter);
                }
            }

            Logger.PrintTimeStamp();

            Console.WriteLine("Fix Bloom filter");
            var initKmers = bloomFilter.Cutoff(param.KmerCutoff);

            Logger.PrintTimeStamp();

            Console.WriteLine("Make unitigs");

            var unitigMaker = new UnitigMaker(bloomFilter);
            HashSet<string> unitigSet = [];
            var id = 1;
            using (var writer = new FastaWriter(Consts.UnitigFileName))
            {
                foreach (var kmer in initKmers)
                {
                    var unitig = unitigMaker.MakeUnitig(kmer);
                    if (unitigSet.Contains(unitig.Sequence) || unitigSet.Contains(Util.ReverseComprement(unitig.Sequence)))
                    {
                        continue;
                    }
                    _ = unitigSet.Add(unitig.Sequence);
                    _ = unitigSet.Add(Util.ReverseComprement(unitig.Sequence));
                    writer.Write(id++, unitig.Sequence);
                    if (id > Consts.MaximumUnitigCount)
                    {
                        break;
                    }
                }
            }

            Logger.PrintTimeStamp();

            if (id > Consts.MaximumUnitigCount)
            {
                Console.WriteLine("""
                    
                    This genome is too complex to assembly...
                    Please adjust the parameters!
                    
                    """);
                Environment.Exit(0);
            }

            Console.WriteLine("Map reads to unitigs");

            var contigMaker = new ContigMaker(Consts.UnitigFileName);

            if (string.IsNullOrWhiteSpace(param.ReadPath2))
            {
                Console.WriteLine(param.ReadPath1);
                contigMaker.MappingRead(param.ReadPath1);
            }
            else
            {
                // ペアエンドの場合、read1/read2 を同時に読み進めて
                // インサートサイズによる隣接検出も行う。
                Console.WriteLine(param.ReadPath1);
                Console.WriteLine(param.ReadPath2);
                contigMaker.MappingPairedReads(param.ReadPath1, param.ReadPath2);
            }

            Logger.PrintTimeStamp();

            Console.WriteLine("unite unitigs");

            contigMaker.UniteContigs(Consts.ContigFileName, 0.8m, 10);

            Console.WriteLine("Maked contigs");

            Logger.PrintTimeStamp();

            Console.WriteLine("開発中！");

            Logger.PrintTimeStamp();
        }

        // 曖昧塩基を許容する経路は現状シングルスレッドのまま(workerIndex固定)。
        // 呼ばれる頻度が低い想定のため、並列化の優先度を下げている。
        private static void LoadReadFileToBloomFilterWithAmbiguity(string filePath, CountingBloomFilter bloomFilter)
        {
            ulong count = 0;
            ulong mult = 0;

            using var reader = new FastqReader(filePath);
            while (reader.HasNext())
            {
                var readData = reader.NextRead();
                if (readData.Read!.Count < ConfigurationManager.Arguments.Kmer)
                {
                    continue;
                }
                var readSpan = CollectionsMarshal.AsSpan(readData.Read);
                for (var i = 0; i < readData.Quality.Length; i++)
                {
                    if (readData.Quality[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        readSpan[i] = [Consts.NucleotideID.A, Consts.NucleotideID.C, Consts.NucleotideID.G, Consts.NucleotideID.T];
                    }
                }
                bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer], 0);
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
                {
                    bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer), 0);
                }
                readSpan = Util.ReverseComprement(readSpan);
                bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer], 0);
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
                {
                    bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer), 0);
                }
                if (++count == Consts.ProgressLogInterval)
                {
                    Console.WriteLine((++mult * Consts.ProgressLogInterval) + " reads Loaded");
                    count = 0;
                }
            }
            var fileName = Path.GetFileName(filePath);
            Console.WriteLine($"Loaded {(mult * Consts.ProgressLogInterval) + count} reads from {fileName}");
        }

        /// <summary>
        /// FASTQ を1本のスレッドで順に読み進めつつ(ディスクI/Oはシーケンシャルなまま)、
        /// 読み取ったリードを BlockingCollection 経由でワーカースレッド群に配る
        /// プロデューサー/コンシューマ方式。各ワーカーは自分専用の workerIndex を使って
        /// CountingBloomFilter.Add を呼ぶため、ロックなしで並列に k-mer を登録できる。
        /// (CPU 側の処理 = k-mer 分解・品質判定・Dictionary 更新 が重い場合に効果が出る。
        ///  ディスクI/O自体が律速の場合は改善が小さい点に注意。)
        /// </summary>
        private static void LoadReadFileToBloomFilterIgnoreAmbiguity(string filePath, CountingBloomFilter bloomFilter)
        {
            var threadCount = Math.Max(1, ConfigurationManager.Arguments.ThreadCount);
            ulong totalCount = 0;
            ulong mult = 0;
            var countLock = new object();

            // キューの深さはスレッド数に応じて適当な余裕を持たせる。
            // 大きすぎるとメモリを圧迫するため、ある程度で背圧をかける。
            using var queue = new BlockingCollection<ReadData>(boundedCapacity: threadCount * 64);

            var workers = new Task[threadCount];
            for (var w = 0; w < threadCount; w++)
            {
                var workerIndex = w;
                workers[w] = Task.Run(() =>
                {
                    foreach (var readData in queue.GetConsumingEnumerable())
                    {
                        ProcessRead(readData, bloomFilter, workerIndex);

                        var shouldLog = false;
                        ulong logValue = 0;
                        lock (countLock)
                        {
                            totalCount++;
                            if (totalCount % Consts.ProgressLogInterval == 0)
                            {
                                mult++;
                                shouldLog = true;
                                logValue = mult * Consts.ProgressLogInterval;
                            }
                        }
                        if (shouldLog)
                        {
                            Console.WriteLine(logValue + " reads Loaded");
                        }
                    }
                });
            }

            using (var reader = new FastqReader(filePath))
            {
                while (reader.HasNext())
                {
                    var readData = reader.NextReadSimple();
                    queue.Add(readData);
                }
            }
            queue.CompleteAdding();

            Task.WaitAll(workers);

            var fileName = Path.GetFileName(filePath);
            Console.WriteLine($"Loaded {totalCount} reads from {fileName}");
        }

        /// <summary>
        /// 1リード分の k-mer 抽出・品質フィルタリング・Bloom filter 登録を行う。
        /// LoadReadFileToBloomFilterIgnoreAmbiguity の元の逐次実装と同一のロジックを、
        /// ワーカースレッドから呼び出せる形に切り出したもの。
        /// </summary>
        private static void ProcessRead(ReadData readData, CountingBloomFilter bloomFilter, int workerIndex)
        {
            var simpleRead = readData.SimpleRead!;
            if (simpleRead.Length < ConfigurationManager.Arguments.Kmer)
            {
                return;
            }
            var badQualityCount = 0;
            var qualitySpan = readData.Quality.ToCharArray().AsSpan();
            var readSpan = simpleRead.AsSpan();
            for (var i = 0; i < ConfigurationManager.Arguments.Kmer; i++)
            {
                if (readSpan[i] == Consts.InvalidBase ||
                    qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                {
                    badQualityCount++;
                }
            }
            if (badQualityCount == 0)
            {
                bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer], workerIndex);
            }
            for (var i = ConfigurationManager.Arguments.Kmer; i < simpleRead.Length; i++)
            {
                if (readSpan[i - ConfigurationManager.Arguments.Kmer] == Consts.InvalidBase ||
                    qualitySpan[i - ConfigurationManager.Arguments.Kmer] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                {
                    badQualityCount--;
                }
                if (readSpan[i] == Consts.InvalidBase ||
                    qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                {
                    badQualityCount++;
                }
                if (badQualityCount == 0)
                {
                    bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer), workerIndex);
                }
            }
            badQualityCount = 0;
            readSpan = Util.ReverseComprement(readSpan);
            qualitySpan.Reverse();
            for (var i = 0; i < ConfigurationManager.Arguments.Kmer; i++)
            {
                if (readSpan[i] == Consts.InvalidBase ||
                    qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                {
                    badQualityCount++;
                }
            }
            if (badQualityCount == 0)
            {
                bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer], workerIndex);
            }
            for (var i = ConfigurationManager.Arguments.Kmer; i < simpleRead.Length; i++)
            {
                if (readSpan[i - ConfigurationManager.Arguments.Kmer] == Consts.InvalidBase ||
                    qualitySpan[i - ConfigurationManager.Arguments.Kmer] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                {
                    badQualityCount--;
                }
                if (readSpan[i] == Consts.InvalidBase ||
                    qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                {
                    badQualityCount++;
                }
                if (badQualityCount == 0)
                {
                    bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer), workerIndex);
                }
            }
        }
    }
}