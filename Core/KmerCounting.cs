using System.Collections.Concurrent;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// FASTQ ファイルを読み進めて TrustedKmerIndex へ k-mer を登録する処理
    /// (曖昧塩基を無視する既定経路)。元は Program.cs の
    /// LoadReadFileToBloomFilterIgnoreAmbiguity / ProcessRead だったものを、
    /// 本パイプラインと ErrorCorrector の事前カウントパスの両方から
    /// 呼べるよう切り出したもの。ロジック自体は変更していない。
    /// </summary>
    internal static class KmerCounting
    {
        /// <summary>
        /// FASTQ を1本のスレッドで順に読み進めつつ(ディスクI/Oはシーケンシャルなまま)、
        /// 読み取ったリードを BlockingCollection 経由でワーカースレッド群に配る
        /// プロデューサー/コンシューマ方式。各ワーカーは自分専用の workerIndex を使って
        /// TrustedKmerIndex.Add を呼ぶため、ロックなしで並列に k-mer を登録できる。
        /// </summary>
        public static void LoadReadFile(string filePath, TrustedKmerIndex index)
        {
            var threadCount = Math.Max(1, ConfigurationManager.Arguments.ThreadCount);
            ulong totalCount = 0;
            ulong mult = 0;
            var countLock = new object();

            using var queue = new BlockingCollection<ReadData>(boundedCapacity: threadCount * 64);

            var workers = new Task[threadCount];
            for (var w = 0; w < threadCount; w++)
            {
                var workerIndex = w;
                workers[w] = Task.Run(() =>
                {
                    foreach (var readData in queue.GetConsumingEnumerable())
                    {
                        ProcessRead(readData, index, workerIndex);

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
        /// 1リード分の k-mer 抽出・品質フィルタリング・TrustedKmerIndex 登録を行う。
        /// </summary>
        private static void ProcessRead(ReadData readData, TrustedKmerIndex index, int workerIndex)
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
                index.Add(readSpan[..ConfigurationManager.Arguments.Kmer], workerIndex);
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
                    index.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer), workerIndex);
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
                index.Add(readSpan[..ConfigurationManager.Arguments.Kmer], workerIndex);
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
                    index.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer), workerIndex);
                }
            }
        }
    }
}
