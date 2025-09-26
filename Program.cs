using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Utility;

namespace Tsumiki
{
    internal class Program
    {
        private static void Main(string[] args)
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

            Console.WriteLine("Start construction Bloom filter");
            var bloomFilter = new CountingBloomFilter(param.RowBitSize);
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
            using (var writer = new FastaWriter(Consts.UnitigFileName))
            {
                foreach (var kmer in initKmers)
                {
                    var unitig = unitigMaker.MakeUnitig(kmer);
                    writer.Write(unitig.Id, unitig.Sequence);
                }
            }

            Logger.PrintTimeStamp();

            Console.WriteLine("開発中！");

            Logger.PrintTimeStamp();
        }

        private static void LoadReadFileToBloomFilterWithAmbiguity(string filePath, CountingBloomFilter bloomFilter)
        {
            const ulong ProgressLogInterval = 100000;

            ulong count = 0;
            ulong mult = 0;

            using var reader = new FastqReader(filePath);
            while (reader.HasNext())
            {
                var readData = reader.NextRead();
                if (readData.Read.Count < ConfigurationManager.Arguments.Kmer)
                {
                    continue;
                }
                var badQualityCount = 0;
                var qualitySpan = readData.Quality.ToCharArray().AsSpan();
                for (var i = 0; i < ConfigurationManager.Arguments.Kmer; i++)
                {
                    if (qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                }
                var readSpan = CollectionsMarshal.AsSpan(readData.Read);
                if (badQualityCount == 0)
                {
                    bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer]);
                }
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
                {
                    if (qualitySpan[i - ConfigurationManager.Arguments.Kmer] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount--;
                    }
                    if (qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                    if (badQualityCount == 0)
                    {
                        bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
                    }
                }
                badQualityCount = 0;
                readSpan = Util.ReverseComprement(readSpan);
                qualitySpan.Reverse();
                for (var i = 0; i < ConfigurationManager.Arguments.Kmer; i++)
                {
                    if (qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                }
                if (badQualityCount == 0)
                {
                    bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer]);
                }
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
                {
                    if (qualitySpan[i - ConfigurationManager.Arguments.Kmer] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount--;
                    }
                    if (qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                    if (badQualityCount == 0)
                    {
                        bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
                    }
                }
                if (++count == ProgressLogInterval)
                {
                    Console.WriteLine((++mult * ProgressLogInterval) + " reads Loaded");
                    count = 0;
                }
            }
            var fileName = Path.GetFileName(filePath);
            Console.WriteLine($"Loaded {(mult * ProgressLogInterval) + count} reads from {fileName}");
        }

        private static void LoadReadFileToBloomFilterIgnoreAmbiguity(string filePath, CountingBloomFilter bloomFilter)
        {
            const ulong ProgressLogInterval = 100000;
            const byte InvalidBase = 5;

            ulong count = 0;
            ulong mult = 0;

            using var reader = new FastqReader(filePath);
            while (reader.HasNext())
            {
                var readData = reader.NextRead();
                if (readData.Read.Count < ConfigurationManager.Arguments.Kmer)
                {
                    continue;
                }
                var badQualityCount = 0;
                var qualitySpan = readData.Quality.ToCharArray().AsSpan();
                var readSpan = CollectionsMarshal.AsSpan(readData.Read.Select(x => x.Length > 1 ? InvalidBase : x[0]).ToList());
                for (var i = 0; i < ConfigurationManager.Arguments.Kmer; i++)
                {
                    if (readSpan[i] == InvalidBase ||
                        qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                }
                if (badQualityCount == 0)
                {
                    bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer]);
                }
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
                {
                    if (readSpan[i - ConfigurationManager.Arguments.Kmer] == InvalidBase ||
                        qualitySpan[i - ConfigurationManager.Arguments.Kmer] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount--;
                    }
                    if (readSpan[i] == InvalidBase ||
                        qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                    if (badQualityCount == 0)
                    {
                        bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
                    }
                }
                badQualityCount = 0;
                readSpan = Util.ReverseComprement(readSpan);
                qualitySpan.Reverse();
                for (var i = 0; i < ConfigurationManager.Arguments.Kmer; i++)
                {
                    if (readSpan[i] == InvalidBase ||
                        qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                }
                if (badQualityCount == 0)
                {
                    bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer]);
                }
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
                {
                    if (readSpan[i - ConfigurationManager.Arguments.Kmer] == InvalidBase ||
                        qualitySpan[i - ConfigurationManager.Arguments.Kmer] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount--;
                    }
                    if (readSpan[i] == InvalidBase ||
                        qualitySpan[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                    if (badQualityCount == 0)
                    {
                        bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
                    }
                }
                if (++count == ProgressLogInterval)
                {
                    Console.WriteLine((++mult * ProgressLogInterval) + " reads Loaded");
                    count = 0;
                }
            }
            var fileName = Path.GetFileName(filePath);
            Console.WriteLine($"Loaded {(mult * ProgressLogInterval) + count} reads from {fileName}");
        }
    }
}
