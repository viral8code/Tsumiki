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
            try
            {
                Run(args);
            }
            catch (Exception ex)
            {
                Logger.PrintError("Unhandled Tsumiki's method", ex);
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
                    writer.Write(unitig.Id, unitig.Sequence);
                }
            }

            Logger.PrintTimeStamp();

            Console.WriteLine("開発中！");

            Directory.Delete(tempDir);

            Logger.PrintTimeStamp();
        }

        private static void LoadReadFileToBloomFilterWithAmbiguity(string filePath, CountingBloomFilter bloomFilter)
        {
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
                var readSpan = CollectionsMarshal.AsSpan(readData.Read);
                for (var i = 0; i < readData.Quality.Length; i++)
                {
                    if (readData.Quality[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        readSpan[i] = [Consts.NucleotideID.A, Consts.NucleotideID.C, Consts.NucleotideID.G, Consts.NucleotideID.T];
                    }
                }
                bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer]);
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
                {
                    bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
                }
                readSpan = Util.ReverseComprement(readSpan);
                bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer]);
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
                {
                    bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
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

        private static void LoadReadFileToBloomFilterIgnoreAmbiguity(string filePath, CountingBloomFilter bloomFilter)
        {
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
                var readSpan = CollectionsMarshal.AsSpan(readData.Read.Select(x => x.Length > 1 ? Consts.InvalidBase : x[0]).ToList());
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
                    bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer]);
                }
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
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
                        bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
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
                    bloomFilter.Add(readSpan[..ConfigurationManager.Arguments.Kmer]);
                }
                for (var i = ConfigurationManager.Arguments.Kmer; i < readData.Read.Count; i++)
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
                        bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
                    }
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
    }
}
