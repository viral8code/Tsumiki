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
                Console.WriteLine(Consts.DETAILS_TEXT);
                Environment.Exit(0);
            }

            Logger.PrintTimeStamp();

            var param = ArgumentsReader.ReadArguments(args);
            ConfigurationManager.Arguments = param;

            Console.WriteLine(param);

            Console.WriteLine("Start construction Bloom filter");
            var bloomFilter = new CountingBloomFilter(param.RowBitSize);
            ConfigurationManager.BloomFilter = bloomFilter;

            Console.WriteLine("Loading Read1");
            LoadReadFileToBloomFilter(param.ReadPath1, bloomFilter);

            Console.WriteLine("Loading Read2");
            LoadReadFileToBloomFilter(param.ReadPath2, bloomFilter);

            Console.WriteLine("Fix Bloom filter");
            var initKmers = bloomFilter.Cutoff(param.KmerCutoff);

            Logger.PrintTimeStamp();

            Console.WriteLine("Make unitigs");
            var unitigMaker = new UnitigMaker(bloomFilter);
            using (var writer = new FastaWriter("unitigs.fasta"))
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

        private static void LoadReadFileToBloomFilter(string fileName, CountingBloomFilter bloomFilter)
        {
            const int ProgressLogInterval = 100000;
            using var reader = new FastqReader(fileName);
            ulong count = 0;
            while (reader.HasNext())
            {
                var readData = reader.NextRead();
                if (readData.Read.Count < ConfigurationManager.Arguments.Kmer)
                {
                    continue;
                }
                var badQualityCount = 0;
                for (var i = 0; i < ConfigurationManager.Arguments.Kmer; i++)
                {
                    if (readData.Quality[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
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
                    if (readData.Quality[i - ConfigurationManager.Arguments.Kmer] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount--;
                    }
                    if (readData.Quality[i] - ConfigurationManager.Arguments.Phred - ConfigurationManager.Arguments.QualityCutoff < 0)
                    {
                        badQualityCount++;
                    }
                    if (badQualityCount == 0)
                    {
                        bloomFilter.Add(readSpan.Slice(i - ConfigurationManager.Arguments.Kmer + 1, ConfigurationManager.Arguments.Kmer));
                    }
                }
                if (++count % ProgressLogInterval == 0)
                {
                    Console.WriteLine(count + " reads Loaded");
                }
            }
        }
    }
}
