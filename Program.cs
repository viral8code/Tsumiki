using Tsumiki.Common;
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
            using (var reader = new FastqReader(param.ReadPath1))
            {
                while (reader.HasNext())
                {
                    var readData = reader.NextRead();
                    if (readData.Read.Length < param.Kmer)
                    {
                        continue;
                    }
                    var badQualityCount = 0;
                    for (var i = 0; i < param.Kmer; i++)
                    {
                        if (readData.Quality[i] - param.Phred - param.QualityCutoff < 0)
                        {
                            badQualityCount++;
                        }
                    }
                    if (badQualityCount == 0)
                    {
                        bloomFilter.Add(readData.Read[..param.Kmer]);
                    }
                    for (var i = param.Kmer; i < readData.Read.Length; i++)
                    {
                        if (readData.Quality[i - param.Kmer] - param.Phred - param.QualityCutoff < 0)
                        {
                            badQualityCount--;
                        }
                        if (readData.Quality[i] - param.Phred - param.QualityCutoff < 0)
                        {
                            badQualityCount++;
                        }
                        if (badQualityCount == 0)
                        {
                            bloomFilter.Add(readData.Read.Substring(i - param.Kmer + 1, param.Kmer));
                        }
                    }
                }
            }

            Console.WriteLine("Loading Read2");
            using (var reader = new FastqReader(param.ReadPath2))
            {
                while (reader.HasNext())
                {
                    var readData = reader.NextRead();
                    if (readData.Read.Length < param.Kmer)
                    {
                        continue;
                    }
                    var badQualityCount = 0;
                    for (var i = 0; i < param.Kmer; i++)
                    {
                        if (readData.Quality[i] - param.Phred - param.QualityCutoff < 0)
                        {
                            badQualityCount++;
                        }
                    }
                    if (badQualityCount == 0)
                    {
                        bloomFilter.Add(readData.Read[..param.Kmer]);
                    }
                    for (var i = param.Kmer; i < readData.Read.Length; i++)
                    {
                        if (readData.Quality[i - param.Kmer] - param.Phred - param.QualityCutoff < 0)
                        {
                            badQualityCount--;
                        }
                        if (readData.Quality[i] - param.Phred - param.QualityCutoff < 0)
                        {
                            badQualityCount++;
                        }
                        if (badQualityCount == 0)
                        {
                            bloomFilter.Add(readData.Read.Substring(i - param.Kmer + 1, param.Kmer));
                        }
                    }
                }
            }
            Console.WriteLine("Fix Bloom filter");
            bloomFilter.Cutoff(param.KmerCutoff);

            Console.WriteLine("開発中！");

            Logger.PrintTimeStamp();
        }
    }
}
