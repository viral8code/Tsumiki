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

            PhredSniffer.ResolveOffset(param, param.ReadPath1, param.ReadPath2);

            Logger.PrintTimeStamp();

            var tempDir = Path.Combine(Environment.CurrentDirectory, param.TempDirectory);

            if (Path.Exists(tempDir))
            {
                Console.WriteLine($"{param.TempDirectory} already exists");
                Console.WriteLine("Please check path!");
                Environment.Exit(0);
            }

            _ = Directory.CreateDirectory(tempDir);

            if (param.EnableErrorCorrection)
            {
                Console.WriteLine("Correcting reads before assembly");

                var correctedPath1 = Path.Combine(tempDir, "corrected.1.fq");
                var hasRead2 = !string.IsNullOrWhiteSpace(param.ReadPath2);
                var correctedPath2 = hasRead2 ? Path.Combine(tempDir, "corrected.2.fq") : null;

                ErrorCorrector.CorrectReadFiles(
                    param.ReadPath1,
                    hasRead2 ? param.ReadPath2 : null,
                    tempDir,
                    correctedPath1,
                    correctedPath2);

                // 以降の全処理(k-merカウント・グラフ構築・リードの再マッピング)は
                // 訂正済みファイルを見るようにする。
                param.ReadPath1 = correctedPath1;
                if (hasRead2)
                {
                    param.ReadPath2 = correctedPath2!;
                }

                Logger.PrintTimeStamp();
            }

            if (param.RowBitSize != int.MaxValue)
            {
                Console.WriteLine($"[Info] {Consts.ArgumentKey.BloomFilterSize} is deprecated and no longer has any effect " +
                    "(k-mer membership is now tracked with an exact set built from the trusted k-mer count, not a Bloom filter).");
            }

            Console.WriteLine("Start construction k-mer index");
            var bloomFilter = new TrustedKmerIndex(tempDir);
            ConfigurationManager.TrustedKmerIndex = bloomFilter;

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
                KmerCounting.LoadReadFile(param.ReadPath1, bloomFilter);
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
                    KmerCounting.LoadReadFile(param.ReadPath2, bloomFilter);
                }
            }

            Logger.PrintTimeStamp();

            Console.WriteLine("Applying k-mer cutoff");
            _ = bloomFilter.Cutoff(param.KmerCutoff);

            Logger.PrintTimeStamp();

            Console.WriteLine("Clipping short tips");
            var initKmers = GraphSimplifier.ClipTips(bloomFilter, param.Kmer);

            Logger.PrintTimeStamp();

            Console.WriteLine("Make unitigs");

            var unitigMaker = new UnitigMaker(bloomFilter);
            HashSet<string> unitigSet = [];
            Dictionary<int, string> unitigSequences = [];
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
                    unitigSequences[id] = unitig.Sequence;
                    writer.Write(id++, unitig.Sequence);
                    if (id > Consts.MaximumUnitigCount)
                    {
                        break;
                    }
                }
            }

            AssemblyStatsReporter.Report("unitigs", Consts.UnitigFileName);

            // 各 unitig のカバレッジからコピー数を推定する。反復配列かどうかを
            // グラフの形ではなく量的な根拠で判定でき、後段の経路探索では
            // 「この unitig は何回まで使ってよいか」という予算になる。
            // k-mer インデックスがまだ生きているこの時点でしか計算できない。
            var unitigLengthMap = unitigSequences.ToDictionary(kv => kv.Key, kv => kv.Value.Length);
            var unitigCoverage = CopyNumberEstimator.ComputeCoverage(bloomFilter, unitigSequences, param.Kmer);
            var copyNumbers = CopyNumberEstimator.Estimate(unitigCoverage, unitigLengthMap);
            CopyNumberEstimator.Report(copyNumbers, unitigLengthMap);

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

            contigMaker.UniteContigs(Consts.ContigFileName, param.PairUniteThreshold, param.PairCountThreshold, copyNumbers.CopyNumber);

            Console.WriteLine("Maked contigs");

            AssemblyStatsReporter.Report("contigs", Consts.ContigFileName);

            Logger.PrintTimeStamp();

            // スキャフォールディングはペアエンド情報(pairPath)を前提とするため、
            // read2 が指定されている(=ペアエンドで実行された)場合のみ行う。
            if (!string.IsNullOrWhiteSpace(param.ReadPath2))
            {
                Console.WriteLine("Scaffolding contigs");

                var scaffolder = new Scaffolder(contigMaker, Consts.ContigFileName);
                scaffolder.Run(Consts.ScaffoldFileName);

                AssemblyStatsReporter.Report("scaffolds", Consts.ScaffoldFileName);

                // スキャフォールドの N を、グラフ上で両端を繋ぐ経路を探して
                // 実配列に置き換える。contig が途切れたのは配列が無いからではなく
                // 分岐で決められなかったからであることが多く、その場合
                // ギャップを埋める配列はグラフ上に実在する。
                Console.WriteLine("Filling scaffold gaps");
                var gapStats = GapFiller.Run(Consts.ScaffoldFileName, bloomFilter, param.Kmer);
                GapFiller.Report(gapStats);
                if (gapStats.FilledGaps > 0)
                {
                    AssemblyStatsReporter.Report("scaffolds (gaps filled)", Consts.ScaffoldFileName);
                }

                Logger.PrintTimeStamp();
            }

            Console.WriteLine("開発中！");

            Logger.PrintTimeStamp();
        }

        // 曖昧塩基を許容する経路は現状シングルスレッドのまま(workerIndex固定)。
        // 呼ばれる頻度が低い想定のため、並列化の優先度を下げている。
        // (Core.KmerCounting.LoadReadFile は既定の「曖昧塩基を無視する」経路のみ
        //  切り出したもので、こちらは対象外。)
        private static void LoadReadFileToBloomFilterWithAmbiguity(string filePath, TrustedKmerIndex bloomFilter)
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
    }
}