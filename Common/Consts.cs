namespace Tsumiki.Common
{
    internal class Consts
    {
        public const string Version = "1.0";

        public static readonly List<string> AuthorList = [
            "viral",
            ];

        public static readonly string DetailsText = $"""

            Tsumiki is a genome assembler.
            author: {string.Join(", ", AuthorList)}
            version: {Version}

            """;

        public static class ArgumentKey
        {
            public const string ReadPath1 = "-1";

            public const string ReadPath2 = "-2";

            public const string Kmer = "-k";

            public const string KmerCutoff = "-kc";

            public const string Phred = "-p";

            public const string QualityCutoff = "-q";

            public const string BloomFilterSize = "-b";

            public const string InsertSize = "-i";

            public const string AllowAmbiguousBases = "-ab";

            public const string Help = "-h";

            public const string TempDirectory = "-t";

            public const string ThreadCount = "-th";

            public const string PairUniteThreshold = "-pu";

            public const string PairCountThreshold = "-pc";

            public const string ErrorCorrection = "-ec";
        }

        public const string NullInsertSizeText = "unspecified";

        public const int DefaultKmerValue = 31;

        public const int DefaultKmerCutoffValue = 2;

        public const int DefaultPhredValue = 33;

        public const int DefaultQualityCutoffValue = 1;

        public const ulong DefaultBloolFilterSize = int.MaxValue;

        public const string DefaultTempFolder = "temp";

        public static readonly int[] AllowedPhredValue = [33, 64];

        public const decimal DefaultPairUniteThreshold = 0.8m;

        public const ulong DefaultPairCountThreshold = 10;

        public static readonly string HelpText = $"""
            {DetailsText}

            # Arguments
            {ArgumentKey.ReadPath1} [path] : forward fastq(.gz) path (required) (when using single reads, set the path using this argument)
            {ArgumentKey.ReadPath2} [path] : backward fastq(.gz) path
            {ArgumentKey.Kmer} [integer] : length of k-mer (default : {DefaultKmerValue})
            {ArgumentKey.KmerCutoff} [integer] : threshold of k-mer count (use kmers with this value or higher) (default : {DefaultKmerCutoffValue})
            {ArgumentKey.Phred} [integer] : base of phred score ({string.Join(" or ", AllowedPhredValue)}) (default : {DefaultPhredValue})
            {ArgumentKey.QualityCutoff} [integer] : threshold of base quality (use kmers with this value or higher) (default : {DefaultQualityCutoffValue})
            {ArgumentKey.BloomFilterSize} [decimal] : memory allocation for the Bloom Filter (e.g. 300M, 1.2G) (default : 200M)
            {ArgumentKey.InsertSize} : excepted insert size of pair-end reads (default : {NullInsertSizeText}, auto-estimated from mapped pairs when possible)
            {ArgumentKey.TempDirectory} [path] : temp directory (default : {DefaultTempFolder})
            {ArgumentKey.ThreadCount} [integer] : number of worker threads used for loading reads (default : number of logical processors)
            {ArgumentKey.PairUniteThreshold} [decimal] : minimum ratio of the best-supported pair-end scaffold edge among all candidates for a node (default : {DefaultPairUniteThreshold})
            {ArgumentKey.PairCountThreshold} [integer] : minimum read-pair support required for a pair-end scaffold edge (default : {DefaultPairCountThreshold})
            {ArgumentKey.ErrorCorrection} : run k-mer-spectrum-based read error correction before assembly (default : false)
            {ArgumentKey.Help} : output this text (default : false)

            """;

        public const string UnitigFileName = "unitigs.fasta";

        public const string ContigFileName = "contigs.fasta";

        public const string ScaffoldFileName = "scaffolds.fasta";

        public static class NucleotideID
        {
            public const byte A = 1;
            public const byte C = 2;
            public const byte G = 3;
            public const byte T = 4;
        }

        public const ulong ProgressLogInterval = 100_000;

        public const byte InvalidBase = 5;

        public const int MaximumUnitigCount = 100_000;

        /// <summary>
        /// InsertSize が未指定の場合の自動推定に使う、単一unitigへ両リードが
        /// マップされたペアの最小サンプル数。これに満たない場合は推定を諦め、
        /// ペアエンド由来のスキャフォールディングをスキップする。
        /// </summary>
        public const int MinInsertSizeSampleCount = 30;

        /// <summary>
        /// ギャップ長が推定上0以下になった場合に最低限挿入するNの数。
        /// </summary>
        public const int MinimumGapLength = 1;
    }
}