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

            public const string KmerCutOff = "-kc";

            public const string Phred = "-p";

            public const string QualityCutOff = "-q";

            public const string BloomFilterSize = "-b";

            public const string AllowAmbiguousBases = "-ab";

            public const string Help = "-h";
        }

        public static readonly string HelpText = $"""
            {DetailsText}

            # Arguments
            {ArgumentKey.ReadPath1} : forward fastq(.gz) path
            {ArgumentKey.ReadPath2} : backward fastq(.gz) path
            {ArgumentKey.Kmer} : length of k-mer
            {ArgumentKey.KmerCutOff} : threshold of k-mer count (use kmers with this value or higher)
            {ArgumentKey.Phred} : base of phred score (33 or 64)
            {ArgumentKey.QualityCutOff} : threshold of base quality (use kmers with this value or higher)
            {ArgumentKey.BloomFilterSize} : memory allocation for the Bloom Filter
            {ArgumentKey.AllowAmbiguousBases} : allow ambiguous bases as valid bases (default : false)
            {ArgumentKey.Help} : output this text
            """;

        public const string UnitigFileName = "unitigs.fasta";
    }
}
