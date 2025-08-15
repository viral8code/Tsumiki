namespace Tsumiki.Common
{
    internal class Consts
    {
        public const string VERSION = "1.0";

        public static readonly List<string> AUTHOR_LIST = [
            "viral",
            ];

        public static readonly string DETAILS_TEXT = $"""
            Tsumiki is a genome assebly.
            author: {string.Join(", ", AUTHOR_LIST)}
            version: {VERSION}
            """;

        public static class ArgumentKey
        {
            public const string READPATH1 = "-1";

            public const string READPATH2 = "-2";

            public const string KMER = "-k";

            public const string KMERCUTOFF = "-kc";

            public const string PHRED = "-p";

            public const string QUARITYCUTOFF = "-q";
        }
    }
}
