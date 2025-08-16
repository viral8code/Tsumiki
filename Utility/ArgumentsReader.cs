using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Utility
{
    internal class ArgumentsReader
    {

        private static class ArgumentKey
        {
            public const string READPATH1 = "-1";

            public const string READPATH2 = "-2";

            public const string KMER = "-k";

            public const string KMERCUTOFF = "-kc";

            public const string PHRED = "-p";

            public const string QUARITYCUTOFF = "-q";

            public const string BITSIZE = "-b";
        }

        public static Parameters ReadArguments(string[] args)
        {
            var param = new Parameters();
            try
            {
                var index = 0;
                while (index < args.Length)
                {
                    var key = args[index++];
                    switch (key)
                    {
                        case ArgumentKey.READPATH1:
                            param.ReadPath1 = args[index++];
                            break;

                        case ArgumentKey.READPATH2:
                            param.ReadPath2 = args[index++];
                            break;

                        case ArgumentKey.KMER:
                            param.Kmer = int.Parse(args[index++]);
                            break;

                        case ArgumentKey.KMERCUTOFF:
                            param.KmerCutoff = int.Parse(args[index++]);
                            break;

                        case ArgumentKey.PHRED:
                            param.Phred = int.Parse(args[index++]);
                            break;

                        case ArgumentKey.QUARITYCUTOFF:
                            param.QualityCutoff = int.Parse(args[index++]);
                            break;

                        case ArgumentKey.BITSIZE:
                            param.BitSize = args[index++];
                            break;

                        default:
                            Logger.PrintWarning(Logger.GetMethodName(), new ArgumentException($"Unknown argment: {key}"));
                            break;
                    }
                }
            }
            catch (Exception e)
            {
                Logger.PrintError(Logger.GetMethodName(), e);
                Environment.Exit(1);
            }

            if (string.IsNullOrWhiteSpace(param.ReadPath1) || string.IsNullOrWhiteSpace(param.ReadPath2))
            {
                Logger.PrintError(Logger.GetMethodName(), new ArgumentException("Please set pair end read's pathes"));
                Environment.Exit(1);
            }

            return param;
        }
    }
}
