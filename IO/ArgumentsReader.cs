using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.IO
{
    internal class ArgumentsReader
    {
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
                        case Consts.ArgumentKey.ReadPath1:
                            param.ReadPath1 = args[index++];
                            break;

                        case Consts.ArgumentKey.ReadPath2:
                            param.ReadPath2 = args[index++];
                            break;

                        case Consts.ArgumentKey.Kmer:
                            param.Kmer = int.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.KmerCutOff:
                            param.KmerCutoff = ulong.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.Phred:
                            param.Phred = int.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.QualityCutOff:
                            param.QualityCutoff = int.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.BloomFilterSize:
                            param.BitSize = args[index++];
                            break;

                        case Consts.ArgumentKey.Help:
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
