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

                        case Consts.ArgumentKey.KmerCutoff:
                            param.KmerCutoff = ulong.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.Phred:
                            param.Phred = int.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.QualityCutoff:
                            param.QualityCutoff = int.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.MemoryBudget:
                            param.MemoryBudgetMB = int.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.BloomFilterSize:
                            param.BitSize = args[index++];
                            break;

                        case Consts.ArgumentKey.InsertSize:
                            param.InsertSize = int.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.TempDirectory:
                            param.TempDirectory = args[index++];
                            break;

                        case Consts.ArgumentKey.ThreadCount:
                            param.ThreadCount = int.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.PairUniteThreshold:
                            param.PairUniteThreshold = decimal.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.PairCountThreshold:
                            param.PairCountThreshold = ulong.Parse(args[index++]);
                            break;

                        case Consts.ArgumentKey.Help:
                            param.IsHelpMode = true;
                            break;

                        case Consts.ArgumentKey.AllowAmbiguousBases:
                            param.AllowAmbiguousBases = true;
                            break;

                        case Consts.ArgumentKey.ErrorCorrection:
                            param.EnableErrorCorrection = true;
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

            // ヘルプ表示だけを求められている場合は、リードパスの必須チェックを行わない。
            // (以前は -h のみを指定してもここで「Please set read path」エラーになり
            //  ヘルプが表示できなかった)
            if (param.IsHelpMode)
            {
                return param;
            }

            if (string.IsNullOrWhiteSpace(param.ReadPath1))
            {
                param.ReadPath1 = param.ReadPath2;
                param.ReadPath2 = string.Empty;
            }

            if (string.IsNullOrWhiteSpace(param.ReadPath1))
            {
                Logger.PrintError(Logger.GetMethodName(), new ArgumentException("Please set read path"));
                Environment.Exit(0);
            }

            return param;
        }
    }
}