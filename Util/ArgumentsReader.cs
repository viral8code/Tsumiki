using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Tsumiki.Common;
using Tsumiki.Model;
using static Tsumiki.Common.Consts;

namespace Tsumiki.Util
{
    internal class ArgumentsReader
    {
        public static Parameters ReadArguments(string[] args)
        {
            var index = 0;
            var param = new Parameters();
            try
            {
                while (index < args.Length)
                {
                    var key = args[index++];
                    switch (key)
                    {
                        case ArgumentKey.READPATH1:
                            var path1 = args[index++];
                            param.ReadPath1 = path1;
                            break;
                        case ArgumentKey.READPATH2:
                            var path2 = args[index++];
                            param.ReadPath2 = path2;
                            break;
                        case ArgumentKey.KMER:
                            var kmer = args[index++];
                            param.Kmer = int.Parse(kmer);
                            break;
                        case ArgumentKey.KMERCUTOFF:
                            var kmerCutoff = args[index++];
                            param.KmerCutoff = int.Parse(kmerCutoff);
                            break;
                        case ArgumentKey.PHRED:
                            var phred = args[index++];
                            param.Phred = int.Parse(phred);
                            break;
                        case ArgumentKey.QUARITYCUTOFF:
                            var qualityCutoff = args[index++];
                            param.QualityCutoff = int.Parse(qualityCutoff);
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
