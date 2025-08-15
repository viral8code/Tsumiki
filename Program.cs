using System.Reflection.Metadata;
using Tsumiki.Common;
using Tsumiki.Util;

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

            var param = ArgumentsReader.ReadArguments(args);
            Console.WriteLine("開発中！");
        }
    }
}
