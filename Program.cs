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

            Console.WriteLine("開発中！");

            Logger.PrintTimeStamp();
        }
    }
}
