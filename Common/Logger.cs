using System.Runtime.CompilerServices;

namespace Tsumiki.Common
{
    internal class Logger
    {
        public static string GetMethodName([CallerMemberName] string methodName = "")
        {
            return methodName;
        }

        public static void PrintWarning(string methodName, Exception ex)
        {
            Console.Error.WriteLine($"[Warning] The following exception will be ignored.");
            Console.Error.WriteLine($"Method: {methodName}");
            Console.Error.WriteLine(ex.ToString());
        }

        public static void PrintError(string methodName, Exception ex)
        {
            Console.Error.WriteLine($"[Error] Program was stopped.");
            Console.Error.WriteLine($"Method: {methodName}");
            Console.Error.WriteLine(ex.ToString());
        }

        public static void PrintTimeStamp()
        {
            Console.WriteLine($"[Log] {DateTime.Now.ToString("yyyy/MM/dd HH:mm:ss")}");
        }
    }
}
