using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace Tsumiki.Common
{
    internal class Logger
    {
        public static string GetMethodName([CallerMemberName] string methodName = "")
        {
            return methodName;
        }

        public static void PrintInnerException(string methodName, Exception ex)
        {
            Console.Error.WriteLine($"InnerException");
            Console.Error.WriteLine($"Method: {methodName}");
            Console.Error.WriteLine($"Message: {ex}");
        }
    }
}
