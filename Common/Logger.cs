using System.Runtime.CompilerServices;

namespace Tsumiki.Common
{
    internal class Logger
    {
        public static string Get_メソッド名([CallerMemberName] string p_メソッド名 = "")
        {
            return p_メソッド名;
        }

        public static void V_出力_警告(string p_メソッド名, Exception p_例外)
        {
            Console.Error.WriteLine($"[Warning] The following exception will be ignored.");
            Console.Error.WriteLine($"          Method: {p_メソッド名}");
            Console.Error.WriteLine(p_例外.ToString());
        }

        public static void V_出力_エラー(string p_メソッド名, Exception p_例外)
        {
            Console.Error.WriteLine($"[Error] Program was stopped.");
            Console.Error.WriteLine($"        Method: {p_メソッド名}");
            Console.Error.WriteLine(p_例外.ToString());
        }

        public static void V_出力_タイムスタンプ()
        {
            Console.WriteLine($"[Log] {DateTime.Now:yyyy/MM/dd HH:mm:ss}");
        }
    }
}
