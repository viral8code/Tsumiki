using System.Runtime.CompilerServices;
using Tsumiki.Model;

namespace Tsumiki.Common
{
    internal class Logger
    {
        public static string Get_メソッド名([CallerMemberName] string p_メソッド名 = "")
        {
            return p_メソッド名;
        }

        /// <summary>標準出力へ1行出す。文言は言語ごとのカタログから引く。</summary>
        public static void V_出力(メッセージID p_ID, params object?[] p_引数)
        {
            Console.WriteLine(Messages.Get_文言(p_ID, p_引数));
        }

        /// <summary>標準エラーへ1行出す。</summary>
        public static void V_出力_標準エラー(メッセージID p_ID, params object?[] p_引数)
        {
            Console.Error.WriteLine(Messages.Get_文言(p_ID, p_引数));
        }

        public static void V_出力_警告(string p_メソッド名, Exception p_例外)
        {
            Logger.V_出力_標準エラー(メッセージID.例外を無視_見出し);
            Logger.V_出力_標準エラー(メッセージID.例外を無視_メソッド, p_メソッド名);

            // 例外の内容そのものは訳す対象ではない。
            Console.Error.WriteLine(p_例外.ToString());
        }

        public static void V_出力_エラー(string p_メソッド名, Exception p_例外)
        {
            Logger.V_出力_標準エラー(メッセージID.停止_見出し);
            Logger.V_出力_標準エラー(メッセージID.停止_メソッド, p_メソッド名);
            Console.Error.WriteLine(p_例外.ToString());
        }

        public static void V_出力_タイムスタンプ()
        {
            Logger.V_出力(メッセージID.タイムスタンプ, DateTime.Now);
        }

        /// <summary>区切りの空行。文言を持たないのでカタログには載せない。</summary>
        public static void V_出力_空行()
        {
            Console.WriteLine();
        }
    }
}
