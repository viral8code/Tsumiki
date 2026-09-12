using Tsumiki.Cores.Pipeline;

namespace Tsumiki
{
    /// <summary>
    /// プロセスの入口
    /// </summary>
    internal static class Program
    {
        #region 内部メソッド

        /// <summary>
        /// アプリケーションを起動する
        /// </summary>
        /// <param name="p_引数列">コマンドライン引数</param>
        /// <returns>プロセスの終了コード</returns>
        private static int Main(string[] p_引数列)
        {
            return AssemblyApplication.Get_終了コード(p_引数列);
        }

        #endregion
    }
}
