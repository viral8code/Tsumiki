using System.Runtime.InteropServices;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// Win32 の IO_COUNTERS
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    internal struct 入出力計数
    {
        #region 公開変数

        /// <summary>
        /// 読み込み操作の回数
        /// </summary>
        public ulong A_読込回数;

        /// <summary>
        /// 書き込み操作の回数
        /// </summary>
        public ulong A_書込回数;

        /// <summary>
        /// その他の操作の回数
        /// </summary>
        public ulong A_その他回数;

        /// <summary>
        /// 読み込んだバイト数
        /// </summary>
        public ulong A_読込量;

        /// <summary>
        /// 書き込んだバイト数
        /// </summary>
        public ulong A_書込量;

        /// <summary>
        /// その他の操作で転送したバイト数
        /// </summary>
        public ulong A_その他量;

        #endregion
    }
}
