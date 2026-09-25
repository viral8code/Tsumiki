using System.Runtime.InteropServices;

namespace Tsumiki.Utilities
{
    /// <summary>Win32 の IO_COUNTERS</summary>
    [StructLayout(LayoutKind.Sequential)]
    internal struct 入出力計数
    {
        public ulong A_読込回数;
        public ulong A_書込回数;
        public ulong A_その他回数;
        public ulong A_読込量;
        public ulong A_書込量;
        public ulong A_その他量;
    }
}
