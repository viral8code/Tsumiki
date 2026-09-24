using System.Diagnostics;
using System.Runtime.InteropServices;
using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Models.Reporting;

namespace Tsumiki.Utilities
{
    /// <summary>工程の経過時間とプロセス全体の資源使用量を記録する</summary>
    internal sealed class StageTimer : IDisposable
    {
        #region 内部変数

        /// <summary>工程名</summary>
        private readonly string _工程;
        /// <summary>経過時間</summary>
        private readonly Stopwatch _時計 = Stopwatch.StartNew();
        /// <summary>計測するプロセス</summary>
        private readonly Process _プロセス = Process.GetCurrentProcess();
        /// <summary>開始時の CPU 時間</summary>
        private readonly TimeSpan _CPU;
        /// <summary>開始時の累積確保量</summary>
        private readonly long _確保量 = GC.GetTotalAllocatedBytes(false);
        /// <summary>開始時の世代 2 回収回数</summary>
        private readonly int _回収数 = GC.CollectionCount(2);
        /// <summary>開始時のプロセス全体の読み書き量</summary>
        private readonly (ulong A_読込, ulong A_書込)? _入出力 = Get_入出力量();

        #endregion

        #region コンストラクタ

        /// <summary>工程の計測を始める</summary>
        /// <param name="p_工程">配列や入力パスを含まない工程名</param>
        public StageTimer(string p_工程)
        {
            this._工程 = p_工程;
            this._CPU = this._プロセス.TotalProcessorTime;
            Logger.V_出力_そのまま($"[Perf] {p_工程}: start");
        }

        #endregion

        #region 公開メソッド

        /// <summary>計測区間の終了時点の資源使用量を記録する</summary>
        public void Dispose()
        {
            this._プロセス.Refresh();
            var l_経過秒 = this._時計.Elapsed.TotalSeconds;
            var l_CPU秒 = (this._プロセス.TotalProcessorTime - this._CPU).TotalSeconds;
            var l_確保MB = (GC.GetTotalAllocatedBytes(false) - this._確保量) / 1048576D;
            var l_ワーキングセットMB = this._プロセス.WorkingSet64 / 1048576D;
            var l_ピークワーキングセットMB = this._プロセス.PeakWorkingSet64 / 1048576D;
            var l_世代2回収回数 = GC.CollectionCount(2) - this._回収数;
            var l_入出力 = this._入出力 is { } l_開始 && Get_入出力量() is { } l_終了
                ? FormattableString.Invariant($", io_read_mb={(l_終了.A_読込 - l_開始.A_読込) / 1048576D:F1}, io_write_mb={(l_終了.A_書込 - l_開始.A_書込) / 1048576D:F1}")
                : string.Empty;
            Logger.V_出力_そのまま(FormattableString.Invariant($"[Perf] {this._工程}: elapsed_s={l_経過秒:F3}, cpu_s={l_CPU秒:F3}, allocated_mb={l_確保MB:F1}, working_set_mb={l_ワーキングセットMB:F1}, peak_working_set_mb={l_ピークワーキングセットMB:F1}, gen2_gc={l_世代2回収回数}{l_入出力}"));
            PhaseTimingRecorder.V_記録(new フェーズ計測(this._工程, l_経過秒, l_CPU秒, l_確保MB, l_ワーキングセットMB, l_ピークワーキングセットMB, l_世代2回収回数));
            this._プロセス.Dispose();
        }

        #endregion

        #region 内部メソッド

        /// <summary>プロセスがこれまでに読み書きした量</summary>
        /// <remarks>
        /// OS ごとに取り方が違い、macOS には手軽な口が無いので取らない
        /// </remarks>
        /// <returns>取れなければ null</returns>
        private static (ulong A_読込, ulong A_書込)? Get_入出力量()
        {
            if (OperatingSystem.IsWindows())
            {
                using var l_プロセス = Process.GetCurrentProcess();
                return GetProcessIoCounters(l_プロセス.Handle, out var l_計数) ? (l_計数.A_読込量, l_計数.A_書込量) : null;
            }
            return OperatingSystem.IsLinux() ? Get_入出力量_Linux() : null;
        }

        /// <summary>/proc/self/io から読み書きした量を取る</summary>
        /// <remarks>
        /// Windows の値と揃えるため、ディスクまで届いた量 (read_bytes) ではなく読み書きの呼び出しで渡した量 (rchar・wchar) を使う
        /// </remarks>
        /// <returns>取れなければ null</returns>
        private static (ulong A_読込, ulong A_書込)? Get_入出力量_Linux()
        {
            try
            {
                ulong? l_読込 = null;
                ulong? l_書込 = null;
                foreach (var l_行 in File.ReadLines("/proc/self/io"))
                {
                    var l_区切り = l_行.IndexOf(':');
                    if (l_区切り < 0 || !ulong.TryParse(l_行.AsSpan(l_区切り + 1).Trim(), out var l_値))
                    {
                        continue;
                    }
                    switch (l_行[..l_区切り])
                    {
                        case "rchar":
                            l_読込 = l_値;
                            break;
                        case "wchar":
                            l_書込 = l_値;
                            break;
                    }
                }
                return l_読込 is { } l_r && l_書込 is { } l_w ? (l_r, l_w) : null;
            }
            catch (IOException)
            {
                return null;
            }
            catch (UnauthorizedAccessException)
            {
                return null;
            }
        }

        /// <summary>Win32 の IO_COUNTERS</summary>
        [StructLayout(LayoutKind.Sequential)]
        private struct 入出力計数
        {
            public ulong A_読込回数;
            public ulong A_書込回数;
            public ulong A_その他回数;
            public ulong A_読込量;
            public ulong A_書込量;
            public ulong A_その他量;
        }

        /// <summary>プロセスの読み書き量を取る</summary>
        [DllImport("kernel32.dll")]
        [return: MarshalAs(UnmanagedType.Bool)]
        private static extern bool GetProcessIoCounters(IntPtr p_プロセス, out 入出力計数 p_計数);

        #endregion
    }
}
