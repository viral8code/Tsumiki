using System.Diagnostics;
using Tsumiki.Commons;

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
            Logger.V_出力_そのまま(FormattableString.Invariant($"[Perf] {this._工程}: elapsed_s={this._時計.Elapsed.TotalSeconds:F3}, cpu_s={(this._プロセス.TotalProcessorTime - this._CPU).TotalSeconds:F3}, allocated_mb={(GC.GetTotalAllocatedBytes(false) - this._確保量) / 1048576D:F1}, working_set_mb={this._プロセス.WorkingSet64 / 1048576D:F1}, gen2_gc={GC.CollectionCount(2) - this._回収数}"));
            this._プロセス.Dispose();
        }

        #endregion
    }
}
