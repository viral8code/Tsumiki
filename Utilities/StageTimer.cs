using System.Diagnostics;
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
            Logger.V_出力_そのまま(FormattableString.Invariant($"[Perf] {this._工程}: elapsed_s={l_経過秒:F3}, cpu_s={l_CPU秒:F3}, allocated_mb={l_確保MB:F1}, working_set_mb={l_ワーキングセットMB:F1}, peak_working_set_mb={l_ピークワーキングセットMB:F1}, gen2_gc={l_世代2回収回数}"));
            PhaseTimingRecorder.V_記録(new フェーズ計測(this._工程, l_経過秒, l_CPU秒, l_確保MB, l_ワーキングセットMB, l_ピークワーキングセットMB, l_世代2回収回数));
            this._プロセス.Dispose();
        }

        #endregion
    }
}
