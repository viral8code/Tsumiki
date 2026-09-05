using System.Collections.Concurrent;
using System.Runtime.ExceptionServices;

namespace Tsumiki.Utility
{
    /// <summary>
    /// リードを1本のスレッドで順に読み進めつつ、ワーカー群へ配って並列に処理する。
    ///
    /// 素朴に書くと、ワーカーが例外で落ちたときキューを引き取る者がいなくなり、
    /// プロデューサーが満杯のキューへの Add で永久に待つ。Task.WaitAll に
    /// 到達しないため例外も観測されず、無言のハングになる。それを防ぐ。
    /// </summary>
    internal static class ReadPipeline
    {
        /// <summary>
        /// p_供給元 の各要素を p_処理 へ並列に配る。
        /// p_処理 の第2引数はワーカー番号で、ワーカーごとのローカル集計用配列の
        /// 添字として使うことを想定している。
        /// </summary>
        public static void V_実行<T>(
            int p_スレッド数, int p_キュー容量, IEnumerable<T> p_供給元, Action<T, int> p_処理)
        {
            using var l_キュー = new BlockingCollection<T>(p_キュー容量);
            using var l_中断 = new CancellationTokenSource();

            var l_ワーカー = new Task[p_スレッド数];
            for (var w = 0; w < p_スレッド数; w++)
            {
                var l_ワーカー番号 = w;
                l_ワーカー[w] = Task.Run(() =>
                {
                    try
                    {
                        foreach (var l_項目 in l_キュー.GetConsumingEnumerable())
                        {
                            p_処理(l_項目, l_ワーカー番号);
                        }
                    }
                    catch
                    {
                        // 伝えないと、供給側が満杯のキューへの Add で永久に待つ。
                        l_中断.Cancel();
                        throw;
                    }
                });
            }

            ExceptionDispatchInfo? l_供給側の例外 = null;
            try
            {
                foreach (var l_項目 in p_供給元)
                {
                    l_キュー.Add(l_項目, l_中断.Token);
                }
            }
            catch (OperationCanceledException)
            {
                // 真の原因はワーカー側の例外で、下の WaitAll が送出する。
            }
            catch (Exception ex)
            {
                l_供給側の例外 = ExceptionDispatchInfo.Capture(ex);
            }

            // 供給が途中で終わっても、待っているワーカーを必ず解放する。
            l_キュー.CompleteAdding();

            // 供給側の例外よりワーカー側を優先する。供給側の中断は
            // ワーカーが落ちた結果であることが多いため。
            Task.WaitAll(l_ワーカー);

            l_供給側の例外?.Throw();
        }
    }
}
