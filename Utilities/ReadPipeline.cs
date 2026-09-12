using System.Collections.Concurrent;
using System.Runtime.ExceptionServices;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// リードを 1 本のスレッドで順に読み進めつつ、ワーカー群へ配って並列に処理する
    /// </summary>
    /// <remarks>
    /// 素朴に書くと、ワーカーが例外で落ちたときキューを引き取る者がいなくなり、プロデューサーが満杯のキューへの Add で永久に待つ<br/>
    /// Task.WaitAll に到達しないため例外も観測されず、無言のハングになる<br/>
    /// それを防ぐ
    /// </remarks>
    internal static class ReadPipeline
    {
        #region 公開メソッド

        /// <summary>
        /// p_供給元 の各要素を p_処理 へ並列に配る
        /// </summary>
        /// <param name="p_スレッド数"></param>
        /// <param name="p_キュー容量"></param>
        /// <param name="p_供給元"></param>
        /// <param name="p_処理"></param>
        /// <remarks>
        /// p_処理 の第 2 引数はワーカー番号で、ワーカーごとのローカル集計用配列の添字として使うことを想定している
        /// </remarks>
        public static void V_実行<T>(int p_スレッド数, int p_キュー容量, IEnumerable<T> p_供給元, Action<T, int> p_処理)
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
                        // 伝えないと、供給側が満杯のキューへの Add で永久に待つ
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
                // 真の原因はワーカー側の例外で、下の WaitAll が送出する
            }
            catch (Exception l_例外)
            {
                l_供給側の例外 = ExceptionDispatchInfo.Capture(l_例外);
            }

            // 供給が途中で終わっても、待っているワーカーを必ず解放する
            l_キュー.CompleteAdding();

            // 供給側の例外よりワーカー側を優先する
            // 供給側の中断は
            // ワーカーが落ちた結果であることが多いため
            Task.WaitAll(l_ワーカー);

            l_供給側の例外?.Throw();
        }

        #endregion
    }
}
