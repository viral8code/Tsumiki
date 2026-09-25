using System.Collections.Concurrent;
using System.Runtime.ExceptionServices;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// リードを 1 本のスレッドで順に読み進めつつ、ワーカー群へ配って並列に処理する
    /// </summary>
    internal static class ReadPipeline
    {
        #region 定数

        /// <summary>
        /// ワーカーへ 1 回に渡す要素数
        /// </summary>
        private const int C_束の要素数 = 256;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_供給元 の各要素を p_処理 へ並列に配る
        /// </summary>
        /// <param name="p_スレッド数"></param>
        /// <param name="p_キュー容量">溜めておく要素数の目安</param>
        /// <param name="p_供給元"></param>
        /// <param name="p_処理"></param>
        public static void V_実行<T>(int p_スレッド数, int p_キュー容量, IEnumerable<T> p_供給元, Action<T, int> p_処理)
        {
            using var l_キュー = new BlockingCollection<T[]>(Math.Max(2 * p_スレッド数, p_キュー容量 / C_束の要素数));
            using var l_中断 = new CancellationTokenSource();

            var l_ワーカー = new Task[p_スレッド数];
            for (var w = 0; w < p_スレッド数; w++)
            {
                var l_ワーカー番号 = w;
                l_ワーカー[w] = Task.Run(() =>
                {
                    try
                    {
                        foreach (var l_束 in l_キュー.GetConsumingEnumerable())
                        {
                            foreach (var l_項目 in l_束)
                            {
                                p_処理(l_項目, l_ワーカー番号);
                            }
                        }
                    }
                    catch
                    {
                        l_中断.Cancel();
                        throw;
                    }
                });
            }

            ExceptionDispatchInfo? l_供給側の例外 = null;
            try
            {
                foreach (var l_束 in p_供給元.Chunk(C_束の要素数))
                {
                    l_キュー.Add(l_束, l_中断.Token);
                }
            }
            catch (OperationCanceledException)
            {
            }
            catch (Exception l_例外)
            {
                l_供給側の例外 = ExceptionDispatchInfo.Capture(l_例外);
            }

            l_キュー.CompleteAdding();

            Task.WaitAll(l_ワーカー);

            l_供給側の例外?.Throw();
        }

        #endregion
    }
}
