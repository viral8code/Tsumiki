using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// プロデューサー/コンシューマが、ワーカーの例外で無言のハングに
    /// 陥らないことを固定する。
    ///
    /// 素朴に書くとこうなる: キューには容量上限があるため、ワーカーが例外で
    /// 落ちるとキューを引き取る者がいなくなり、プロデューサーは Add で永久に
    /// 待ち続ける。Task.WaitAll に到達しないのでワーカーの例外は誰にも
    /// 観測されず、ログも例外も出ないままプロセスが CPU 0% で止まる。
    /// 実際に GAGE-B のデータで2時間以上まったく無言でハングした。
    /// </summary>
    public class ReadPipelineTests
    {
        [Fact]
        public void Run_ProcessesEveryItemExactlyOnce()
        {
            var 入力 = Enumerable.Range(0, 5000).ToList();
            var 結果 = new System.Collections.Concurrent.ConcurrentBag<int>();

            ReadPipeline.V_実行(4, 32, 入力, (l_項目, _) => 結果.Add(l_項目));

            Assert.Equal(入力.Count, 結果.Count);
            Assert.Equal(入力, 結果.OrderBy(x => x).ToList());
        }

        [Fact]
        public void Run_PassesTheWorkerIndexWithinRange()
        {
            const int スレッド数 = 4;
            var 観測した番号 = new System.Collections.Concurrent.ConcurrentBag<int>();

            ReadPipeline.V_実行(スレッド数, 16, Enumerable.Range(0, 500), (_, l_番号) => 観測した番号.Add(l_番号));

            Assert.All(観測した番号, l_番号 => Assert.InRange(l_番号, 0, スレッド数 - 1));
        }

        /// <summary>
        /// ワーカーが例外を投げたら、供給が残っていても呼び出し元へ伝わること。
        /// 入力数はキュー容量よりずっと多くしてあり、対策が無ければ
        /// プロデューサーが満杯のキューで待ち続けてこのテストはタイムアウトする。
        /// </summary>
        [Fact]
        public void Run_WorkerThrows_SurfacesTheExceptionInsteadOfHanging()
        {
            var 例外 = Assert.Throws<AggregateException>(() =>
                ReadPipeline.V_実行(4, 8, Enumerable.Range(0, 100_000), (l_項目, _) =>
                {
                    if (l_項目 >= 0)
                    {
                        throw new InvalidOperationException("worker failed");
                    }
                }));

            _ = Assert.IsType<InvalidOperationException>(例外.InnerExceptions[0]);
            Assert.Equal("worker failed", 例外.InnerExceptions[0].Message);
        }

        /// <summary>
        /// 一部のワーカーだけが落ちた場合も、放置せずに伝えること。
        /// 残ったワーカーが処理を続けられてしまうと、結果が中途半端なまま
        /// 「成功」として先へ進んでしまう。
        /// </summary>
        [Fact]
        public void Run_OneWorkerThrows_StillSurfacesTheException()
        {
            var 処理数 = 0;

            var 例外 = Assert.Throws<AggregateException>(() =>
                ReadPipeline.V_実行(4, 8, Enumerable.Range(0, 100_000), (l_項目, _) =>
                {
                    if (Interlocked.Increment(ref 処理数) == 50)
                    {
                        throw new InvalidOperationException("one worker failed");
                    }
                }));

            Assert.Contains(例外.InnerExceptions, x => x.Message == "one worker failed");
        }

        /// <summary>
        /// 供給側が例外を投げた場合も、スタックトレースを保ったまま伝わること
        /// (壊れた FASTQ を読んだ場合などがこれに当たる)。
        /// </summary>
        [Fact]
        public void Run_ProducerThrows_SurfacesTheException()
        {
            var 例外 = Assert.Throws<FormatException>(() =>
                ReadPipeline.V_実行(4, 8, Get_途中で壊れる入力(), (_, _) => { }));

            Assert.Equal("broken input", 例外.Message);
        }

        private static IEnumerable<int> Get_途中で壊れる入力()
        {
            for (var i = 0; i < 100; i++)
            {
                yield return i;
            }
            throw new FormatException("broken input");
        }
    }
}
