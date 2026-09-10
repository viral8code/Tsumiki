using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// プロデューサー/コンシューマが、ワーカーの例外で無言のハングに
    /// 陥らないことを固定する
    /// </summary>
    /// <remarks>
    /// 素朴に書くとこうなる: キューには容量上限があるため、ワーカーが例外で
    /// 落ちるとキューを引き取る者がいなくなり、プロデューサーは Add で永久に
    /// 待ち続ける<br/>
    /// Task.WaitAll に到達しないのでワーカーの例外は誰にも
    /// 観測されず、ログも例外も出ないままプロセスが CPU 0% で止まる<br/>
    /// 実際に GAGE-B のデータで 2 時間以上まったく無言でハングした
    /// </remarks>
    public class ReadPipelineTests
    {
        /// <summary>
        /// 全ての項目をちょうど 1 回ずつ処理する
        /// </summary>
        [Fact]
        public void 全ての項目をちょうど1回ずつ処理する()
        {
            var l_入力 = Enumerable.Range(0, 5000).ToList();
            var l_結果 = new System.Collections.Concurrent.ConcurrentBag<int>();

            ReadPipeline.V_実行(4, 32, l_入力, (l_項目, _) => l_結果.Add(l_項目));

            Assert.Equal(l_入力.Count, l_結果.Count);
            Assert.Equal(l_入力, l_結果.OrderBy(x => x).ToList());
        }

        /// <summary>
        /// ワーカー番号を範囲内で渡す
        /// </summary>
        [Fact]
        public void ワーカー番号を範囲内で渡す()
        {
            const int l_スレッド数 = 4;
            var l_観測した番号 = new System.Collections.Concurrent.ConcurrentBag<int>();

            ReadPipeline.V_実行(l_スレッド数, 16, Enumerable.Range(0, 500), (_, l_番号) => l_観測した番号.Add(l_番号));

            Assert.All(l_観測した番号, l_番号 => Assert.InRange(l_番号, 0, l_スレッド数 - 1));
        }

        /// <summary>
        /// ワーカーが例外を投げたら、供給が残っていても呼び出し元へ伝わること
        /// </summary>
        /// <remarks>
        /// 入力数はキュー容量よりずっと多くしてあり、対策が無ければ
        /// プロデューサーが満杯のキューで待ち続けてこのテストはタイムアウトする
        /// </remarks>
        [Fact]
        public void ワーカーが例外を投げると呼び出し元に伝わる()
        {
            var l_例外 = Assert.Throws<AggregateException>(() =>
                ReadPipeline.V_実行(4, 8, Enumerable.Range(0, 100_000), (l_項目, _) =>
                {
                    if (l_項目 >= 0)
                    {
                        throw new InvalidOperationException("worker failed");
                    }
                }));

            _ = Assert.IsType<InvalidOperationException>(l_例外.InnerExceptions[0]);
            Assert.Equal("worker failed", l_例外.InnerExceptions[0].Message);
        }

        /// <summary>
        /// 一部のワーカーだけが落ちた場合も、放置せずに伝えること
        /// </summary>
        /// <remarks>
        /// 残ったワーカーが処理を続けられてしまうと、結果が中途半端なまま
        /// 「成功」として先へ進んでしまう
        /// </remarks>
        [Fact]
        public void 一部のワーカーだけが例外を投げても伝わる()
        {
            var l_処理数 = 0;

            var l_例外 = Assert.Throws<AggregateException>(() =>
                ReadPipeline.V_実行(4, 8, Enumerable.Range(0, 100_000), (l_項目, _) =>
                {
                    if (Interlocked.Increment(ref l_処理数) == 50)
                    {
                        throw new InvalidOperationException("one worker failed");
                    }
                }));

            Assert.Contains(l_例外.InnerExceptions, x => x.Message == "one worker failed");
        }

        /// <summary>
        /// 供給側が例外を投げた場合も、スタックトレースを保ったまま伝わること
        /// (壊れた FASTQ を読んだ場合などがこれに当たる)
        /// </summary>
        [Fact]
        public void 供給側が例外を投げても伝わる()
        {
            var l_例外 = Assert.Throws<FormatException>(() =>
                ReadPipeline.V_実行(4, 8, Get_途中で壊れる入力(), (_, _) => { }));

            Assert.Equal("broken input", l_例外.Message);
        }

        /// <summary>
        /// 途中で例外を投げる入力を返す
        /// </summary>
        /// <returns>途中で壊れる入力</returns>
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
