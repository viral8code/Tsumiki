using System.Collections.Concurrent;
using System.Runtime.ExceptionServices;

namespace Tsumiki.Utility
{
    /// <summary>
    /// リードを1本のスレッドで順に読み進めつつ(ディスクI/Oはシーケンシャルなまま)、
    /// ワーカースレッド群へ配って並列に処理するプロデューサー/コンシューマ。
    ///
    /// このクラスが存在する理由は、素朴に書いた場合の障害の出方が最悪だからである。
    /// キューには容量上限があるため、ワーカーが例外で落ちるとキューを引き取る者が
    /// いなくなり、プロデューサーは Add で永久に待ち続ける。Task.WaitAll に
    /// 到達しないのでワーカーの例外は誰にも観測されず、ログも例外も出ないまま
    /// プロセスが CPU 0% で停止する。
    ///
    /// 実際、リードのマッピングでリード長が k 未満のときに添字が範囲外になる
    /// 不具合があり、GAGE-B のトリミング済みリード(8%以上が k 未満)で
    /// 2時間以上まったく無言でハングした。原因の特定にプロセスの CPU 時間と
    /// ページフォルト数を測るところから始める羽目になっており、
    /// 「落ちるべきときに落ちる」ようにしておく価値は大きい。
    /// </summary>
    internal static class ReadPipeline
    {
        /// <summary>
        /// p_供給元 の各要素を p_処理 へ並列に配る。
        /// p_処理 の第2引数はワーカー番号(0 以上 p_スレッド数 未満)で、
        /// ワーカーごとのローカル集計用配列の添字として使うことを想定している。
        ///
        /// ワーカー・供給元のいずれかが例外を投げた場合は、その例外を
        /// スタックトレースを保ったまま呼び出し元へ送出する。
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
                        // ワーカーが落ちたことを供給側へ伝える。伝えないと、
                        // 満杯のキューへの Add で永久に待つことになる。
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
                // ワーカー側の例外が真の原因。下の WaitAll で送出される。
            }
            catch (Exception ex)
            {
                l_供給側の例外 = ExceptionDispatchInfo.Capture(ex);
            }

            // 供給が途中で終わった場合でも、待っているワーカーを必ず解放する。
            l_キュー.CompleteAdding();

            // ワーカー側の例外があればここで送出される。供給側の例外より
            // ワーカー側を優先するのは、供給側が中断されたのは
            // ワーカーが落ちた結果であることが多いため。
            Task.WaitAll(l_ワーカー);

            l_供給側の例外?.Throw();
        }
    }
}
