using System.Collections.Concurrent;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Correction;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Preprocessing
{
    /// <summary>
    /// k-mer スペクトラムに基づく、Quake/BayesHammer 類似の簡易リードエラー訂正
    /// </summary>
    internal static class ErrorCorrector
    {
        #region 定数

        /// <summary>
        /// 1 バッチあたりのリード数
        /// </summary>
        internal const int C_訂正バッチサイズ = 20_000;

        /// <summary>
        /// 読み込み・訂正・書き出しの間に溜めておくバッチ数
        /// </summary>
        private const int C_先読みするバッチ数 = 2;

        #endregion

        #region 内部変数

        /// <summary>
        /// スレッドごとの訂正の作業域
        /// </summary>
        [ThreadStatic]
        private static 訂正作業域? _作業域;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// リードファイルを読み込んでエラー訂正を行い、結果を出力先へ書き出す
        /// </summary>
        /// <param name="p_リード1のパス"></param>
        /// <param name="p_リード2のパス"></param>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_出力先1"></param>
        /// <param name="p_出力先2"></param>
        public static void V_訂正_リードファイル(string p_リード1のパス, string? p_リード2のパス, string p_一時ディレクトリ, string p_出力先1, string? p_出力先2, int p_Phredオフセット)
        {
            using var l_計測 = new StageTimer("error-correction");
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            var l_訂正用一時ディレクトリ = Path.Combine(p_一時ディレクトリ, "error_correction");
            _ = Directory.CreateDirectory(l_訂正用一時ディレクトリ);

            Logger.V_出力(メッセージID.エラー訂正_スペクトル構築);
            using (var l_kmerインデックス = new TrustedKmerIndex(l_訂正用一時ディレクトリ))
            {
                if (p_リード2のパス != null)
                {
                    Parallel.Invoke(
                        () => KmerCounting.V_読込_リードファイル(p_リード1のパス, l_kmerインデックス, p_Phredオフセット),
                        () => KmerCounting.V_読込_リードファイル(p_リード2のパス, l_kmerインデックス, p_Phredオフセット));
                }
                else
                {
                    KmerCounting.V_読込_リードファイル(p_リード1のパス, l_kmerインデックス, p_Phredオフセット);
                }

                KmerCutoffSelector.V_解決_kmerカットオフ(ConfigurationManager.A_実行時引数, l_kmerインデックス);
                l_kmerインデックス.V_適用_カットオフ(ConfigurationManager.A_実行時引数.A_kmerカットオフ);

                Logger.V_出力(メッセージID.エラー訂正_訂正開始);
                var l_統計1 = Get_訂正統計_ファイル(p_リード1のパス, p_出力先1, l_kmerインデックス, l_k長);
                Logger.V_出力(メッセージID.エラー訂正_ファイル別統計, Path.GetFileName(p_リード1のパス), l_統計1.A_訂正されたリード数, l_統計1.A_総リード数, l_統計1.A_総訂正塩基数);

                if (p_リード2のパス != null && p_出力先2 != null)
                {
                    var l_統計2 = Get_訂正統計_ファイル(p_リード2のパス, p_出力先2, l_kmerインデックス, l_k長);
                    Logger.V_出力(メッセージID.エラー訂正_ファイル別統計, Path.GetFileName(p_リード2のパス), l_統計2.A_訂正されたリード数, l_統計2.A_総リード数, l_統計2.A_総訂正塩基数);
                }
            }

            Directory.Delete(l_訂正用一時ディレクトリ, recursive: true);
        }

        /// <summary>
        /// 1 リード (塩基 ID 空間のバイト列、曖昧塩基は <see cref="Consts.無効な塩基"/>) を貪欲法で訂正する
        /// </summary>
        /// <param name="p_リード"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_最大反復数"></param>
        /// <returns></returns>
        public static 訂正結果 Get_訂正結果(ReadOnlySpan<byte> p_リード, TrustedKmerIndex p_kmerインデックス, int p_k長, int p_最大反復数 = 10)
        {
            return p_リード.Length < p_k長
                ? new 訂正結果(p_リード.ToArray(), 0)
                : p_k長 <= 64
                ? Get_訂正結果_パック(p_リード.ToArray(), p_kmerインデックス, p_k長, p_最大反復数)
                : Get_訂正結果_逐次(p_リード.ToArray(), p_kmerインデックス, p_k長, p_最大反復数);
        }

        /// <summary>
        /// k が 64 を超える場合の経路
        /// </summary>
        /// <param name="p_塩基列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_最大反復数"></param>
        /// <returns></returns>
        internal static 訂正結果 Get_訂正結果_逐次(byte[] p_塩基列, TrustedKmerIndex p_kmerインデックス, int p_k長, int p_最大反復数)
        {
            var l_塩基列 = p_塩基列;
            var l_窓数 = l_塩基列.Length - p_k長 + 1;
            var l_訂正数 = 0;

            for (var l_反復 = 0; l_反復 < p_最大反復数; l_反復++)
            {
                var l_信頼状況 = Get_窓別信頼状況(l_塩基列, p_k長, p_kmerインデックス);
                if (Array.TrueForAll(l_信頼状況, x => x))
                {
                    break;
                }

                var l_最良位置 = -1;
                byte l_最良塩基 = 0;
                var l_最良改善数 = 0;

                for (var l_位置 = 0; l_位置 < l_塩基列.Length; l_位置++)
                {
                    if (l_塩基列[l_位置] == Consts.無効な塩基)
                    {
                        continue;
                    }

                    var l_窓開始 = Math.Max(0, l_位置 - p_k長 + 1);
                    var l_窓終了 = Math.Min(l_窓数 - 1, l_位置);

                    var l_Has未信頼窓 = false;
                    for (var w = l_窓開始; w <= l_窓終了; w++)
                    {
                        if (!l_信頼状況[w])
                        {
                            l_Has未信頼窓 = true;
                            break;
                        }
                    }

                    if (!l_Has未信頼窓)
                    {
                        continue;
                    }

                    var l_現在の塩基 = l_塩基列[l_位置];
                    for (var l_候補 = Consts.塩基ID.A; l_候補 <= Consts.塩基ID.T; l_候補++)
                    {
                        if (l_候補 == l_現在の塩基)
                        {
                            continue;
                        }

                        var l_改善数 = Get_置換改善数(l_塩基列, l_位置, l_候補, l_窓開始, l_窓終了, p_k長, l_信頼状況, p_kmerインデックス);
                        if (l_改善数 > l_最良改善数)
                        {
                            l_最良改善数 = l_改善数;
                            l_最良位置 = l_位置;
                            l_最良塩基 = l_候補;
                        }
                    }
                }

                if (l_最良位置 < 0)
                {
                    break;
                }

                l_塩基列[l_最良位置] = l_最良塩基;
                l_訂正数++;
            }

            return new 訂正結果(l_塩基列, l_訂正数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 ファイルぶんのリードを訂正して書き出し、その集計を返す
        /// </summary>
        /// <param name="p_入力パス">訂正するリードのパス</param>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>訂正の集計</returns>
        private static ファイル訂正統計 Get_訂正統計_ファイル(string p_入力パス, string p_出力パス, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_総リード数 = 0;
            var l_訂正されたリード数 = 0;
            var l_総訂正塩基数 = 0;

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            using var l_中断 = new CancellationTokenSource();
            using var l_読込済み = new BlockingCollection<訂正バッチ>(C_先読みするバッチ数);
            using var l_訂正済み = new BlockingCollection<訂正バッチ>(C_先読みするバッチ数);

            var l_読み込み = Task.Run(() => V_読込_バッチ群(p_入力パス, l_読込済み, l_中断));
            var l_書き出し = Task.Run(() => V_書出_バッチ群(p_出力パス, l_訂正済み, l_中断));

            try
            {
                foreach (var l_バッチ in l_読込済み.GetConsumingEnumerable(l_中断.Token))
                {
                    l_総リード数 += l_バッチ.A_件数;

                    _ = Parallel.For(0, l_バッチ.A_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, i =>
                    {
                        l_バッチ.A_結果群[i] = Get_訂正結果(l_バッチ.A_塩基列群[i], p_kmerインデックス, p_k長);
                    });

                    for (var i = 0; i < l_バッチ.A_件数; i++)
                    {
                        if (l_バッチ.A_結果群[i].A_訂正数 > 0)
                        {
                            l_訂正されたリード数++;
                            l_総訂正塩基数 += l_バッチ.A_結果群[i].A_訂正数;
                        }
                    }

                    l_訂正済み.Add(l_バッチ, l_中断.Token);
                }
            }
            catch (OperationCanceledException) when (l_中断.IsCancellationRequested)
            {
            }
            catch
            {
                l_中断.Cancel();
                throw;
            }
            finally
            {
                l_訂正済み.CompleteAdding();
            }

            l_読み込み.GetAwaiter().GetResult();
            l_書き出し.GetAwaiter().GetResult();

            return new ファイル訂正統計(l_総リード数, l_訂正されたリード数, l_総訂正塩基数);
        }

        /// <summary>
        /// リードをバッチに詰めて順に渡す
        /// </summary>
        /// <param name="p_入力パス">訂正するリードのパス</param>
        /// <param name="p_渡し先">読み込んだバッチの渡し先</param>
        /// <param name="p_中断">どこかで失敗したときに止める合図</param>
        private static void V_読込_バッチ群(string p_入力パス, BlockingCollection<訂正バッチ> p_渡し先, CancellationTokenSource p_中断)
        {
            try
            {
                using var l_読み込み = new FastqReader(p_入力パス);
                while (l_読み込み.Has続き())
                {
                    var l_バッチ = new 訂正バッチ();
                    while (l_バッチ.A_件数 < C_訂正バッチサイズ && l_読み込み.Has続き())
                    {
                        var l_リード = l_読み込み.Get_次のリード_軽量();
                        l_バッチ.A_ID群[l_バッチ.A_件数] = l_リード.A_ID;
                        l_バッチ.A_クオリティ群[l_バッチ.A_件数] = l_リード.A_クオリティ;
                        l_バッチ.A_塩基列群[l_バッチ.A_件数] = l_リード.A_塩基列!;
                        l_バッチ.A_件数++;
                    }

                    p_渡し先.Add(l_バッチ, p_中断.Token);
                }
            }
            catch (OperationCanceledException) when (p_中断.IsCancellationRequested)
            {
            }
            catch
            {
                p_中断.Cancel();
                throw;
            }
            finally
            {
                p_渡し先.CompleteAdding();
            }
        }

        /// <summary>
        /// 訂正を終えたバッチを、受け取った順に書き出す
        /// </summary>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_受け取り元">訂正を終えたバッチ</param>
        /// <param name="p_中断">どこかで失敗したときに止める合図</param>
        private static void V_書出_バッチ群(string p_出力パス, BlockingCollection<訂正バッチ> p_受け取り元, CancellationTokenSource p_中断)
        {
            try
            {
                using var l_書き込み = new FastqWriter(p_出力パス);
                foreach (var l_バッチ in p_受け取り元.GetConsumingEnumerable(p_中断.Token))
                {
                    for (var i = 0; i < l_バッチ.A_件数; i++)
                    {
                        l_書き込み.V_書き込み(l_バッチ.A_ID群[i], string.Create(l_バッチ.A_結果群[i].A_塩基列.Length, l_バッチ.A_結果群[i].A_塩基列, static (l_文字, l_塩基列) =>
                        {
                            for (var j = 0; j < l_塩基列.Length; j++)
                            {
                                l_文字[j] = Util.Get_塩基文字(l_塩基列[j]);
                            }
                        }), l_バッチ.A_クオリティ群[i]);
                    }
                }
            }
            catch (OperationCanceledException) when (p_中断.IsCancellationRequested)
            {
            }
            catch
            {
                p_中断.Cancel();
                throw;
            }
        }

        /// <summary>
        /// パック経路の訂正本体
        /// </summary>
        /// <param name="p_塩基列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_最大反復数"></param>
        /// <returns></returns>
        private static 訂正結果 Get_訂正結果_パック(byte[] p_塩基列, TrustedKmerIndex p_kmerインデックス, int p_k長, int p_最大反復数)
        {
            var l_窓数 = p_塩基列.Length - p_k長 + 1;
            var l_作業域 = _作業域 ??= new 訂正作業域();
            l_作業域.V_確保(l_窓数);
            var l_パック = l_作業域.A_パック.AsSpan(0, l_窓数);
            var l_逆相補 = l_作業域.A_逆相補.AsSpan(0, l_窓数);
            var l_無効数 = l_作業域.A_無効数.AsSpan(0, l_窓数);
            var l_信頼状況 = l_作業域.A_信頼状況.AsSpan(0, l_窓数);

            var l_未信頼累積 = l_作業域.A_未信頼累積.AsSpan(0, l_窓数 + 1);
            var l_訂正数 = 0;

            for (var l_反復 = 0; l_反復 < p_最大反復数; l_反復++)
            {
                V_計算_窓状態(p_塩基列, p_k長, p_kmerインデックス, l_パック, l_逆相補, l_無効数, l_信頼状況);

                l_未信頼累積[0] = 0;
                for (var w = 0; w < l_窓数; w++)
                {
                    l_未信頼累積[w + 1] = l_未信頼累積[w] + (l_信頼状況[w] ? 0 : 1);
                }

                if (l_未信頼累積[l_窓数] == 0)
                {
                    break;
                }

                var l_最良位置 = -1;
                byte l_最良塩基 = 0;
                var l_最良改善数 = 0;

                for (var l_位置 = 0; l_位置 < p_塩基列.Length; l_位置++)
                {
                    if (p_塩基列[l_位置] == Consts.無効な塩基)
                    {
                        continue;
                    }

                    var l_窓開始 = Math.Max(0, l_位置 - p_k長 + 1);
                    var l_窓終了 = Math.Min(l_窓数 - 1, l_位置);
                    if (l_未信頼累積[l_窓終了 + 1] - l_未信頼累積[l_窓開始] == 0)
                    {
                        continue;
                    }

                    var l_現在の塩基 = p_塩基列[l_位置];
                    for (var l_候補 = Consts.塩基ID.A; l_候補 <= Consts.塩基ID.T; l_候補++)
                    {
                        if (l_候補 == l_現在の塩基)
                        {
                            continue;
                        }

                        var l_改善数 = Get_置換改善数_パック(l_位置, l_候補, l_窓開始, l_窓終了, p_k長, p_kmerインデックス, l_パック, l_逆相補, l_無効数, l_信頼状況, l_未信頼累積, l_最良改善数);
                        if (l_改善数 > l_最良改善数)
                        {
                            l_最良改善数 = l_改善数;
                            l_最良位置 = l_位置;
                            l_最良塩基 = l_候補;
                        }
                    }
                }

                if (l_最良位置 < 0)
                {
                    break;
                }

                p_塩基列[l_最良位置] = l_最良塩基;
                l_訂正数++;
            }

            return new 訂正結果(p_塩基列, l_訂正数);
        }

        /// <summary>
        /// 全窓のパック値・逆相補・曖昧塩基の数・信頼状況を、隣の窓から転がして求める
        /// </summary>
        /// <param name="p_塩基列"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_パック"></param>
        /// <param name="p_逆相補"></param>
        /// <param name="p_無効数"></param>
        /// <param name="p_信頼状況"></param>
        private static void V_計算_窓状態(byte[] p_塩基列, int p_k長, TrustedKmerIndex p_kmerインデックス, Span<UInt128> p_パック, Span<UInt128> p_逆相補, Span<int> p_無効数, Span<bool> p_信頼状況)
        {
            var l_マスク = Get_マスク(p_k長);
            var l_最上位への移動 = 2 * (p_k長 - 1);

            UInt128 l_パック = 0;
            UInt128 l_逆相補 = 0;
            var l_無効数 = 0;
            for (var i = 0; i < p_k長; i++)
            {
                l_パック = ((l_パック << 2) | Get_コドン(p_塩基列[i])) & l_マスク;
                l_逆相補 = (l_逆相補 >> 2) | (Get_相補コドン(p_塩基列[i]) << l_最上位への移動);
                if (p_塩基列[i] == Consts.無効な塩基)
                {
                    l_無効数++;
                }
            }

            for (var w = 0; w < p_パック.Length; w++)
            {
                if (w > 0)
                {
                    var l_出る塩基 = p_塩基列[w - 1];
                    var l_入る塩基 = p_塩基列[w + p_k長 - 1];
                    l_パック = ((l_パック << 2) | Get_コドン(l_入る塩基)) & l_マスク;
                    l_逆相補 = (l_逆相補 >> 2) | (Get_相補コドン(l_入る塩基) << l_最上位への移動);
                    if (l_出る塩基 == Consts.無効な塩基)
                    {
                        l_無効数--;
                    }

                    if (l_入る塩基 == Consts.無効な塩基)
                    {
                        l_無効数++;
                    }
                }

                p_パック[w] = l_パック;
                p_逆相補[w] = l_逆相補;
                p_無効数[w] = l_無効数;
                p_信頼状況[w] = l_無効数 == 0
                    && Haskmer(p_kmerインデックス, p_k長, l_パック, l_逆相補);
            }
        }

        /// <summary>
        /// p_位置 を p_候補 に置換したときの、信頼できる窓の純増数
        /// </summary>
        /// <param name="p_位置"></param>
        /// <param name="p_候補"></param>
        /// <param name="p_窓開始"></param>
        /// <param name="p_窓終了"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_パック"></param>
        /// <param name="p_逆相補"></param>
        /// <param name="p_無効数"></param>
        /// <param name="p_信頼状況"></param>
        /// <param name="p_未信頼累積"></param>
        /// <param name="p_最良改善数"></param>
        /// <returns></returns>
        private static int Get_置換改善数_パック(int p_位置, byte p_候補, int p_窓開始, int p_窓終了, int p_k長, TrustedKmerIndex p_kmerインデックス, ReadOnlySpan<UInt128> p_パック, ReadOnlySpan<UInt128> p_逆相補, ReadOnlySpan<int> p_無効数, ReadOnlySpan<bool> p_信頼状況, ReadOnlySpan<int> p_未信頼累積, int p_最良改善数)
        {
            var l_コドン = Get_コドン(p_候補);
            var l_相補コドン = Get_相補コドン(p_候補);
            var l_改善数 = 0;

            for (var w = p_窓開始; w <= p_窓終了; w++)
            {
                var l_残りの上界 = p_未信頼累積[p_窓終了 + 1] - p_未信頼累積[w];
                if (l_改善数 + l_残りの上界 <= p_最良改善数)
                {
                    return l_改善数;
                }

                bool l_Is信頼可能;
                if (p_無効数[w] > 0)
                {
                    l_Is信頼可能 = false;
                }
                else
                {
                    var l_窓内の位置 = p_位置 - w;
                    var l_移動 = 2 * (p_k長 - 1 - l_窓内の位置);
                    var l_逆相補の移動 = 2 * l_窓内の位置;
                    var l_パック = (p_パック[w] & ~((UInt128)3 << l_移動)) | (l_コドン << l_移動);
                    var l_逆相補 =
                        (p_逆相補[w] & ~((UInt128)3 << l_逆相補の移動)) | (l_相補コドン << l_逆相補の移動);
                    l_Is信頼可能 = Haskmer(p_kmerインデックス, p_k長, l_パック, l_逆相補);
                }

                if (l_Is信頼可能 && !p_信頼状況[w])
                {
                    l_改善数++;
                }
                else if (!l_Is信頼可能 && p_信頼状況[w])
                {
                    l_改善数--;
                }
            }

            return l_改善数;
        }

        /// <summary>
        /// パック済みの順鎖・逆鎖から、信頼できる k-mer 集合に含まれるかを引く
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_パック"></param>
        /// <param name="p_逆相補"></param>
        /// <returns></returns>
        private static bool Haskmer(TrustedKmerIndex p_kmerインデックス, int p_k長, UInt128 p_パック, UInt128 p_逆相補)
        {
            var l_正規形 = p_パック < p_逆相補 ? p_パック : p_逆相補;
            return p_k長 <= 32
                ? p_kmerインデックス.Haskmer_小((ulong)l_正規形)
                : p_kmerインデックス.Haskmer_中(l_正規形);
        }

        /// <summary>
        /// k 塩基ぶんだけを残すマスクを返す
        /// </summary>
        /// <param name="p_k長">k 長</param>
        /// <returns>マスク</returns>
        private static UInt128 Get_マスク(int p_k長)
        {
            return p_k長 >= 64 ? UInt128.MaxValue : ((UInt128)1 << (2 * p_k長)) - 1;
        }

        /// <summary>
        /// 塩基 ID の 2 bit 表現
        /// </summary>
        /// <param name="p_塩基ID"></param>
        /// <returns></returns>
        private static UInt128 Get_コドン(byte p_塩基ID)
        {
            return p_塩基ID == Consts.無効な塩基 ? 0 : (UInt128)(p_塩基ID - 1);
        }

        /// <summary>
        /// 塩基 ID に対応する相補の 2 bit を返す
        /// </summary>
        /// <param name="p_塩基ID">塩基 ID</param>
        /// <returns>相補の 2 bit</returns>
        private static UInt128 Get_相補コドン(byte p_塩基ID)
        {
            return p_塩基ID == Consts.無効な塩基 ? 3 : (UInt128)(3 - (p_塩基ID - 1));
        }

        /// <summary>
        /// リード上の各窓が信頼できる k-mer かを返す
        /// </summary>
        /// <param name="p_塩基列">リードの塩基 ID 列</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <returns>窓ごとに信頼できるか</returns>
        private static bool[] Get_窓別信頼状況(byte[] p_塩基列, int p_k長, TrustedKmerIndex p_kmerインデックス)
        {
            var l_窓数 = p_塩基列.Length - p_k長 + 1;
            var l_信頼状況 = new bool[l_窓数];
            for (var w = 0; w < l_窓数; w++)
            {
                l_信頼状況[w] = Is信頼窓(p_塩基列, w, p_k長, p_kmerインデックス);
            }

            return l_信頼状況;
        }

        /// <summary>
        /// その窓の k-mer が信頼できる集合にあるか
        /// </summary>
        /// <param name="p_塩基列">リードの塩基 ID 列</param>
        /// <param name="p_窓開始">窓の開始位置</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <returns>信頼できれば true</returns>
        private static bool Is信頼窓(byte[] p_塩基列, int p_窓開始, int p_k長, TrustedKmerIndex p_kmerインデックス)
        {
            for (var i = p_窓開始; i < p_窓開始 + p_k長; i++)
            {
                if (p_塩基列[i] == Consts.無効な塩基)
                {
                    return false;
                }
            }

            return p_kmerインデックス.Haskmer(p_塩基列.AsSpan(p_窓開始, p_k長));
        }

        /// <summary>
        /// p_位置 を p_候補 に置換した場合の「信頼できる窓の純増数」を計算する (p_窓開始..p_窓終了 の範囲、すなわち その位置を含みうる窓のみが影響を受けるため、その範囲だけを再評価すれば十分)
        /// </summary>
        /// <param name="p_塩基列"></param>
        /// <param name="p_位置"></param>
        /// <param name="p_候補"></param>
        /// <param name="p_窓開始"></param>
        /// <param name="p_窓終了"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_置換前の信頼状況"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <returns></returns>
        private static int Get_置換改善数(byte[] p_塩基列, int p_位置, byte p_候補, int p_窓開始, int p_窓終了, int p_k長, bool[] p_置換前の信頼状況, TrustedKmerIndex p_kmerインデックス)
        {
            var l_元の塩基 = p_塩基列[p_位置];
            p_塩基列[p_位置] = p_候補;

            var l_改善数 = 0;
            for (var w = p_窓開始; w <= p_窓終了; w++)
            {
                var l_Is信頼可能 = Is信頼窓(p_塩基列, w, p_k長, p_kmerインデックス);
                if (l_Is信頼可能 && !p_置換前の信頼状況[w])
                {
                    l_改善数++;
                }
                else if (!l_Is信頼可能 && p_置換前の信頼状況[w])
                {
                    l_改善数--;
                }
            }

            p_塩基列[p_位置] = l_元の塩基;
            return l_改善数;
        }

        #endregion

    }
}
