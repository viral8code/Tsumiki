using System.Collections.Concurrent;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Preprocessing;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Preprocessing
{
    /// <summary>
    /// ペアエンドの 2 本 (R1 と RC (R2)) を重ね合わせて、アダプタリードスルーのトリムと、高信頼と低信頼が明確に分かれる位置でのペア相互訂正を行う
    /// </summary>
    internal static class Preprocessor
    {
        #region 定数

        /// <summary>
        /// 項目 N0
        /// </summary>
        private const string C_項目_N0 = "N0";

        /// <summary>
        /// 項目 F1
        /// </summary>
        private const string C_項目_F1 = "F1";

        /// <summary>
        /// ペアの重なりとみなすために要求する最小長
        /// </summary>
        private const int C_最小オーバーラップ長 = 30;

        /// <summary>
        /// この不一致率までは同一断片から重なって読んだものとみなす
        /// </summary>
        private const double C_許容不一致率 = 0.2D;

        /// <summary>
        /// 相互訂正で高信頼とみなす最小 Phred スコア
        /// </summary>
        private const int C_高信頼スコア = 30;

        /// <summary>
        /// 相互訂正で低信頼とみなす最大 Phred スコア
        /// </summary>
        private const int C_低信頼スコア = 14;

        /// <summary>
        /// 1 バッチあたりのペア数
        /// </summary>
        internal const int C_前処理バッチサイズ = 20_000;

        /// <summary>
        /// 読み込みと書き出しで先に溜めておくバッチの数
        /// </summary>
        private const int C_先読みするバッチ数 = 2;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// ペアの FASTQ を読み込んで前処理し、結果を出力先へ書き出す
        /// </summary>
        /// <param name="p_順リードのパス">リード 1 のパス</param>
        /// <param name="p_逆リードのパス">リード 2 のパス</param>
        /// <param name="p_出力先1">前処理したリード 1 の書き出し先</param>
        /// <param name="p_出力先2">前処理したリード 2 の書き出し先</param>
        /// <param name="p_Phredオフセット"></param>
        /// <returns>前処理の集計</returns>
        public static 前処理統計 V_前処理_リードファイル(string p_順リードのパス, string p_逆リードのパス, string p_出力先1, string p_出力先2, int p_Phredオフセット)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_トリム閾値 = ConfigurationManager.A_実行時引数.A_品質トリム閾値;
            var l_総ペア数 = 0;
            var l_アダプタ検出ペア数 = 0;
            var l_訂正塩基数 = 0;
            var l_総塩基数 = 0L;
            var l_トリム塩基数 = 0L;
            using var l_中断 = new CancellationTokenSource();
            using var l_読込済み = new BlockingCollection<前処理バッチ>(C_先読みするバッチ数);
            using var l_処理済み = new BlockingCollection<前処理バッチ>(C_先読みするバッチ数);
            var l_読み込み = Task.Run(() => V_読込_バッチ群(p_順リードのパス, p_逆リードのパス, l_読込済み, l_中断));
            var l_書き出し = Task.Run(() => V_書出_バッチ群(p_出力先1, p_出力先2, l_処理済み, l_中断));
            try
            {
                foreach (var l_バッチ in l_読込済み.GetConsumingEnumerable(l_中断.Token))
                {
                    l_総ペア数 += l_バッチ.A_件数;
                    _ = Parallel.For(0, l_バッチ.A_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, i =>
                    {
                        l_バッチ.A_結果群[i] = Get_品質トリム済み(Get_前処理結果(l_バッチ.A_順リード配列群[i], l_バッチ.A_順リード品質群[i], l_バッチ.A_逆リード配列群[i], l_バッチ.A_逆リード品質群[i], p_Phredオフセット), p_Phredオフセット, l_トリム閾値);
                    });
                    for (var i = 0; i < l_バッチ.A_件数; i++)
                    {
                        var l_結果 = l_バッチ.A_結果群[i];
                        if (l_結果.A_Hasアダプタ検出)
                        {
                            l_アダプタ検出ペア数++;
                        }

                        l_訂正塩基数 += l_結果.A_訂正塩基数;
                        l_総塩基数 += l_バッチ.A_順リード配列群[i].Length + l_バッチ.A_逆リード配列群[i].Length;
                        l_トリム塩基数 += l_結果.A_品質トリム塩基数;
                    }

                    l_処理済み.Add(l_バッチ, l_中断.Token);
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
                l_処理済み.CompleteAdding();
            }

            l_読み込み.GetAwaiter().GetResult();
            l_書き出し.GetAwaiter().GetResult();
            return new 前処理統計(l_総ペア数, l_アダプタ検出ペア数, l_訂正塩基数, l_総塩基数, l_トリム塩基数, l_トリム閾値);
        }

        /// <summary>
        /// 前処理の集計をログへ出力する
        /// </summary>
        /// <param name="p_統計">前処理の集計</param>
        public static void V_出力_前処理統計(前処理統計 p_統計)
        {
            Logger.V_出力(メッセージID.前処理統計, p_統計.A_アダプタ検出ペア数, p_統計.A_総ペア数, p_統計.A_訂正塩基数);
            if (p_統計.A_品質トリム閾値 > 0)
            {
                var l_割合 = p_統計.A_総塩基数 > 0L ? 100D * p_統計.A_品質トリム塩基数 / p_統計.A_総塩基数 : 0D;
                Logger.V_出力(メッセージID.品質トリム統計, p_統計.A_品質トリム閾値, p_統計.A_品質トリム塩基数.ToString(C_項目_N0, System.Globalization.CultureInfo.InvariantCulture), l_割合.ToString(C_項目_F1, System.Globalization.CultureInfo.InvariantCulture));
            }
        }

        /// <summary>
        /// ペアの両方の 3' 末端から、品質の崩れた区間を切る
        /// </summary>
        /// <param name="p_結果">重なり処理まで済んだペア</param>
        /// <param name="p_Phredオフセット">クオリティ文字から Phred スコアを引くためのオフセット</param>
        /// <param name="p_閾値">これより低い品質が続く区間を切る (0 なら切らない)</param>
        /// <returns>切った後のペア (切った塩基数を含む)</returns>
        public static ペア前処理結果 Get_品質トリム済み(ペア前処理結果 p_結果, int p_Phredオフセット, int p_閾値)
        {
            if (p_閾値 <= 0)
            {
                return p_結果;
            }

            var l_基準長 = Get_品質トリム後の長さ(p_結果.A_順リード品質, p_Phredオフセット, p_閾値);
            var l_比較長 = Get_品質トリム後の長さ(p_結果.A_逆リード品質, p_Phredオフセット, p_閾値);
            var l_切った数 = p_結果.A_順リード配列.Length - l_基準長 + p_結果.A_逆リード配列.Length - l_比較長;
            return p_結果 with
            {
                A_順リード配列 = p_結果.A_順リード配列[..l_基準長],
                A_順リード品質 = p_結果.A_順リード品質[..l_基準長],
                A_逆リード配列 = p_結果.A_逆リード配列[..l_比較長],
                A_逆リード品質 = p_結果.A_逆リード品質[..l_比較長],
                A_品質トリム塩基数 = p_結果.A_品質トリム塩基数 + l_切った数,
            };
        }

        /// <summary>
        /// 3' 末端から品質の崩れた区間を切った後の長さ
        /// </summary>
        /// <param name="p_クオリティ">クオリティ文字列</param>
        /// <param name="p_Phredオフセット">クオリティ文字から Phred スコアを引くためのオフセット</param>
        /// <param name="p_閾値">品質の閾値</param>
        /// <returns>残す長さ</returns>
        public static int Get_品質トリム後の長さ(string p_クオリティ, int p_Phredオフセット, int p_閾値)
        {
            var l_和 = 0;
            var l_最大 = 0;
            var l_長さ = p_クオリティ.Length;
            for (var i = p_クオリティ.Length - 1; i >= 1; i--)
            {
                l_和 += p_閾値 - (p_クオリティ[i] - p_Phredオフセット);
                if (l_和 < 0)
                {
                    break;
                }

                if (l_和 > l_最大)
                {
                    l_最大 = l_和;
                    l_長さ = i;
                }
            }

            return l_長さ;
        }

        /// <summary>
        /// 1 ペア分の前処理を行う
        /// </summary>
        /// <param name="p_基準配列">read1 の配列</param>
        /// <param name="p_クオリティ1">read1 のクオリティ</param>
        /// <param name="p_比較配列">read2 の配列</param>
        /// <param name="p_クオリティ2">read2 のクオリティ</param>
        /// <param name="p_Phredオフセット">クオリティ文字から Phred スコアを引くためのオフセット</param>
        /// <returns>前処理の結果</returns>
        public static ペア前処理結果 Get_前処理結果(string p_基準配列, string p_クオリティ1, string p_比較配列, string p_クオリティ2, int p_Phredオフセット)
        {
            var l_塩基列1 = Util.V_変換_塩基列(p_基準配列);
            var l_RC比較配列 = Util.V_逆相補_曖昧塩基あり(p_比較配列);
            var l_塩基列2RC = Util.V_変換_塩基列(l_RC比較配列);
            var l_オーバーラップ = Get_最適オーバーラップ(l_塩基列1, l_塩基列2RC, C_最小オーバーラップ長, C_許容不一致率, out _);
            if (l_オーバーラップ is not { } l_重なり)
            {
                return new ペア前処理結果(p_基準配列, p_クオリティ1, p_比較配列, p_クオリティ2, false, 0);
            }

            var l_フラグメント長 = l_重なり.A_offset + p_比較配列.Length;
            var l_開始1 = Math.Max(0, l_重なり.A_offset);
            var l_開始2RC = Math.Max(0, -l_重なり.A_offset);
            var l_基準配列文字 = p_基準配列.ToCharArray();
            var l_比較配列文字 = p_比較配列.ToCharArray();
            var l_訂正数 = 0;
            for (var i = 0; i < l_重なり.A_重なり長; i++)
            {
                var l_位置1 = l_開始1 + i;
                var l_位置2RC = l_開始2RC + i;
                var l_塩基1 = l_塩基列1[l_位置1];
                var l_塩基2RC = l_塩基列2RC[l_位置2RC];
                if (l_塩基1 == l_塩基2RC || l_塩基1 == Consts.無効な塩基 || l_塩基2RC == Consts.無効な塩基)
                {
                    continue;
                }

                var l_位置2 = p_比較配列.Length - 1 - l_位置2RC;
                var l_スコア1 = p_クオリティ1[l_位置1] - p_Phredオフセット;
                var l_スコア2 = p_クオリティ2[l_位置2] - p_Phredオフセット;
                if (l_スコア1 >= C_高信頼スコア && l_スコア2 <= C_低信頼スコア)
                {
                    l_比較配列文字[l_位置2] = Util.Get_相補塩基(l_RC比較配列[l_位置2RC]);
                    l_訂正数++;
                }
                else if (l_スコア2 >= C_高信頼スコア && l_スコア1 <= C_低信頼スコア)
                {
                    l_基準配列文字[l_位置1] = l_RC比較配列[l_位置2RC];
                    l_訂正数++;
                }
            }

            var l_新基準長 = Math.Min(p_基準配列.Length, l_フラグメント長);
            var l_新比較長 = Math.Min(p_比較配列.Length, l_フラグメント長);
            var l_Hasアダプタ検出 = l_新基準長 < p_基準配列.Length || l_新比較長 < p_比較配列.Length;
            return new ペア前処理結果(new string(l_基準配列文字, 0, l_新基準長), p_クオリティ1[..l_新基準長], new string(l_比較配列文字, 0, l_新比較長), p_クオリティ2[..l_新比較長], l_Hasアダプタ検出, l_訂正数);
        }

        /// <summary>
        /// R1 と RC (R2) の最良の重なり位置を探す
        /// </summary>
        /// <param name="p_塩基列1">R1 の塩基 ID 列</param>
        /// <param name="p_塩基列2RC">RC (R2) の塩基 ID 列</param>
        /// <param name="p_最小重なり長">重なりとみなすために要求する最小長</param>
        /// <param name="p_許容不一致率">重なりとみなすために許す不一致率</param>
        /// <param name="p_Has対抗馬">条件を満たすオフセットが 2 つ以上あったか</param>
        /// <param name="p_最小オフセット">試すオフセットの下限、null なら重なりが取れる全域</param>
        /// <param name="p_最大オフセット">試すオフセットの上限、null なら重なりが取れる全域</param>
        /// <returns>条件を満たすもののなかで重なりが最長のもの (同点なら不一致数が少ないもの) <br/> 見つからなければ null (通常の、フラグメント長がリード長を超える場合)</returns>
        internal static オーバーラップ結果? Get_最適オーバーラップ(byte[] p_塩基列1, byte[] p_塩基列2RC, int p_最小重なり長, double p_許容不一致率, out bool p_Has対抗馬, int? p_最小オフセット = null, int? p_最大オフセット = null)
        {
            p_Has対抗馬 = false;
            var l_n1 = p_塩基列1.Length;
            var l_n2 = p_塩基列2RC.Length;
            var l_詰め1 = PackedBases.Get_作る(p_塩基列1);
            var l_詰め2 = PackedBases.Get_作る(p_塩基列2RC);
            オーバーラップ結果? l_最良 = null;
            var l_下限 = Math.Max(-(l_n2 - p_最小重なり長), p_最小オフセット ?? int.MinValue);
            var l_上限 = Math.Min(l_n1 - p_最小重なり長, p_最大オフセット ?? int.MaxValue);
            for (var l_offset = l_下限; l_offset <= l_上限; l_offset++)
            {
                var l_重なり長 = l_offset >= 0 ? Math.Min(l_n1 - l_offset, l_n2) : Math.Min(l_n1, l_n2 + l_offset);
                if (l_重なり長 < p_最小重なり長)
                {
                    continue;
                }

                var l_開始1 = Math.Max(0, l_offset);
                var l_開始2 = Math.Max(0, -l_offset);
                var l_許容不一致数 = (int)(l_重なり長 * p_許容不一致率);
                var l_不一致数 = l_詰め1 is { } l_詰めA && l_詰め2 is { } l_詰めB ? Get_不一致数_語単位(l_詰めA, l_詰めB, l_開始1, l_開始2, l_重なり長, l_許容不一致数) : Get_不一致数_1塩基ずつ(p_塩基列1, p_塩基列2RC, l_開始1, l_開始2, l_重なり長, l_許容不一致数);
                if (l_不一致数 > l_許容不一致数)
                {
                    continue;
                }

                if (l_最良 is not { } l_現在最良)
                {
                    l_最良 = new オーバーラップ結果(l_offset, l_重なり長, l_不一致数);
                    continue;
                }

                p_Has対抗馬 = true;
                if (l_重なり長 > l_現在最良.A_重なり長 || (l_重なり長 == l_現在最良.A_重なり長 && l_不一致数 < l_現在最良.A_不一致数))
                {
                    l_最良 = new オーバーラップ結果(l_offset, l_重なり長, l_不一致数);
                }
            }

            return l_最良;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// ペアをバッチに詰めて順に渡す (どちらかのファイルが尽きたところで止める)
        /// </summary>
        /// <param name="p_順リードのパス">リード 1 のパス</param>
        /// <param name="p_逆リードのパス">リード 2 のパス</param>
        /// <param name="p_渡し先">読み込んだバッチの渡し先</param>
        /// <param name="p_中断">どこかで失敗したときに止める合図</param>
        private static void V_読込_バッチ群(string p_順リードのパス, string p_逆リードのパス, BlockingCollection<前処理バッチ> p_渡し先, CancellationTokenSource p_中断)
        {
            try
            {
                using var l_読み込み1 = new FastqReader(p_順リードのパス);
                using var l_読み込み2 = new FastqReader(p_逆リードのパス);
                while (l_読み込み1.Has続き() && l_読み込み2.Has続き())
                {
                    var l_バッチ = new 前処理バッチ();
                    while (l_バッチ.A_件数 < C_前処理バッチサイズ && l_読み込み1.Has続き() && l_読み込み2.Has続き())
                    {
                        var l_順リード = l_読み込み1.Get_次のリード_軽量();
                        var l_逆リード = l_読み込み2.Get_次のリード_軽量();
                        l_バッチ.A_順リードID群[l_バッチ.A_件数] = l_順リード.A_ID;
                        l_バッチ.A_逆リードID群[l_バッチ.A_件数] = l_逆リード.A_ID;
                        l_バッチ.A_順リード配列群[l_バッチ.A_件数] = l_順リード.A_生リード;
                        l_バッチ.A_逆リード配列群[l_バッチ.A_件数] = l_逆リード.A_生リード;
                        l_バッチ.A_順リード品質群[l_バッチ.A_件数] = l_順リード.A_クオリティ;
                        l_バッチ.A_逆リード品質群[l_バッチ.A_件数] = l_逆リード.A_クオリティ;
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
        /// 前処理を終えたバッチを、受け取った順に書き出す
        /// </summary>
        /// <param name="p_出力先1">前処理したリード 1 の書き出し先</param>
        /// <param name="p_出力先2">前処理したリード 2 の書き出し先</param>
        /// <param name="p_受け取り元">前処理を終えたバッチ</param>
        /// <param name="p_中断">どこかで失敗したときに止める合図</param>
        private static void V_書出_バッチ群(string p_出力先1, string p_出力先2, BlockingCollection<前処理バッチ> p_受け取り元, CancellationTokenSource p_中断)
        {
            try
            {
                using var l_書き込み1 = new FastqWriter(p_出力先1);
                using var l_書き込み2 = new FastqWriter(p_出力先2);
                foreach (var l_バッチ in p_受け取り元.GetConsumingEnumerable(p_中断.Token))
                {
                    for (var i = 0; i < l_バッチ.A_件数; i++)
                    {
                        var l_結果 = l_バッチ.A_結果群[i];
                        l_書き込み1.V_書き込み(l_バッチ.A_順リードID群[i], l_結果.A_順リード配列, l_結果.A_順リード品質);
                        l_書き込み2.V_書き込み(l_バッチ.A_逆リードID群[i], l_結果.A_逆リード配列, l_結果.A_逆リード品質);
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
        /// 重なり区間の不一致数を 32 塩基ずつ数える
        /// </summary>
        /// <param name="p_詰め1">R1 を詰めたもの</param>
        /// <param name="p_詰め2">RC (R2) を詰めたもの</param>
        /// <param name="p_開始1">R1 側で比べ始める位置</param>
        /// <param name="p_開始2">RC (R2) 側で比べ始める位置</param>
        /// <param name="p_重なり長">比べる長さ</param>
        /// <param name="p_許容不一致数">許容する不一致の数</param>
        /// <returns>不一致の数、許容数を超えた時点で打ち切るのでその場合は許容数より大きい値</returns>
        private static int Get_不一致数_語単位(PackedBases p_詰め1, PackedBases p_詰め2, int p_開始1, int p_開始2, int p_重なり長, int p_許容不一致数)
        {
            var l_不一致数 = 0;
            for (var i = 0; i < p_重なり長; i += Consts.ワードあたりの塩基数)
            {
                var l_今回 = Math.Min(Consts.ワードあたりの塩基数, p_重なり長 - i);
                l_不一致数 += PackedBases.Get_不一致数(p_詰め1.Get_窓(p_開始1 + i), p_詰め2.Get_窓(p_開始2 + i), l_今回);
                if (l_不一致数 > p_許容不一致数)
                {
                    return l_不一致数;
                }
            }

            return l_不一致数;
        }

        /// <summary>
        /// 曖昧塩基を含んで詰められない場合に、1 塩基ずつ不一致の数を数える
        /// </summary>
        /// <param name="p_塩基列1">R1 の塩基 ID 列</param>
        /// <param name="p_塩基列2RC">RC (R2) の塩基 ID 列</param>
        /// <param name="p_開始1">R1 側で比べ始める位置</param>
        /// <param name="p_開始2">RC (R2) 側で比べ始める位置</param>
        /// <param name="p_重なり長">比べる長さ</param>
        /// <param name="p_許容不一致数">許容する不一致の数</param>
        /// <returns>不一致の数、許容数を超えた時点で打ち切るのでその場合は許容数より大きい値</returns>
        private static int Get_不一致数_1塩基ずつ(byte[] p_塩基列1, byte[] p_塩基列2RC, int p_開始1, int p_開始2, int p_重なり長, int p_許容不一致数)
        {
            var l_不一致数 = 0;
            for (var i = 0; i < p_重なり長; i++)
            {
                if (p_塩基列1[p_開始1 + i] != p_塩基列2RC[p_開始2 + i])
                {
                    l_不一致数++;
                    if (l_不一致数 > p_許容不一致数)
                    {
                        return l_不一致数;
                    }
                }
            }

            return l_不一致数;
        }

        #endregion
    }
}
