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
    /// <remarks>
    /// k-mer カウントより前段の、fastp 型の前処理層になる (ErrorCorrector のさらに前段) <br/>
    /// ErrorCorrector (k-mer スペクトルに基づく訂正) とは独立な証拠源、つまり同じ断片を両端から 2 回読んだという事実そのものを使うため、両方を通して初めて捕まえられる誤りがある<br/>
    /// アダプタ配列は全リードで共通のため、カットオフでは絶対に落ちない高カバレッジな偽 k-mer としてグラフに定着するので、ここで先に取り除く
    /// </remarks>
    internal static class Preprocessor
    {
        #region 定数

        /// <summary>
        /// ペアの重なりとみなすために要求する最小長
        /// </summary>
        /// <remarks>
        /// これより短い一致は偶然の一致と区別できない
        /// </remarks>
        private const int 最小オーバーラップ長 = 30;

        /// <summary>
        /// この不一致率までは同一断片から重なって読んだものとみなす
        /// </summary>
        /// <remarks>
        /// シーケンシングエラー由来の不一致を許容しつつ、無関係な配列同士の偶然の一致を弾くための閾値になる
        /// </remarks>
        private const double 許容不一致率 = 0.2D;

        /// <summary>
        /// 相互訂正で高信頼とみなす最小 Phred スコア
        /// </summary>
        private const int 高信頼スコア = 30;

        /// <summary>
        /// 相互訂正で低信頼とみなす最大 Phred スコア
        /// </summary>
        /// <remarks>
        /// 相方が高信頼スコア以上のときに限り、この閾値以下の側だけを上書きする<br/>
        /// 片方だけが明確に高品質という非対称なケースに限定することで、誤訂正を避ける
        /// </remarks>
        private const int 低信頼スコア = 14;

        /// <summary>
        /// 1 バッチあたりのペア数
        /// </summary>
        /// <remarks>
        /// ErrorCorrector と同じ理由 (出力の行順を入力と厳密に一致させる必要がある) で、まとめて読み並列に処理し順番通りに書く形にしている
        /// </remarks>
        private const int 前処理バッチサイズ = 20_000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// ペアの FASTQ を読み込んで前処理し、結果を出力先へ書き出す
        /// </summary>
        /// <param name="p_リード1のパス">リード 1 のパス</param>
        /// <param name="p_リード2のパス">リード 2 のパス</param>
        /// <param name="p_出力先1">前処理したリード 1 の書き出し先</param>
        /// <param name="p_出力先2">前処理したリード 2 の書き出し先</param>
        /// <returns>前処理の集計</returns>
        public static 前処理統計 V_前処理_リードファイル(string p_リード1のパス, string p_リード2のパス, string p_出力先1, string p_出力先2)
        {
            var l_Phredオフセット = ConfigurationManager.A_実行時引数.A_Phredオフセット;
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            var l_総ペア数 = 0;
            var l_アダプタ検出ペア数 = 0;
            var l_訂正塩基数 = 0;

            using var l_読み込み1 = new FastqReader(p_リード1のパス);
            using var l_読み込み2 = new FastqReader(p_リード2のパス);
            using var l_書き込み1 = new FastqWriter(p_出力先1);
            using var l_書き込み2 = new FastqWriter(p_出力先2);

            var l_ID1群 = new string[前処理バッチサイズ];
            var l_ID2群 = new string[前処理バッチサイズ];
            var l_配列1群 = new string[前処理バッチサイズ];
            var l_配列2群 = new string[前処理バッチサイズ];
            var l_クオリティ1群 = new string[前処理バッチサイズ];
            var l_クオリティ2群 = new string[前処理バッチサイズ];
            var l_結果群 = new ペア前処理結果[前処理バッチサイズ];

            while (l_読み込み1.Has続き() && l_読み込み2.Has続き())
            {
                var l_件数 = 0;
                while (l_件数 < 前処理バッチサイズ && l_読み込み1.Has続き() && l_読み込み2.Has続き())
                {
                    var l_リード1 = l_読み込み1.Get_次のリード_軽量();
                    var l_リード2 = l_読み込み2.Get_次のリード_軽量();
                    l_ID1群[l_件数] = l_リード1.A_ID;
                    l_ID2群[l_件数] = l_リード2.A_ID;
                    l_配列1群[l_件数] = l_リード1.A_生リード;
                    l_配列2群[l_件数] = l_リード2.A_生リード;
                    l_クオリティ1群[l_件数] = l_リード1.A_クオリティ;
                    l_クオリティ2群[l_件数] = l_リード2.A_クオリティ;
                    l_件数++;
                }
                l_総ペア数 += l_件数;

                _ = Parallel.For(0, l_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, i =>
                {
                    l_結果群[i] = Get_前処理結果(l_配列1群[i], l_クオリティ1群[i], l_配列2群[i], l_クオリティ2群[i], l_Phredオフセット);
                });

                for (var i = 0; i < l_件数; i++)
                {
                    var l_結果 = l_結果群[i];
                    if (l_結果.A_Hasアダプタ検出)
                    {
                        l_アダプタ検出ペア数++;
                    }
                    l_訂正塩基数 += l_結果.A_訂正塩基数;
                    l_書き込み1.V_書き込み(l_ID1群[i], l_結果.A_配列1, l_結果.A_クオリティ1);
                    l_書き込み2.V_書き込み(l_ID2群[i], l_結果.A_配列2, l_結果.A_クオリティ2);
                }
            }

            return new 前処理統計(l_総ペア数, l_アダプタ検出ペア数, l_訂正塩基数);
        }

        /// <summary>
        /// 前処理の集計をログへ出力する
        /// </summary>
        /// <param name="p_統計">前処理の集計</param>
        public static void V_出力_前処理統計(前処理統計 p_統計)
        {
            Logger.V_出力(メッセージID.前処理統計, p_統計.A_アダプタ検出ペア数, p_統計.A_総ペア数, p_統計.A_訂正塩基数);
        }

        /// <summary>
        /// 1 ペア分の前処理を行う
        /// </summary>
        /// <remarks>
        /// 副作用のない純粋関数<br/>
        /// R1 と RC (R2) を重ね合わせ、重なりが見つかった場合のみアダプタ読み過ごし分のトリムと、重なり領域内の相互訂正を行う<br/>
        /// 重なりが見つからない (フラグメント長がリード長を上回る、通常の) 場合は元のリードをそのまま返す
        /// </remarks>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_クオリティ1">read1 のクオリティ</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_クオリティ2">read2 のクオリティ</param>
        /// <param name="p_Phredオフセット">クオリティ文字から Phred スコアを引くためのオフセット</param>
        /// <returns>前処理の結果</returns>
        internal static ペア前処理結果 Get_前処理結果(string p_配列1, string p_クオリティ1, string p_配列2, string p_クオリティ2, int p_Phredオフセット)
        {
            var l_塩基列1 = Util.V_変換_塩基列(p_配列1);
            var l_RC配列2 = Util.V_逆相補_曖昧塩基あり(p_配列2);
            var l_塩基列2RC = Util.V_変換_塩基列(l_RC配列2);

            var l_オーバーラップ = Get_最適オーバーラップ(l_塩基列1, l_塩基列2RC, 最小オーバーラップ長, 許容不一致率, out _);
            if (l_オーバーラップ is not { } l_重なり)
            {
                return new ペア前処理結果(p_配列1, p_クオリティ1, p_配列2, p_クオリティ2, false, 0);
            }

            // R1 は断片の先頭から、RC (R2) は断片上の (offset) 位置から始まる
            // ([[オーバーラップ結果]] 参照) したがってフラグメント長は offset + R2 の長さ で求まる
            var l_フラグメント長 = l_重なり.A_offset + p_配列2.Length;
            var l_開始1 = Math.Max(0, l_重なり.A_offset);
            var l_開始2RC = Math.Max(0, -l_重なり.A_offset);

            var l_配列1文字 = p_配列1.ToCharArray();
            var l_配列2文字 = p_配列2.ToCharArray();
            var l_訂正数 = 0;

            for (var i = 0; i < l_重なり.A_重なり長; i++)
            {
                var l_位置1 = l_開始1 + i;
                var l_位置2RC = l_開始2RC + i;
                var l_塩基1 = l_塩基列1[l_位置1];
                var l_塩基2RC = l_塩基列2RC[l_位置2RC];

                if (l_塩基1 == l_塩基2RC
                    || l_塩基1 == Consts.無効な塩基 || l_塩基2RC == Consts.無効な塩基)
                {
                    // 一致しているか、片方が曖昧塩基で真の塩基が確認できない
                    // 曖昧塩基の位置は書き換えない (ErrorCorrector と同じ方針)
                    continue;
                }

                var l_位置2 = p_配列2.Length - 1 - l_位置2RC;
                var l_スコア1 = p_クオリティ1[l_位置1] - p_Phredオフセット;
                var l_スコア2 = p_クオリティ2[l_位置2] - p_Phredオフセット;

                if (l_スコア1 >= 高信頼スコア && l_スコア2 <= 低信頼スコア)
                {
                    // RC (R2) 側は既に R1 と同じ向きの塩基になっているので、そのまま写す
                    l_配列2文字[l_位置2] = Util.Get_相補塩基(l_RC配列2[l_位置2RC]);
                    l_訂正数++;
                }
                else if (l_スコア2 >= 高信頼スコア && l_スコア1 <= 低信頼スコア)
                {
                    l_配列1文字[l_位置1] = l_RC配列2[l_位置2RC];
                    l_訂正数++;
                }
                // どちらも高信頼側が無い (質の差が明確でない) 場合は、fastp と同じ方針でどちらも訂正しない
            }

            var l_新長さ1 = Math.Min(p_配列1.Length, l_フラグメント長);
            var l_新長さ2 = Math.Min(p_配列2.Length, l_フラグメント長);
            var l_Hasアダプタ検出 = l_新長さ1 < p_配列1.Length || l_新長さ2 < p_配列2.Length;

            return new ペア前処理結果(new string(l_配列1文字, 0, l_新長さ1), p_クオリティ1[..l_新長さ1], new string(l_配列2文字, 0, l_新長さ2), p_クオリティ2[..l_新長さ2], l_Hasアダプタ検出, l_訂正数);
        }

        /// <summary>
        /// R1 と RC (R2) の最良の重なり位置を探す
        /// </summary>
        /// <remarks>
        /// 反復配列の中では周期のぶんだけずれた位置が同じくらい良く合うため、最良を 1 つ選ぶだけでは取り違えに気づけない<br/>
        /// 重ねた結果を 1 本の配列として下流へ渡す用途では、この曖昧さを見て捨てる必要がある
        /// </remarks>
        /// <param name="p_塩基列1">R1 の塩基 ID 列</param>
        /// <param name="p_塩基列2RC">RC (R2) の塩基 ID 列</param>
        /// <param name="p_最小重なり長">重なりとみなすために要求する最小長</param>
        /// <param name="p_許容不一致率">重なりとみなすために許す不一致率</param>
        /// <param name="p_Has対抗馬">条件を満たすオフセットが 2 つ以上あったか</param>
        /// <returns>
        /// 条件を満たすもののなかで重なりが最長のもの (同点なら不一致数が少ないもの) <br/>
        /// 見つからなければ null (通常の、フラグメント長がリード長を超える場合)
        /// </returns>
        internal static オーバーラップ結果? Get_最適オーバーラップ(byte[] p_塩基列1, byte[] p_塩基列2RC, int p_最小重なり長, double p_許容不一致率, out bool p_Has対抗馬)
        {
            p_Has対抗馬 = false;
            var l_n1 = p_塩基列1.Length;
            var l_n2 = p_塩基列2RC.Length;

            // 語単位で比べるための詰め直し
            // 曖昧塩基を含む配列は詰められないので、その場合だけ 1 塩基ずつ比べる経路へ落ちる
            var l_詰め1 = PackedBases.Get_作る(p_塩基列1);
            var l_詰め2 = PackedBases.Get_作る(p_塩基列2RC);

            オーバーラップ結果? l_最良 = null;
            for (var l_offset = -(l_n2 - p_最小重なり長); l_offset <= l_n1 - p_最小重なり長; l_offset++)
            {
                var l_重なり長 = l_offset >= 0 ? Math.Min(l_n1 - l_offset, l_n2) : Math.Min(l_n1, l_n2 + l_offset);
                if (l_重なり長 < p_最小重なり長)
                {
                    continue;
                }

                var l_開始1 = Math.Max(0, l_offset);
                var l_開始2 = Math.Max(0, -l_offset);
                var l_許容不一致数 = (int)(l_重なり長 * p_許容不一致率);

                // 早期打ち切り: 大半のオフセットは無関係な配列同士の比較になるため、
                // 閾値を超えた時点でやめないと全ペア × 全オフセットが O (read 長) のままになる
                var l_不一致数 = l_詰め1 is { } l_詰めA && l_詰め2 is { } l_詰めB
                    ? Get_不一致数_語単位(l_詰めA, l_詰めB, l_開始1, l_開始2, l_重なり長, l_許容不一致数)
                    : Get_不一致数_1塩基ずつ(p_塩基列1, p_塩基列2RC, l_開始1, l_開始2, l_重なり長, l_許容不一致数);
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
                if (l_重なり長 > l_現在最良.A_重なり長
                    || (l_重なり長 == l_現在最良.A_重なり長 && l_不一致数 < l_現在最良.A_不一致数))
                {
                    l_最良 = new オーバーラップ結果(l_offset, l_重なり長, l_不一致数);
                }
            }

            return l_最良;
        }

        #endregion

        #region 内部メソッド

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
            for (var i = 0; i < p_重なり長; i += Consts.語あたりの塩基数)
            {
                var l_今回 = Math.Min(Consts.語あたりの塩基数, p_重なり長 - i);
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
        /// <remarks>
        /// 曖昧塩基同士は一致として扱う (2 bit に落とすとこの区別ができないため、詰められる場合と結果が変わりうる)
        /// </remarks>
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
