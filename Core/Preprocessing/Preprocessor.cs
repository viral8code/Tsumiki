using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// ペアエンドの2本を重ね合わせて(R1 と RC(R2))、アダプタリードスルーの
    /// トリムと、高信頼/低信頼が明確に分かれる位置でのペア相互訂正を行う。
    /// k-mer カウントより前段の、fastp 型の前処理層(ErrorCorrector のさらに前段)。
    ///
    /// ErrorCorrector(k-mer スペクトルに基づく訂正)とは独立な証拠源
    /// ―― 同じ断片を両端から2回読んだという事実そのもの ―― を使うため、
    /// 両方を通して初めて捕まえられる誤りがある。
    /// アダプタ配列は全リードで共通のため、カットオフでは絶対に落ちない
    /// 高カバレッジな偽k-merとしてグラフに定着する。ここで先に取り除く。
    /// </summary>
    internal static class Preprocessor
    {
        /// <summary>
        /// ペアの重なりとみなすために要求する最小長。これより短い一致は
        /// 偶然の一致と区別できない。
        /// </summary>
        private const int 最小オーバーラップ長 = 30;

        /// <summary>
        /// この不一致率までは同一断片から重なって読んだものとみなす。
        /// シーケンシングエラー由来の不一致を許容しつつ、無関係な配列同士の
        /// 偶然の一致を弾くための閾値。
        /// </summary>
        private const double 許容不一致率 = 0.2;

        /// <summary>相互訂正で「高信頼」とみなす最小Phredスコア。</summary>
        private const int 高信頼スコア = 30;

        /// <summary>
        /// 相互訂正で「低信頼」とみなす最大Phredスコア。相方が高信頼スコア以上の
        /// ときに限り、この閾値以下の側だけを上書きする。片方だけが明確に
        /// 高品質という非対称なケースに限定することで、誤訂正を避ける。
        /// </summary>
        private const int 低信頼スコア = 14;

        /// <summary>
        /// 1バッチあたりのペア数。ErrorCorrector と同じ理由(出力の行順を
        /// 入力と厳密に一致させる必要がある)で、まとめて読み並列に処理し
        /// 順番通りに書く形にしている。
        /// </summary>
        private const int 前処理バッチサイズ = 20000;

        /// <summary>
        /// ペアの FASTQ を読み込んで前処理し、結果を出力先へ書き出す。
        /// </summary>
        public static 前処理統計 V_前処理_リードファイル(
            string p_リード1のパス, string p_リード2のパス, string p_出力先1, string p_出力先2)
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

            while (l_読み込み1.Get_続きがあるか() && l_読み込み2.Get_続きがあるか())
            {
                var l_件数 = 0;
                while (l_件数 < 前処理バッチサイズ
                    && l_読み込み1.Get_続きがあるか() && l_読み込み2.Get_続きがあるか())
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
                    l_結果群[i] = Get_前処理結果(
                        l_配列1群[i], l_クオリティ1群[i], l_配列2群[i], l_クオリティ2群[i], l_Phredオフセット);
                });

                for (var i = 0; i < l_件数; i++)
                {
                    var l_結果 = l_結果群[i];
                    if (l_結果.A_アダプタを検出したか)
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

        public static void V_出力_前処理統計(前処理統計 p_統計)
        {
            Logger.V_出力(メッセージID.前処理統計, p_統計.A_アダプタ検出ペア数, p_統計.A_総ペア数, p_統計.A_訂正塩基数);
        }

        /// <summary>
        /// 1ペア分の前処理。副作用のない純粋関数。
        /// R1 と RC(R2) を重ね合わせ、重なりが見つかった場合のみ
        /// (1) アダプタ読み過ごし分のトリムと (2) 重なり領域内の相互訂正を行う。
        /// 重なりが見つからない(=フラグメント長がリード長を上回る、通常の)場合は
        /// 元のリードをそのまま返す。
        /// </summary>
        internal static ペア前処理結果 Get_前処理結果(
            string p_配列1, string p_クオリティ1, string p_配列2, string p_クオリティ2, int p_Phredオフセット)
        {
            var l_塩基列1 = Util.V_変換_塩基列(p_配列1);
            var l_RC配列2 = Util.V_逆相補_曖昧塩基あり(p_配列2);
            var l_塩基列2RC = Util.V_変換_塩基列(l_RC配列2);

            var l_オーバーラップ = Get_最適オーバーラップ(l_塩基列1, l_塩基列2RC);
            if (l_オーバーラップ is not { } l_重なり)
            {
                return new ペア前処理結果(p_配列1, p_クオリティ1, p_配列2, p_クオリティ2, false, 0);
            }

            // R1 は断片の先頭から、RC(R2) は断片上の (offset) 位置から始まる
            // ([[オーバーラップ結果]] 参照)。したがってフラグメント長は
            // offset + R2 の長さ で求まる。
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
                    // 一致している、または片方が曖昧塩基で真の塩基が確認できない。
                    // 曖昧塩基の位置は書き換えない(ErrorCorrector と同じ方針)。
                    continue;
                }

                var l_位置2 = p_配列2.Length - 1 - l_位置2RC;
                var l_スコア1 = p_クオリティ1[l_位置1] - p_Phredオフセット;
                var l_スコア2 = p_クオリティ2[l_位置2] - p_Phredオフセット;

                if (l_スコア1 >= 高信頼スコア && l_スコア2 <= 低信頼スコア)
                {
                    // RC(R2) 側は既に R1 と同じ向きの塩基になっているので、そのまま写す。
                    l_配列2文字[l_位置2] = Get_相補文字(l_RC配列2[l_位置2RC]);
                    l_訂正数++;
                }
                else if (l_スコア2 >= 高信頼スコア && l_スコア1 <= 低信頼スコア)
                {
                    l_配列1文字[l_位置1] = l_RC配列2[l_位置2RC];
                    l_訂正数++;
                }
                // どちらも高信頼側が無い(質の差が明確でない)場合は、
                // fastp と同じ方針でどちらも訂正しない。
            }

            var l_新長さ1 = Math.Min(p_配列1.Length, l_フラグメント長);
            var l_新長さ2 = Math.Min(p_配列2.Length, l_フラグメント長);
            var l_アダプタを検出したか = l_新長さ1 < p_配列1.Length || l_新長さ2 < p_配列2.Length;

            return new ペア前処理結果(
                new string(l_配列1文字, 0, l_新長さ1),
                p_クオリティ1[..l_新長さ1],
                new string(l_配列2文字, 0, l_新長さ2),
                p_クオリティ2[..l_新長さ2],
                l_アダプタを検出したか,
                l_訂正数);
        }

        /// <summary>
        /// R1 と RC(R2) の最良の重なり位置を探す。全オフセットのうち、
        /// 最小オーバーラップ長以上・不一致率が閾値以下のもののなかで、
        /// 重なりが最長のもの(同点なら不一致数が少ないもの)を返す。
        /// 見つからなければ null(=通常の、フラグメント長がリード長を超える場合)。
        /// </summary>
        private static オーバーラップ結果? Get_最適オーバーラップ(byte[] p_塩基列1, byte[] p_塩基列2RC)
        {
            var l_n1 = p_塩基列1.Length;
            var l_n2 = p_塩基列2RC.Length;

            オーバーラップ結果? l_最良 = null;
            for (var l_offset = -(l_n2 - 最小オーバーラップ長); l_offset <= l_n1 - 最小オーバーラップ長; l_offset++)
            {
                var l_重なり長 = l_offset >= 0 ? Math.Min(l_n1 - l_offset, l_n2) : Math.Min(l_n1, l_n2 + l_offset);
                if (l_重なり長 < 最小オーバーラップ長)
                {
                    continue;
                }

                var l_開始1 = Math.Max(0, l_offset);
                var l_開始2 = Math.Max(0, -l_offset);
                var l_許容不一致数 = (int)(l_重なり長 * 許容不一致率);

                var l_不一致数 = 0;
                for (var i = 0; i < l_重なり長; i++)
                {
                    if (p_塩基列1[l_開始1 + i] != p_塩基列2RC[l_開始2 + i])
                    {
                        l_不一致数++;
                        if (l_不一致数 > l_許容不一致数)
                        {
                            // 早期打ち切り: 大半のオフセットは無関係な配列同士の比較になるため、
                            // 閾値を超えた時点でやめないと全ペア×全オフセットが O(read長) のままになる。
                            break;
                        }
                    }
                }
                if (l_不一致数 > l_許容不一致数)
                {
                    continue;
                }

                if (l_最良 is not { } l_現在最良
                    || l_重なり長 > l_現在最良.A_重なり長
                    || (l_重なり長 == l_現在最良.A_重なり長 && l_不一致数 < l_現在最良.A_不一致数))
                {
                    l_最良 = new オーバーラップ結果(l_offset, l_重なり長, l_不一致数);
                }
            }

            return l_最良;
        }

        private static char Get_相補文字(char p_塩基)
        {
            return p_塩基 switch
            {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                var l_文字 => l_文字,
            };
        }
    }
}
