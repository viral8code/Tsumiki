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
    /// <remarks>
    /// 信頼できない k-mer 窓を最も多く信頼状態へ変える 1 塩基置換を貪欲に選び、改善が見込めなくなるまで反復する<br/>
    /// 曖昧塩基の位置は書き換えず、それを含む窓は評価からも除外する
    /// </remarks>
    internal static class ErrorCorrector
    {
        #region 定数

        /// <summary>
        /// 1 バッチあたりのリード数
        /// </summary>
        /// <remarks>
        /// 訂正自体は独立に並列化できるが、出力の行順はペアの対応付けを保つため入力と厳密に一致させる必要があり、「まとめて読む → 並列に訂正 → 順番通りに書く」形にしている
        /// </remarks>
        private const int 訂正バッチサイズ = 20_000;

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
        /// <remarks>
        /// 「信頼できる k-mer」の判定には、本アセンブリと同じ -kc のカットオフ値を使って構築した専用の k-mer インデックス (このメソッド内で完結し、本パイプライン用のインデックスとは独立) を用いる
        /// </remarks>
        public static void V_訂正_リードファイル(string p_リード1のパス, string? p_リード2のパス, string p_一時ディレクトリ, string p_出力先1, string? p_出力先2)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            var l_訂正用一時ディレクトリ = Path.Combine(p_一時ディレクトリ, "error_correction");
            _ = Directory.CreateDirectory(l_訂正用一時ディレクトリ);

            Logger.V_出力(メッセージID.エラー訂正_スペクトル構築);
            using (var l_kmerインデックス = new TrustedKmerIndex(l_訂正用一時ディレクトリ))
            {
                KmerCounting.V_読込_リードファイル(p_リード1のパス, l_kmerインデックス);
                if (p_リード2のパス != null)
                {
                    KmerCounting.V_読込_リードファイル(p_リード2のパス, l_kmerインデックス);
                }

                // 訂正の判定はこのカットオフが全て
                // 既定値のままだと
                // エラー由来の k-mer まで信頼扱いになり、訂正が起きない
                KmerCutoffSelector.V_解決_kmerカットオフ(ConfigurationManager.A_実行時引数, l_kmerインデックス);
                _ = l_kmerインデックス.V_カットオフ(ConfigurationManager.A_実行時引数.A_kmerカットオフ);

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
        /// 1 リード (塩基 ID 空間のバイト列、曖昧塩基は Consts.無効な塩基) を貪欲法で訂正する
        /// </summary>
        /// <param name="p_リード"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_最大反復数"></param>
        /// <remarks>
        /// 副作用のない純粋関数 (入力は変更しない) <br/>
        /// k が 64 以下なら、窓を 2 bit パックして転がす経路を使う<br/>
        /// 判定内容も選ぶ置換も逐次経路と同じで、1 窓あたりの手間だけが O (k) から O (1) に変わる (逐次経路は窓を見るたびにパックと逆相補を取り直していた)
        /// </remarks>
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
        /// <remarks>
        /// パックできないので窓ごとに評価する
        /// </remarks>
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
                    // これ以上、信頼できる窓を純増させる置換が見つからない
                    // (=残った信頼できない窓は、単発の置換では解決できない)
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

            // 訂正処理は副作用のない純粋関数で、k-mer インデックスも
            // カットオフ後は読み取り専用なので、リード単位で安全に並列化できる
            // 実データ (35 x, 800 k ペア) で単一スレッドだと 40 分以上かかっており、
            // パイプライン全体の律速になっていた
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            using var l_読み込み = new FastqReader(p_入力パス);
            using var l_書き込み = new FastqWriter(p_出力パス);

            var l_ID群 = new string[訂正バッチサイズ];
            var l_クオリティ群 = new string[訂正バッチサイズ];
            var l_塩基列群 = new byte[訂正バッチサイズ][];
            var l_結果群 = new 訂正結果[訂正バッチサイズ];

            while (l_読み込み.Has続き())
            {
                var l_件数 = 0;
                while (l_件数 < 訂正バッチサイズ && l_読み込み.Has続き())
                {
                    var l_リード = l_読み込み.Get_次のリード_軽量();
                    l_ID群[l_件数] = l_リード.A_ID;
                    l_クオリティ群[l_件数] = l_リード.A_クオリティ;
                    l_塩基列群[l_件数] = l_リード.A_塩基列!;
                    l_件数++;
                }
                l_総リード数 += l_件数;

                _ = Parallel.For(0, l_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, i =>
                {
                    l_結果群[i] = Get_訂正結果(l_塩基列群[i], p_kmerインデックス, p_k長);
                });

                for (var i = 0; i < l_件数; i++)
                {
                    var l_結果 = l_結果群[i];
                    if (l_結果.A_訂正数 > 0)
                    {
                        l_訂正されたリード数++;
                        l_総訂正塩基数 += l_結果.A_訂正数;
                    }
                    l_書き込み.V_書き込み(l_ID群[i], string.Join(string.Empty, l_結果.A_塩基列.Select(Util.V_変換_塩基文字)), l_クオリティ群[i]);
                }
            }

            return new ファイル訂正統計(l_総リード数, l_訂正されたリード数, l_総訂正塩基数);
        }

        /// <summary>
        /// パック経路の訂正本体
        /// </summary>
        /// <param name="p_塩基列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_最大反復数"></param>
        /// <remarks>
        /// 逐次経路と同じ貪欲法で、窓の評価だけをパック値の更新で済ませる
        /// </remarks>
        /// <returns></returns>
        private static 訂正結果 Get_訂正結果_パック(byte[] p_塩基列, TrustedKmerIndex p_kmerインデックス, int p_k長, int p_最大反復数)
        {
            var l_窓数 = p_塩基列.Length - p_k長 + 1;
            var l_パック = new UInt128[l_窓数];
            var l_逆相補 = new UInt128[l_窓数];
            var l_無効数 = new int[l_窓数];
            var l_信頼状況 = new bool[l_窓数];

            // 「位置 p を含む窓の中に信頼できないものがあるか」と、候補を
            // 打ち切ってよいかの上界計算に使う、信頼できない窓の累積数
            var l_未信頼累積 = new int[l_窓数 + 1];
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
        /// <remarks>
        /// 曖昧塩基はコドン 0 として詰めておき、判定では曖昧塩基の数で弾く (窓から出れば残りのコドンはそのまま正しい)
        /// </remarks>
        private static void V_計算_窓状態(byte[] p_塩基列, int p_k長, TrustedKmerIndex p_kmerインデックス, UInt128[] p_パック, UInt128[] p_逆相補, int[] p_無効数, bool[] p_信頼状況)
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
        /// <remarks>
        /// 置換で変わるのは各窓のうち 1 コドンだけなので、窓ごとにパック値を詰め直さず、その 1 コドンを差し替えて引く<br/>
        /// 残りの窓が全て改善に転じても現在の最良に届かないと分かった時点で打ち切る (打ち切っても選ばれる置換は変わらない<br/>
        /// 改善数が同じ候補は元から採用されない)
        /// </remarks>
        /// <returns></returns>
        private static int Get_置換改善数_パック(int p_位置, byte p_候補, int p_窓開始, int p_窓終了, int p_k長, TrustedKmerIndex p_kmerインデックス, UInt128[] p_パック, UInt128[] p_逆相補, int[] p_無効数, bool[] p_信頼状況, int[] p_未信頼累積, int p_最良改善数)
        {
            var l_コドン = Get_コドン(p_候補);
            var l_相補コドン = Get_相補コドン(p_候補);
            var l_改善数 = 0;

            for (var w = p_窓開始; w <= p_窓終了; w++)
            {
                // まだ見ていない窓が全て改善しても最良に届かないなら、
                // これ以上引く意味がない
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
        /// <remarks>
        /// 曖昧塩基は 0 として詰める (判定は無効数で弾く)
        /// </remarks>
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
        /// <remarks>
        /// 塩基列は評価後、呼び出し前の状態に戻す (副作用を残さない)
        /// </remarks>
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
