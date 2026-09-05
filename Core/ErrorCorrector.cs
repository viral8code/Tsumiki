using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// k-merスペクトラムに基づく、Quake/BayesHammer類似の簡易リードエラー訂正。
    ///
    /// 考え方: あるk-merが十分な回数(カットオフ以上)出現していれば「信頼できる」
    /// (真のゲノム配列由来である可能性が高い)とみなし、逆に出現回数が少ない
    /// k-merは配列決定エラーに由来する可能性が高いとみなす。1本のリード中に
    /// 信頼できないk-merが含まれる場合、その周辺の1塩基を置換することで
    /// 信頼できるk-merに変わるかどうかを試し、最も多くの「信頼できないk-mer窓」を
    /// 信頼できる状態に変える置換を貪欲に選んで適用する(1リードにつき複数回、
    /// 改善が見込めなくなるまで反復する)。
    ///
    /// 曖昧塩基(N等)を含む位置は書き換えない(置換候補にも含めない)。
    /// そのような位置を含む窓は常に「信頼できない」扱いとし、
    /// 評価の対象からも除外する。
    /// </summary>
    internal static class ErrorCorrector
    {
        /// <summary>
        /// 1バッチあたりのリード数。訂正はリードごとに独立なので並列化できるが、
        /// 出力の行順はペアエンドの対応付け(read1のn番目とread2のn番目が同じ
        /// フラグメント)を保つため入力と厳密に一致させる必要がある。
        /// そのため「まとめて読む → 並列に訂正 → 順番通りに書く」形にする。
        /// バッチサイズはメモリ使用量(1リードあたり数百バイト)と
        /// 並列化の粒度のバランスで決めた値。
        /// </summary>
        private const int 訂正バッチサイズ = 20000;

        /// <summary>
        /// リードファイルを読み込んでエラー訂正を行い、結果を出力先へ書き出す。
        /// 「信頼できるk-mer」の判定には、本アセンブリと同じ -kc のカットオフ値を
        /// 使って構築した専用の k-mer インデックス(このメソッド内で完結し、
        /// 本パイプライン用のインデックスとは独立)を用いる。
        /// </summary>
        public static void V_訂正_リードファイル(
            string p_リード1のパス, string? p_リード2のパス, string p_一時ディレクトリ,
            string p_出力先1, string? p_出力先2)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            var l_訂正用一時ディレクトリ = Path.Combine(p_一時ディレクトリ, "error_correction");
            _ = Directory.CreateDirectory(l_訂正用一時ディレクトリ);

            Console.WriteLine("[ErrorCorrection] Building k-mer spectrum...");
            using (var l_kmerインデックス = new TrustedKmerIndex(l_訂正用一時ディレクトリ))
            {
                KmerCounting.V_読込_リードファイル(p_リード1のパス, l_kmerインデックス);
                if (p_リード2のパス != null)
                {
                    KmerCounting.V_読込_リードファイル(p_リード2のパス, l_kmerインデックス);
                }
                // 訂正の判定はこのカットオフが全て。既定値のままだと
                // エラー由来の k-mer まで信頼扱いになり、訂正が起きない。
                KmerCutoffSelector.V_解決_kmerカットオフ(
                    ConfigurationManager.A_実行時引数, l_kmerインデックス);
                _ = l_kmerインデックス.V_カットオフ(ConfigurationManager.A_実行時引数.A_kmerカットオフ);

                Console.WriteLine("[ErrorCorrection] Correcting reads...");
                var l_統計1 = Get_訂正統計_ファイル(p_リード1のパス, p_出力先1, l_kmerインデックス, l_k長);
                Console.WriteLine($"[ErrorCorrection] {Path.GetFileName(p_リード1のパス)}: " +
                    $"{l_統計1.A_訂正されたリード数}/{l_統計1.A_総リード数} reads corrected ({l_統計1.A_総訂正塩基数} base corrections total).");

                if (p_リード2のパス != null && p_出力先2 != null)
                {
                    var l_統計2 = Get_訂正統計_ファイル(p_リード2のパス, p_出力先2, l_kmerインデックス, l_k長);
                    Console.WriteLine($"[ErrorCorrection] {Path.GetFileName(p_リード2のパス)}: " +
                        $"{l_統計2.A_訂正されたリード数}/{l_統計2.A_総リード数} reads corrected ({l_統計2.A_総訂正塩基数} base corrections total).");
                }
            }

            Directory.Delete(l_訂正用一時ディレクトリ, recursive: true);
        }

        private static ファイル訂正統計 Get_訂正統計_ファイル(
            string p_入力パス, string p_出力パス, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_総リード数 = 0;
            var l_訂正されたリード数 = 0;
            var l_総訂正塩基数 = 0;

            // 訂正処理は副作用のない純粋関数で、k-mer インデックスも
            // カットオフ後は読み取り専用なので、リード単位で安全に並列化できる。
            // 実データ(35x, 800k ペア)で単一スレッドだと 40 分以上かかっており、
            // パイプライン全体の律速になっていた。
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            using var l_読み込み = new FastqReader(p_入力パス);
            using var l_書き込み = new FastqWriter(p_出力パス);

            var l_ID群 = new string[訂正バッチサイズ];
            var l_クオリティ群 = new string[訂正バッチサイズ];
            var l_塩基列群 = new byte[訂正バッチサイズ][];
            var l_結果群 = new 訂正結果[訂正バッチサイズ];

            while (l_読み込み.Get_続きがあるか())
            {
                var l_件数 = 0;
                while (l_件数 < 訂正バッチサイズ && l_読み込み.Get_続きがあるか())
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
                    l_書き込み.V_書き込み(
                        l_ID群[i],
                        string.Join(string.Empty, l_結果.A_塩基列.Select(Util.V_変換_塩基文字)),
                        l_クオリティ群[i]);
                }
            }

            return new ファイル訂正統計(l_総リード数, l_訂正されたリード数, l_総訂正塩基数);
        }

        /// <summary>
        /// 1リード(塩基ID空間のバイト列、曖昧塩基は Consts.無効な塩基)を
        /// 貪欲法で訂正する。副作用のない純粋関数(入力は変更しない)。
        /// </summary>
        public static 訂正結果 Get_訂正結果(
            ReadOnlySpan<byte> p_リード, TrustedKmerIndex p_kmerインデックス, int p_k長, int p_最大反復数 = 10)
        {
            if (p_リード.Length < p_k長)
            {
                return new 訂正結果(p_リード.ToArray(), 0);
            }

            var l_塩基列 = p_リード.ToArray();
            var l_窓数 = l_塩基列.Length - p_k長 + 1;
            var l_訂正数 = 0;

            for (var l_反復 = 0; l_反復 < p_最大反復数; l_反復++)
            {
                var l_信頼状況 = Get_窓ごとの信頼状況(l_塩基列, p_k長, p_kmerインデックス);
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

                    var l_信頼できない窓があるか = false;
                    for (var w = l_窓開始; w <= l_窓終了; w++)
                    {
                        if (!l_信頼状況[w])
                        {
                            l_信頼できない窓があるか = true;
                            break;
                        }
                    }
                    if (!l_信頼できない窓があるか)
                    {
                        continue;
                    }

                    var l_現在の塩基 = l_塩基列[l_位置];
                    for (byte l_候補 = Consts.塩基ID.A; l_候補 <= Consts.塩基ID.T; l_候補++)
                    {
                        if (l_候補 == l_現在の塩基)
                        {
                            continue;
                        }

                        var l_改善数 = Get_置換の改善数(
                            l_塩基列, l_位置, l_候補, l_窓開始, l_窓終了, p_k長, l_信頼状況, p_kmerインデックス);
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
                    // (=残った信頼できない窓は、単発の置換では解決できない)。
                    break;
                }

                l_塩基列[l_最良位置] = l_最良塩基;
                l_訂正数++;
            }

            return new 訂正結果(l_塩基列, l_訂正数);
        }

        private static bool[] Get_窓ごとの信頼状況(byte[] p_塩基列, int p_k長, TrustedKmerIndex p_kmerインデックス)
        {
            var l_窓数 = p_塩基列.Length - p_k長 + 1;
            var l_信頼状況 = new bool[l_窓数];
            for (var w = 0; w < l_窓数; w++)
            {
                l_信頼状況[w] = Get_窓が信頼できるか(p_塩基列, w, p_k長, p_kmerインデックス);
            }
            return l_信頼状況;
        }

        private static bool Get_窓が信頼できるか(
            byte[] p_塩基列, int p_窓開始, int p_k長, TrustedKmerIndex p_kmerインデックス)
        {
            for (var i = p_窓開始; i < p_窓開始 + p_k長; i++)
            {
                if (p_塩基列[i] == Consts.無効な塩基)
                {
                    return false;
                }
            }
            return p_kmerインデックス.Get_含まれるか(p_塩基列.AsSpan(p_窓開始, p_k長));
        }

        /// <summary>
        /// p_位置 を p_候補 に置換した場合の「信頼できる窓の純増数」を計算する
        /// (p_窓開始..p_窓終了 の範囲、すなわち その位置を含みうる窓のみが
        /// 影響を受けるため、その範囲だけを再評価すれば十分)。塩基列は評価後、
        /// 呼び出し前の状態に戻す(副作用を残さない)。
        /// </summary>
        private static int Get_置換の改善数(
            byte[] p_塩基列, int p_位置, byte p_候補, int p_窓開始, int p_窓終了, int p_k長,
            bool[] p_置換前の信頼状況, TrustedKmerIndex p_kmerインデックス)
        {
            var l_元の塩基 = p_塩基列[p_位置];
            p_塩基列[p_位置] = p_候補;

            var l_改善数 = 0;
            for (var w = p_窓開始; w <= p_窓終了; w++)
            {
                var l_信頼できるか = Get_窓が信頼できるか(p_塩基列, w, p_k長, p_kmerインデックス);
                if (l_信頼できるか && !p_置換前の信頼状況[w])
                {
                    l_改善数++;
                }
                else if (!l_信頼できるか && p_置換前の信頼状況[w])
                {
                    l_改善数--;
                }
            }

            p_塩基列[p_位置] = l_元の塩基;
            return l_改善数;
        }
    }
}
