using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer の出現回数の分布から、エラー由来と真のゲノム由来を分ける「谷」と、
    /// 1コピーあたりのカバレッジに相当する「山」を推定する。
    /// </summary>
    internal static class KmerHistogram
    {
        /// <summary>
        /// 推奨カットオフの下限。出現回数1の k-mer はほぼ全てエラー由来で、
        /// 残すとメモリを食ったうえでグラフが偽の枝だらけになる。
        /// </summary>
        public const ulong 推奨カットオフの下限 = 2;

        /// <summary>
        /// 残す k-mer の種類数が推定ゲノムサイズの何倍までなら許容できるか。
        /// ゲノム由来の種類数はゲノムサイズをやや下回る(反復が1種類に潰れる)ため、
        /// この比を超えたぶんはほぼエラー由来の混入とみなせる。
        /// </summary>
        private const double 許容するエラー混入比 = 1.2;

        /// <summary>
        /// ゲノムサイズ推定に含める出現回数の上限(山の位置の倍数)。
        /// これを超えるものはアダプタやコンタミ由来である公算が高く、
        /// 足し込むとゲノムサイズが大きく水増しされる。
        /// </summary>
        private const int ゲノムサイズ推定に含める倍率の上限 = 100;

        /// <summary>
        /// 「山」と認めるために必要な、谷の頻度に対する比。
        /// これを下回る場合は二峰性がはっきりしないとみなして推定を諦める。
        /// </summary>
        private const double 山とみなす頻度比 = 1.5;

        /// <summary>
        /// ヒストグラムを解析して、谷・山・推定ゲノムサイズを求める。
        /// 二峰性がはっきりしない(カバレッジが低すぎる等)場合は null を返す。
        /// </summary>
        public static スペクトル解析結果? Get_解析結果(
            IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限 = 10_000)
        {
            if (p_ヒストグラム.Count == 0)
            {
                return null;
            }

            var l_最大キー = p_ヒストグラム.Keys.Max();
            var l_走査上限 = Math.Min(l_最大キー, p_走査上限);
            if (l_走査上限 < 3)
            {
                return null;
            }

            if (Get_粗い谷(p_ヒストグラム, l_走査上限) is not { } l_粗い谷)
            {
                return null;
            }

            var l_ピーク = Get_ピーク(p_ヒストグラム, l_粗い谷 + 1, l_走査上限);

            // 粗い谷はノイズに引きずられるため、山が分かった時点で取り直す。
            var l_谷 = Get_谷(p_ヒストグラム, l_ピーク);

            var l_谷の頻度 = p_ヒストグラム.GetValueOrDefault(l_谷, 0L);
            var l_ピークの頻度 = p_ヒストグラム.GetValueOrDefault(l_ピーク, 0L);
            if (l_ピーク <= l_谷 || l_ピークの頻度 < l_谷の頻度 * 山とみなす頻度比)
            {
                return null;
            }

            var l_加算上限 = Math.Min(l_最大キー, l_ピーク * ゲノムサイズ推定に含める倍率の上限);
            long l_ゲノム由来の延べ数 = 0;
            long l_延べ数の総和 = 0;
            foreach (var (l_出現回数, l_頻度) in p_ヒストグラム)
            {
                if (l_出現回数 > l_加算上限)
                {
                    continue;
                }
                var l_延べ数 = (long)l_出現回数 * l_頻度;
                l_延べ数の総和 += l_延べ数;
                if (l_出現回数 >= l_谷)
                {
                    l_ゲノム由来の延べ数 += l_延べ数;
                }
            }

            return new スペクトル解析結果(
                A_谷: l_谷,
                A_ピーク出現回数: l_ピーク,
                A_谷の頻度: l_谷の頻度,
                A_ピークの頻度: l_ピークの頻度,
                A_ゲノム由来の延べ数: l_ゲノム由来の延べ数,
                A_延べ数の総和: l_延べ数の総和,
                A_推定ゲノムサイズ: l_ゲノム由来の延べ数 / (long)l_ピーク);
        }

        /// <summary>
        /// 頻度が下げ止まって上がり始めた最初の位置。単調減少のままなら null。
        /// 1段だけの増加はノイズでも起きるため、2つ先まで見て上昇の継続を確かめる。
        /// </summary>
        private static ulong? Get_粗い谷(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限)
        {
            for (var l_出現回数 = 1UL; l_出現回数 + 2 <= p_走査上限; l_出現回数++)
            {
                var l_頻度 = p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L);
                var l_次 = p_ヒストグラム.GetValueOrDefault(l_出現回数 + 1, 0L);
                var l_次の次 = p_ヒストグラム.GetValueOrDefault(l_出現回数 + 2, 0L);
                if (l_次 > l_頻度 && l_次の次 > l_頻度)
                {
                    return l_出現回数;
                }
            }
            return null;
        }

        private static ulong Get_ピーク(
            IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_開始, ulong p_終了)
        {
            var l_ピーク = p_開始;
            var l_最大頻度 = -1L;
            for (var l_出現回数 = p_開始; l_出現回数 <= p_終了; l_出現回数++)
            {
                var l_頻度 = p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L);
                if (l_頻度 > l_最大頻度)
                {
                    l_最大頻度 = l_頻度;
                    l_ピーク = l_出現回数;
                }
            }
            return l_ピーク;
        }

        /// <summary>
        /// 1から山までで頻度が最小になる出現回数。観測された出現回数だけを
        /// 候補にする(疎なヒストグラムでは「データが無いだけ」の穴が
        /// 最小値として選ばれ、谷が山の直前まで押し上げられるため)。
        /// </summary>
        private static ulong Get_谷(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_ピーク)
        {
            var l_谷 = 1UL;
            var l_最小頻度 = long.MaxValue;
            for (var l_出現回数 = 1UL; l_出現回数 <= p_ピーク; l_出現回数++)
            {
                if (!p_ヒストグラム.TryGetValue(l_出現回数, out var l_頻度))
                {
                    continue;
                }
                if (l_頻度 < l_最小頻度)
                {
                    l_最小頻度 = l_頻度;
                    l_谷 = l_出現回数;
                }
            }
            return l_谷;
        }

        /// <summary>
        /// エラー由来の k-mer が集合を支配しない範囲で、できるだけ低いカットオフを返す。
        /// 判定できない場合は null。
        ///
        /// 谷をそのまま使ってはいけない。谷はエラー由来とゲノム由来の曲線が
        /// 交わる点なので、そこで切るとゲノム側の左裾まで削り落とす。欠けた
        /// k-mer の箇所すべてでグラフが切れる一方、偽の枝は tip 除去と
        /// バブル除去が落とせる。両者は対称ではないので低く切るのが原則。
        ///
        /// それでも下限に貼り付けにしないのは、高カバレッジではエラー由来の
        /// k-mer が絶対数として増え、品質を落とさずメモリを減らせるため。
        /// </summary>
        public static ulong? Get_推奨カットオフ(
            IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限 = 10_000)
        {
            if (Get_解析結果(p_ヒストグラム, p_走査上限) is not { } l_解析)
            {
                return null;
            }

            var l_許容種類数 = (long)(l_解析.A_推定ゲノムサイズ * 許容するエラー混入比);

            var l_残る種類数 = p_ヒストグラム.Values.Sum();
            for (var l_出現回数 = 1UL; l_出現回数 < 推奨カットオフの下限; l_出現回数++)
            {
                l_残る種類数 -= p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L);
            }

            // 谷を超えたら残るのはゲノム由来だけなので、それより上げない。
            for (var l_出現回数 = 推奨カットオフの下限; l_出現回数 <= l_解析.A_谷; l_出現回数++)
            {
                if (l_残る種類数 <= l_許容種類数)
                {
                    return l_出現回数;
                }
                l_残る種類数 -= p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L);
            }
            return Math.Max(推奨カットオフの下限, l_解析.A_谷);
        }

        /// <summary>
        /// k-mer スペクトルの解析結果を出力する。推定ゲノムサイズとカバレッジは、
        /// 自動選択された k と -kc の妥当性を利用者が確かめる材料になる。
        /// </summary>
        public static void V_出力_スペクトル(
            IReadOnlyDictionary<ulong, long> p_ヒストグラム, int p_k長, int? p_リード長)
        {
            Console.WriteLine($"[Info] k-mer count histogram (count:#distinct kmers): {Get_要約(p_ヒストグラム)}");

            if (Get_解析結果(p_ヒストグラム) is not { } l_解析)
            {
                Console.WriteLine(
                    "[Info] Could not identify a clear histogram valley " +
                    "(the spectrum may not be bimodal at this coverage).");
                return;
            }

            Console.WriteLine(
                $"[Info] k-mer spectrum: valley at count {l_解析.A_谷} ({l_解析.A_谷の頻度} distinct kmers), " +
                $"single-copy peak at count {l_解析.A_ピーク出現回数} ({l_解析.A_ピークの頻度} distinct kmers)");

            var l_カバレッジ表記 = $"{l_解析.A_ピーク出現回数}x (k-mer)";
            if (p_リード長 is { } l_リード長 && l_リード長 > p_k長)
            {
                // リードのカバレッジを山の位置から逆算してはいけない。エラーを含む
                // k-mer は山ではなく低頻度側へ落ちるため、大きく過小評価になる。
                // 延べ数 = Σ(リードごとの L-k+1) からリード本数を復元すれば、
                // エラーを含む k-mer も勘定に入る。
                var l_リード本数 = l_解析.A_延べ数の総和 / (double)(l_リード長 - p_k長 + 1);
                var l_リードカバレッジ = l_リード本数 * l_リード長 / l_解析.A_推定ゲノムサイズ;
                l_カバレッジ表記 += $" / {l_リードカバレッジ:F1}x (read)";
            }

            Console.WriteLine(
                $"[Info] Estimated genome size: {l_解析.A_推定ゲノムサイズ:N0} bp, estimated coverage: {l_カバレッジ表記}");
        }

        /// <summary>
        /// 出現回数1から上限までのヒストグラムを1行にまとめた要約文字列を作る
        /// (ログ表示用)。
        /// </summary>
        public static string Get_要約(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_表示上限 = 20)
        {
            var l_項目 = new List<string>();
            var l_上限 = Math.Min(p_表示上限, p_ヒストグラム.Count == 0 ? 0 : p_ヒストグラム.Keys.Max());
            for (var l_出現回数 = 1UL; l_出現回数 <= l_上限; l_出現回数++)
            {
                l_項目.Add($"{l_出現回数}:{p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L)}");
            }
            return string.Join(", ", l_項目);
        }
    }
}
