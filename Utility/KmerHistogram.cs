using Tsumiki.Model;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer の出現回数ごとの分布(出現回数 -> 何種類のユニークk-merがその
    /// 回数を持つか)から、エラー由来の低頻度k-merと真のゲノム由来k-merを
    /// 分ける「谷」と、1コピーあたりのカバレッジに相当する「山」を推定する。
    ///
    /// 谷はそのまま k-mer カットオフ(-kc)の適正値であり、
    /// 山の位置とスペクトルの総量からゲノムサイズが見積もれる。
    /// </summary>
    internal static class KmerHistogram
    {
        /// <summary>
        /// 推奨カットオフの下限。出現回数1の k-mer はどのカバレッジ帯でも
        /// ほぼ全てシーケンスエラー由来であり、これを残すとメモリを大量に
        /// 消費したうえでグラフが偽の枝だらけになる。解析上の谷が1と出ても
        /// 推奨値としては2まで引き上げる。
        /// </summary>
        public const ulong 推奨カットオフの下限 = 2;

        /// <summary>
        /// カットオフを通過した k-mer の種類数が、推定ゲノムサイズの何倍までなら
        /// 許容できるか。
        ///
        /// 相異なるゲノム由来 k-mer の種類数はゲノムサイズをやや下回る(反復配列が
        /// 1種類に潰れるため)。したがってこの比を超えたぶんは、ほぼそのまま
        /// エラー由来の混入とみなせる。2割の混入まで許すことで、カットオフを
        /// 必要以上に上げずに済ませつつ、集合がエラーに埋め尽くされるのを防ぐ。
        /// </summary>
        private const double 許容するエラー混入比 = 1.2;

        /// <summary>
        /// ゲノムサイズ推定で延べ数に加算する出現回数の上限(ピーク位置の倍数)。
        /// これを超える出現回数はアダプタ配列やコンタミ由来である公算が高く、
        /// 素直に足し込むとゲノムサイズが大きく水増しされる。
        /// 高コピーの反復配列(rRNA オペロンで7コピー程度)は十分内側に入る。
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

            // 粗い谷はノイズの影響を受けるため、山の位置が分かった時点で
            // 「1から山までの最小値」として取り直す。こちらのほうが
            // 走査の出発点に依存しないぶん安定する。
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
        /// 出現回数1から走査し、「頻度が下げ止まって上がり始めた」最初の位置を返す。
        ///
        /// 1段だけの増加はヒストグラムのノイズでも起きるため、2つ先まで見て
        /// 上昇が続いていることを確認する。単調減少のまま終わった場合は null。
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
        /// 1から山までの間で頻度が最小になる出現回数を返す。
        ///
        /// 実際に観測された出現回数だけを候補にする(ヒストグラムに現れない
        /// 出現回数は飛ばす)。実データのスペクトルは谷と山の間が途切れないが、
        /// 疎なヒストグラムでは「データが無いだけ」の穴が最小値として
        /// 選ばれてしまい、谷が山の直前まで押し上げられてしまうため。
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
        /// 谷をそのままカットオフに使ってはいけない。谷はエラー由来の曲線と
        /// ゲノム由来の曲線が交わる点なので、そこで切るとゲノム側の左裾まで
        /// まとめて削り落とす。6.4Mbp の実データで、谷(4)で切ると本物の k-mer が
        /// 数千個欠け、その数千箇所すべてでグラフが切れて N50 が 175,674 から
        /// 70,492 へ半分以下に落ちた。
        ///
        /// 欠損と偽の枝は対称ではない。偽の枝は tip 除去とバブル除去が落として
        /// くれるが、消えた k-mer はどこからも復元できない。したがって
        /// カットオフは「できるだけ低く」が原則になる。
        ///
        /// 実測(35x / 100x, k=63, 総延長はいずれも約 6.45Mbp):
        ///   35x  -kc 2 -> N50 175,674 (89本) / -kc 3 -> 133,838 / -kc 4 -> 70,492 (181本)
        ///   100x -kc 2 -> N50 176,931 (72本) / -kc 6 -> 176,821 (74本)
        /// 低カバレッジでは低いほど良く、高カバレッジでは差が出ない。
        /// それでも下限に貼り付けにしないのは、高カバレッジではエラー由来の
        /// k-mer が絶対数として増え(100x では出現回数2だけで 177万種類)、
        /// 品質を落とさずにメモリを減らせるため。
        /// </summary>
        public static ulong? Get_推奨カットオフ(
            IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限 = 10_000)
        {
            if (Get_解析結果(p_ヒストグラム, p_走査上限) is not { } l_解析)
            {
                return null;
            }

            var l_許容種類数 = (long)(l_解析.A_推定ゲノムサイズ * 許容するエラー混入比);

            // 出現回数 c 以上の k-mer の種類数を、c を上げながら見ていく。
            var l_残る種類数 = p_ヒストグラム.Values.Sum();
            for (var l_出現回数 = 1UL; l_出現回数 < 推奨カットオフの下限; l_出現回数++)
            {
                l_残る種類数 -= p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L);
            }

            // 谷より上げることは決してしない。谷を超えたら、残っているのは
            // ゲノム由来の k-mer だけであり、削っても損しかしない。
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
        /// k-mer スペクトルの解析結果をコンソールへ出力する。
        ///
        /// 推定ゲノムサイズとカバレッジは、自動選択された k と -kc が
        /// このデータに対して妥当だったかをユーザーが後から確かめるための
        /// 材料でもある(例えば推定カバレッジが10xを切っていれば、
        /// k を下げたほうがよいと判断できる)。
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
                // リードのカバレッジは山の位置からは求められない。1リード(長さ L)から
                // 取れる k-mer は L-k+1 本だが、そのうちエラーを含むものは山ではなく
                // 低頻度側へ落ちるためである。実際、エラー率0.5%・k=63 の合成データでは
                // 63塩基すべてが無傷である確率が 0.995^63 = 0.73 しかなく、
                // 山の位置から素直に逆算すると真の 60x に対して 42.6x と出ていた。
                //
                // 代わりに、数えた k-mer の延べ数からリード本数を復元する。
                // 延べ数 = Σ(リードごとの L-k+1) なので、リード本数と総塩基数が求まり、
                // 推定ゲノムサイズで割ればカバレッジになる。エラーを含む k-mer も
                // 延べ数には入っているので、この経路なら取りこぼさない。
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
