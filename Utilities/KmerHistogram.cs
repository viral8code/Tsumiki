using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// k-mer の出現回数の分布から、エラー由来と真のゲノム由来を分ける「谷」と、1 コピーあたりのカバレッジに相当する「山」を推定する
    /// </summary>
    internal static class KmerHistogram
    {
        #region 定数

        /// <summary>
        /// 推奨カットオフの下限
        /// </summary>
        /// <remarks>
        /// 出現回数 1 の k-mer はほぼ全てエラー由来で、残すとメモリを食ったうえでグラフが偽の枝だらけになる
        /// </remarks>
        public const ulong 推奨カットオフの下限 = 2UL;

        /// <summary>
        /// 残す k-mer の種類数が推定ゲノムサイズの何倍までなら許容できるか
        /// </summary>
        /// <remarks>
        /// ゲノム由来の種類数はゲノムサイズをやや下回る (反復が 1 種類に潰れる) ため、この比を超えたぶんはほぼエラー由来の混入とみなせる
        /// </remarks>
        private const double 許容するエラー混入比 = 1.2D;

        /// <summary>
        /// ゲノムサイズ推定に含める出現回数の上限 (山の位置の倍数)
        /// </summary>
        /// <remarks>
        /// これを超えるものはアダプタやコンタミ由来である公算が高く、足し込むとゲノムサイズが大きく水増しされる
        /// </remarks>
        private const int ゲノムサイズ推定に含める倍率の上限 = 100;

        /// <summary>
        /// 「山」と認めるために必要な、谷の頻度に対する比
        /// </summary>
        /// <remarks>
        /// これを下回る場合は二峰性がはっきりしないとみなして推定を諦める
        /// </remarks>
        private const double 山とみなす頻度比 = 1.5D;

        /// <summary>
        /// 単一コピーとみなすカバレッジの上限を、基準値に対して最低限これだけは取る比
        /// </summary>
        /// <remarks>
        /// ばらつきの小さいデータでは、基準値との比を丸めてコピー数を決めるのと同じになる
        /// </remarks>
        public const double 単一コピー上限の最小比 = 1.5D;

        /// <summary>
        /// 基準値より下側の広がりを測る分位
        /// </summary>
        /// <remarks>
        /// 正規分布の平均 - 1σ に当たる<br/>
        /// 上側には反復配列が混ざるため、混ざらない下側から広がりを測る
        /// </remarks>
        private const double 下側1σの分位 = 0.1587D;

        /// <summary>
        /// 単一コピーとみなす上限を決める片側の z 値 (99%)
        /// </summary>
        private const double 単一コピー上限のz値 = 2.326D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// ヒストグラムを解析して、谷・山・推定ゲノムサイズを求める
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_走査上限">谷・山を探す出現回数の上限</param>
        /// <remarks>
        /// 二峰性がはっきりしない (カバレッジが低すぎる等) 場合は null を返す
        /// </remarks>
        /// <returns></returns>
        public static スペクトル解析結果? Get_解析結果(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限 = 10_000UL)
        {
            if (p_ヒストグラム.Count == 0)
            {
                return null;
            }

            var l_最大キー = p_ヒストグラム.Keys.Max();
            var l_走査上限 = Math.Min(l_最大キー, p_走査上限);
            if (l_走査上限 < 3UL)
            {
                return null;
            }

            if (Get_粗い谷(p_ヒストグラム, l_走査上限) is not { } l_粗い谷)
            {
                return null;
            }

            var l_ピーク = Get_ピーク(p_ヒストグラム, l_粗い谷 + 1UL, l_走査上限);

            // 粗い谷はノイズに引きずられるため、山が分かった時点で取り直す
            var l_谷 = Get_谷(p_ヒストグラム, l_ピーク);

            var l_谷の頻度 = p_ヒストグラム.GetValueOrDefault(l_谷, 0L);
            var l_ピークの頻度 = p_ヒストグラム.GetValueOrDefault(l_ピーク, 0L);
            if (l_ピーク <= l_谷 || l_ピークの頻度 < l_谷の頻度 * 山とみなす頻度比)
            {
                return null;
            }

            var l_加算上限 = Math.Min(l_最大キー, l_ピーク * ゲノムサイズ推定に含める倍率の上限);

            // 延べ数を山の位置で割る推定は、山が単一コピーの平均と一致するときにしか成り立たない
            // GC の偏りでカバレッジが右へ長く裾を引くと、山 (最頻値) は平均より大きく下に来て、ゲノムサイズを倍近くに見積もる
            // k-mer の種類ごとに期待コピー数を数えて足せば、分布の形に依らない
            var (l_単一コピー基準値, l_単一コピー上限) = Get_単一コピーの範囲(p_ヒストグラム, l_谷, l_加算上限, l_ピーク);
            var l_ゲノム由来の延べ数 = 0L;
            var l_延べ数の総和 = 0L;
            var l_推定ゲノムサイズ = 0L;
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
                    l_推定ゲノムサイズ += l_頻度 * Get_期待コピー数(l_出現回数, l_単一コピー基準値, l_単一コピー上限);
                }
            }

            return new スペクトル解析結果(A_谷: l_谷, A_ピーク出現回数: l_ピーク, A_谷の頻度: l_谷の頻度, A_ピークの頻度: l_ピークの頻度, A_ゲノム由来の延べ数: l_ゲノム由来の延べ数, A_延べ数の総和: l_延べ数の総和, A_推定ゲノムサイズ: l_推定ゲノムサイズ, A_単一コピー基準値: l_単一コピー基準値, A_単一コピー上限: l_単一コピー上限);
        }

        /// <summary>
        /// カバレッジから期待されるコピー数
        /// </summary>
        /// <param name="p_カバレッジ"></param>
        /// <param name="p_単一コピー基準値"></param>
        /// <param name="p_単一コピー上限">これ未満を単一コピーとみなす</param>
        /// <returns></returns>
        public static int Get_期待コピー数(double p_カバレッジ, double p_単一コピー基準値, double p_単一コピー上限)
        {
            return p_単一コピー基準値 <= 0D || p_カバレッジ < p_単一コピー上限 ? 1 : Math.Max(2, (int)Math.Round(p_カバレッジ / p_単一コピー基準値));
        }

        /// <summary>
        /// エラー由来の k-mer が集合を支配しない範囲で、できるだけ低いカットオフを返す
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_走査上限">谷・山を探す出現回数の上限</param>
        /// <returns></returns>
        public static ulong? Get_推奨カットオフ(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限 = 10_000UL)
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

            // 谷を超えたら残るのはゲノム由来だけなので、それより上げない
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
        /// k-mer スペクトルの解析結果を出力する
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 推定ゲノムサイズとカバレッジは、自動選択された k と -kc の妥当性を利用者が確かめる材料になる
        /// </remarks>
        public static void V_出力_スペクトル(IReadOnlyDictionary<ulong, long> p_ヒストグラム, int p_k長, int? p_リード長)
        {
            Logger.V_出力(メッセージID.kmerヒストグラム, Get_要約(p_ヒストグラム));

            if (Get_解析結果(p_ヒストグラム) is not { } l_解析)
            {
                Logger.V_出力(メッセージID.スペクトルの谷が不明);
                return;
            }

            Logger.V_出力(メッセージID.スペクトルの谷とピーク, l_解析.A_谷, l_解析.A_谷の頻度, l_解析.A_ピーク出現回数, l_解析.A_ピークの頻度);

            var l_カバレッジ表記 = $"{l_解析.A_ピーク出現回数}x (k-mer)";
            if (p_リード長 is { } l_リード長 && l_リード長 > p_k長)
            {
                // リードのカバレッジを山の位置から逆算してはいけない
                // エラーを含む
                // k-mer は山ではなく低頻度側へ落ちるため、大きく過小評価になる
                // 延べ数 = Σ (リードごとの L-k+1) からリード本数を復元すれば、
                // エラーを含む k-mer も勘定に入る
                var l_リード本数 = l_解析.A_延べ数の総和 / (double)(l_リード長 - p_k長 + 1);
                var l_リードカバレッジ = l_リード本数 * l_リード長 / l_解析.A_推定ゲノムサイズ;
                l_カバレッジ表記 += $" / {l_リードカバレッジ:F1}x (read)";
            }

            Logger.V_出力(メッセージID.ゲノムサイズとカバレッジの推定, l_解析.A_推定ゲノムサイズ, l_カバレッジ表記);
        }

        /// <summary>
        /// 出現回数 1 から上限までのヒストグラムを 1 行にまとめた要約文字列を作る (ログ表示用)
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_表示上限">表示する出現回数の上限</param>
        /// <returns></returns>
        public static string Get_要約(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_表示上限 = 20UL)
        {
            var l_項目 = new List<string>();
            var l_上限 = Math.Min(p_表示上限, p_ヒストグラム.Count == 0 ? 0UL : p_ヒストグラム.Keys.Max());
            for (var l_出現回数 = 1UL; l_出現回数 <= l_上限; l_出現回数++)
            {
                l_項目.Add($"{l_出現回数}:{p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L)}");
            }
            return string.Join(", ", l_項目);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 単一コピーの基準カバレッジと、単一コピーとみなせるカバレッジの上限
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_谷"></param>
        /// <param name="p_加算上限"></param>
        /// <param name="p_ピーク"></param>
        /// <remarks>
        /// 基準は谷より上の k-mer 種類数の中央値を使う<br/>
        /// 細菌ゲノムの反復は種類数では少数に潰れるため、中央値はほぼ単一コピーの水準になる<br/>
        /// 上限は下側の広がりから求めた変動係数で広げ、ばらつきの大きいデータで単一コピーの高カバレッジ側を反復と数えないようにする
        /// </remarks>
        /// <returns></returns>
        private static (double A_基準値, double A_上限) Get_単一コピーの範囲(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_谷, ulong p_加算上限, ulong p_ピーク)
        {
            var l_対象 = p_ヒストグラム
                .Where(x => x.Key >= p_谷 && x.Key <= p_加算上限 && x.Value > 0L)
                .OrderBy(x => x.Key)
                .ToList();
            var l_総数 = l_対象.Sum(x => x.Value);
            if (l_総数 <= 0L)
            {
                return (p_ピーク, 単一コピー上限の最小比 * p_ピーク);
            }

            var l_中央値 = (double)Get_分位の出現回数(l_対象, l_総数, 0.5D);
            var l_下側 = (double)Get_分位の出現回数(l_対象, l_総数, 下側1σの分位);
            var l_変動係数 = (l_中央値 - l_下側) / l_中央値;
            return (l_中央値, l_中央値 * Math.Max(単一コピー上限の最小比, 1D + (単一コピー上限のz値 * l_変動係数)));
        }

        /// <summary>
        /// 種類数で数えた分位に当たる出現回数
        /// </summary>
        /// <param name="p_整列済み">出現回数の昇順に並べたヒストグラム</param>
        /// <param name="p_総数">種類数の総和</param>
        /// <param name="p_分位"></param>
        /// <returns></returns>
        private static ulong Get_分位の出現回数(List<KeyValuePair<ulong, long>> p_整列済み, long p_総数, double p_分位)
        {
            var l_目標 = p_総数 * p_分位;
            var l_累積 = 0L;
            foreach (var (l_出現回数, l_頻度) in p_整列済み)
            {
                l_累積 += l_頻度;
                if (l_累積 >= l_目標)
                {
                    return l_出現回数;
                }
            }
            return p_整列済み[^1].Key;
        }

        /// <summary>
        /// 頻度が下げ止まって上がり始めた最初の位置
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_走査上限">谷を探す出現回数の上限</param>
        /// <remarks>
        /// 単調減少のままなら null<br/>
        /// 1 段だけの増加はノイズでも起きるため、2 つ先まで見て上昇の継続を確かめる
        /// </remarks>
        /// <returns></returns>
        private static ulong? Get_粗い谷(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限)
        {
            for (var l_出現回数 = 1UL; l_出現回数 + 2UL <= p_走査上限; l_出現回数++)
            {
                var l_頻度 = p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L);
                var l_次 = p_ヒストグラム.GetValueOrDefault(l_出現回数 + 1UL, 0L);
                var l_次の次 = p_ヒストグラム.GetValueOrDefault(l_出現回数 + 2UL, 0L);
                if (l_次 > l_頻度 && l_次の次 > l_頻度)
                {
                    return l_出現回数;
                }
            }
            return null;
        }

        /// <summary>
        /// 単一コピーに相当する山の位置を返す
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_開始">探す範囲の下端</param>
        /// <param name="p_終了">探す範囲の上端</param>
        /// <returns>山の位置</returns>
        private static ulong Get_ピーク(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_開始, ulong p_終了)
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
        /// 1 から山までで頻度が最小になる出現回数
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_ピーク">山の位置</param>
        /// <remarks>
        /// 観測された出現回数だけを候補にする (疎なヒストグラムでは「データが無いだけ」の穴が最小値として選ばれ、谷が山の直前まで押し上げられるため)
        /// </remarks>
        /// <returns></returns>
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

        #endregion
    }
}
