using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.UnitigBuilding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// 各 unitig のコピー数をカバレッジから推定する
    /// </summary>
    internal static class CopyNumberEstimator
    {
        #region 定数

        /// <summary>
        /// これを下回るカバレッジ比の unitig は、コピー数を推定できるだけの根拠が無いとみなして 1 として扱う (0 コピーにはしない)
        /// </summary>
        private const double 多コピーとみなす比の下限 = 1.5D;

        /// <summary>
        /// コピー数の上限
        /// </summary>
        /// <remarks>
        /// これを超える比が出た場合、rRNA オペロンのような高コピー反復か、あるいはカバレッジ異常のどちらかで区別がつかない<br/>
        /// 経路探索の予算としては大きすぎると探索が発散するため頭打ちにする
        /// </remarks>
        private const int コピー数の上限 = 12;

        /// <summary>
        /// 「染色体側の確定成分と繋がりが無い、独立した島」を単一の複製単位 (プラスミド等) とみなすために要求する、島の合計長の下限
        /// </summary>
        /// <remarks>
        /// AssemblyStatsReporter の「比較可能」しきい値と同じ 500 bp を使う<br/>
        /// 短い島は偶然の孤立 (flanking 配列が trim で消えた等) や短い反復配列との区別がつきにくいため、この下限より短い場合は判定しない (Unicycler が短いセグメントの単一コピー確定に厳しい条件を課しているのと同じ理由)
        /// </remarks>
        private const int 孤立複製単位とみなす最小合計長 = 500;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// unitig ID (1 始まり) -> その unitig を構成する k-mer の平均カバレッジ、を計算する
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_ユニティグ配列"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        public static Dictionary<int, double> Get_カバレッジ(TrustedKmerIndex p_kmerインデックス, IReadOnlyDictionary<int, string> p_ユニティグ配列, int p_k長)
        {
            Dictionary<int, double> l_カバレッジ = [];
            foreach (var (l_ID, l_配列) in p_ユニティグ配列)
            {
                if (l_配列.Length < p_k長)
                {
                    l_カバレッジ[l_ID] = 0D;
                    continue;
                }

                var l_塩基列 = Util.V_変換_塩基列(l_配列);

                var l_合計 = 0UL;
                var l_件数 = 0;
                for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
                {
                    l_合計 += p_kmerインデックス.Get_カバレッジ(l_塩基列.AsSpan(i, p_k長));
                    l_件数++;
                }
                l_カバレッジ[l_ID] = l_件数 == 0 ? 0D : (double)l_合計 / l_件数;
            }
            return l_カバレッジ;
        }

        /// <summary>
        /// カバレッジからコピー数を推定する
        /// </summary>
        /// <param name="p_カバレッジ"></param>
        /// <param name="p_ユニティグ長"></param>
        /// <param name="p_グラフ"></param>
        /// <returns></returns>
        public static コピー数推定結果 Get_推定結果(IReadOnlyDictionary<int, double> p_カバレッジ, IReadOnlyDictionary<int, int> p_ユニティグ長, UnitigGraph? p_グラフ = null)
        {
            // k-mer スペクトルの 2 成分混合モデルが適合できていれば、その単一コピー平均を
            // 基準値に使う
            // カットオフと同じモデルから導くことで、unitig の長さ加重
            // 中央値という別のヒューリスティックとの食い違いを無くす
            // (KmerCutoffSelector.V_解決_kmerカットオフ 参照)
            // 適合に失敗している場合は
            // 従来どおり unitig カバレッジの長さ加重中央値にフォールバックする
            var l_モデル基準値 = ConfigurationManager.A_スペクトルモデル?.A_単一コピー平均;
            var l_基準値 = l_モデル基準値 is { } l_値 && l_値 > 0D
                ? l_値
                : Get_長さ加重中央値(p_カバレッジ, p_ユニティグ長);

            Dictionary<int, int> l_コピー数 = [];
            foreach (var (l_ID, l_カバレッジ値) in p_カバレッジ)
            {
                if (l_基準値 <= 0D)
                {
                    l_コピー数[l_ID] = 1;
                    continue;
                }

                var l_比 = l_カバレッジ値 / l_基準値;
                if (l_比 < 多コピーとみなす比の下限)
                {
                    // 単一コピー (あるいは低カバレッジで判断できない)
                    // 0 にはしない: 実際に配列は存在しており、経路から
                    // 締め出してしまうと組み立てられなくなる
                    l_コピー数[l_ID] = 1;
                    continue;
                }

                l_コピー数[l_ID] = Math.Clamp((int)Math.Round(l_比), 1, コピー数の上限);
            }

            if (p_グラフ is { } l_グラフ)
            {
                V_修正_孤立複製単位コピー数(l_グラフ, p_カバレッジ, p_ユニティグ長, l_コピー数);
                V_修正_接続による単一コピー再判定(l_グラフ, p_カバレッジ, l_コピー数);
            }

            return new コピー数推定結果(l_基準値, p_カバレッジ, l_コピー数);
        }

        /// <summary>
        /// 推定結果の要約をコンソールへ出力する
        /// </summary>
        /// <param name="p_推定結果"></param>
        /// <param name="p_ユニティグ長"></param>
        /// <remarks>
        /// 「単一コピーが何本・何 bp、2 コピー以上が何本・何 bp」が分かると、反復配列がアセンブリのどれだけを占めているかが把握できる
        /// </remarks>
        public static void V_出力_推定結果(コピー数推定結果 p_推定結果, IReadOnlyDictionary<int, int> p_ユニティグ長)
        {
            Logger.V_出力(メッセージID.単一コピー基準値, p_推定結果.A_単一コピー基準値);

            var l_コピー数別 = p_推定結果.A_コピー数
                .GroupBy(x => x.Value)
                .OrderBy(x => x.Key)
                .Select(x => (A_コピー数: x.Key, A_本数: x.Count(), A_塩基数: x.Sum(y => (long)p_ユニティグ長.GetValueOrDefault(y.Key, 0))))
                .ToList();

            var l_要約 = string.Join(", ", l_コピー数別.Select(x => $"x{x.A_コピー数}: {x.A_本数} unitig(s)/{x.A_塩基数:N0}bp"));
            Logger.V_出力(メッセージID.コピー数の要約, l_要約);

            var l_反復塩基数 = l_コピー数別.Where(x => x.A_コピー数 >= 2).Sum(x => x.A_塩基数);
            var l_総塩基数 = l_コピー数別.Sum(x => x.A_塩基数);
            if (l_総塩基数 > 0L)
            {
                Logger.V_出力(メッセージID.反復配列の割合, l_反復塩基数, l_総塩基数, 100.0D * l_反復塩基数 / l_総塩基数);
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 染色体側と繋がりの無い「島」を、独立した複製単位の単一コピー領域とみなす
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_カバレッジ"></param>
        /// <param name="p_ユニティグ長"></param>
        /// <param name="p_コピー数"></param>
        private static void V_修正_孤立複製単位コピー数(UnitigGraph p_グラフ, IReadOnlyDictionary<int, double> p_カバレッジ, IReadOnlyDictionary<int, int> p_ユニティグ長, Dictionary<int, int> p_コピー数)
        {
            var l_成分ID = Get_連結成分(p_グラフ, p_コピー数.Keys);

            HashSet<int> l_確定済み成分 = [.. p_コピー数
                .Where(x => x.Value <= 1)
                .Select(x => l_成分ID[x.Key])];

            var l_未確定の島一覧 = p_コピー数.Keys
                .GroupBy(l_ID => l_成分ID[l_ID])
                .Where(g => !l_確定済み成分.Contains(g.Key));

            foreach (var l_島 in l_未確定の島一覧)
            {
                var l_島の合計長 = l_島.Sum(l_ID => (long)p_ユニティグ長.GetValueOrDefault(l_ID, 0));
                if (l_島の合計長 < 孤立複製単位とみなす最小合計長)
                {
                    continue;
                }

                var l_島内カバレッジ = l_島
                    .Select(l_ID => p_カバレッジ.GetValueOrDefault(l_ID, 0D))
                    .Where(x => x > 0D)
                    .ToList();
                if (l_島内カバレッジ.Count == 0)
                {
                    continue;
                }

                var l_局所基準値 = StatsUtil.Get_中央値(l_島内カバレッジ);
                if (l_局所基準値 <= 0D)
                {
                    continue;
                }

                var l_Is内部一貫 = l_島.All(l_ID =>
                {
                    var l_値 = p_カバレッジ.GetValueOrDefault(l_ID, 0D);
                    return l_値 <= 0D || l_値 / l_局所基準値 < 多コピーとみなす比の下限;
                });
                if (!l_Is内部一貫)
                {
                    // 島の内部でもカバレッジ水準がばらついている
                    // (=島の中に反復がある)
                    // 触らない
                    continue;
                }

                foreach (var l_ID in l_島)
                {
                    p_コピー数[l_ID] = 1;
                }
            }
        }

        /// <summary>
        /// unitig をグラフ上の連結成分に分ける (向きは無視し、辺があれば繋がっているとみなす)
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_ユニティグID一覧"></param>
        /// <remarks>
        /// 辺 v→w があれば逆鎖対称性より w^1→v^1 もあるため、各 unitig の両頂点 (順鎖・逆鎖) の出辺だけを辿れば入ってくる辺も含めて全方向を辿ったことになる
        /// </remarks>
        /// <returns></returns>
        private static Dictionary<int, int> Get_連結成分(UnitigGraph p_グラフ, IEnumerable<int> p_ユニティグID一覧)
        {
            Dictionary<int, int> l_成分ID = [];
            var l_次の成分ID = 0;

            foreach (var l_開始ID in p_ユニティグID一覧)
            {
                if (l_成分ID.ContainsKey(l_開始ID))
                {
                    continue;
                }

                var l_現在の成分ID = l_次の成分ID++;
                Queue<int> l_キュー = new([l_開始ID]);
                l_成分ID[l_開始ID] = l_現在の成分ID;

                while (l_キュー.Count > 0)
                {
                    var l_ID = l_キュー.Dequeue();
                    foreach (var l_頂点 in (ReadOnlySpan<int>)[2 * l_ID, (2 * l_ID) + 1])
                    {
                        foreach (var l_隣接頂点 in p_グラフ.A_出辺[l_頂点])
                        {
                            var l_隣接ID = l_隣接頂点 >> 1;
                            if (l_成分ID.TryAdd(l_隣接ID, l_現在の成分ID))
                            {
                                l_キュー.Enqueue(l_隣接ID);
                            }
                        }
                    }
                }
            }

            return l_成分ID;
        }

        /// <summary>
        /// 大域基準値との比では多コピーに見える unitig を、接続構造を使って再判定する (Unicycler の copy depth propagation の簡略版)
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_カバレッジ"></param>
        /// <param name="p_コピー数"></param>
        private static void V_修正_接続による単一コピー再判定(UnitigGraph p_グラフ, IReadOnlyDictionary<int, double> p_カバレッジ, Dictionary<int, int> p_コピー数)
        {
            foreach (var l_ID in p_コピー数.Keys.ToList())
            {
                if (p_コピー数[l_ID] <= 1)
                {
                    continue;
                }

                if (!p_カバレッジ.TryGetValue(l_ID, out var l_自身のカバレッジ) || l_自身のカバレッジ <= 0D)
                {
                    continue;
                }

                var l_成分 = Get_排他的成分(p_グラフ, l_ID);
                if (l_成分.Count < 2)
                {
                    // 排他的に繋がる相手がいない
                    // 比較材料が無いので判定しない
                    continue;
                }

                var l_成分内カバレッジ = l_成分
                    .Select(x => p_カバレッジ.GetValueOrDefault(x, 0D))
                    .Where(x => x > 0D)
                    .ToList();
                if (l_成分内カバレッジ.Count < 2)
                {
                    continue;
                }

                var l_局所基準値 = StatsUtil.Get_中央値(l_成分内カバレッジ);
                if (l_局所基準値 <= 0D)
                {
                    continue;
                }

                if (l_自身のカバレッジ / l_局所基準値 < 多コピーとみなす比の下限)
                {
                    p_コピー数[l_ID] = 1;
                }
            }
        }

        /// <summary>
        /// 指定した unitig から、排他的な辺 (出次数 1 かつ行き先の入次数も 1) だけを両方向へ辿って到達できる unitig ID の集合 (自分自身を含む) を返す
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_ユニティグID"></param>
        /// <remarks>
        /// 両方の頂点 (順鎖・逆鎖) から辿ることで、鎖を両方向に伸ばす
        /// </remarks>
        /// <returns></returns>
        private static HashSet<int> Get_排他的成分(UnitigGraph p_グラフ, int p_ユニティグID)
        {
            var l_開始1 = 2 * p_ユニティグID;
            var l_開始2 = l_開始1 ^ 1;

            HashSet<int> l_訪問済み頂点 = [l_開始1, l_開始2];
            HashSet<int> l_結果 = [p_ユニティグID];
            Queue<int> l_キュー = new([l_開始1, l_開始2]);

            while (l_キュー.Count > 0)
            {
                var l_現在 = l_キュー.Dequeue();
                var l_出辺 = p_グラフ.A_出辺[l_現在];
                if (l_出辺.Count != 1)
                {
                    continue;
                }

                var l_次 = l_出辺[0];
                if (p_グラフ.Get_入次数(l_次) != 1 || !l_訪問済み頂点.Add(l_次))
                {
                    continue;
                }

                _ = l_結果.Add(l_次 >> 1);
                if (l_訪問済み頂点.Add(l_次 ^ 1))
                {
                    l_キュー.Enqueue(l_次 ^ 1);
                }
                l_キュー.Enqueue(l_次);
            }

            return l_結果;
        }

        /// <summary>
        /// 長さで重み付けしたカバレッジの中央値
        /// </summary>
        /// <param name="p_カバレッジ"></param>
        /// <param name="p_ユニティグ長"></param>
        /// <remarks>
        /// ゲノムの大部分を占める単一コピー領域の水準を推定するために使う
        /// </remarks>
        /// <returns></returns>
        private static double Get_長さ加重中央値(IReadOnlyDictionary<int, double> p_カバレッジ, IReadOnlyDictionary<int, int> p_ユニティグ長)
        {
            var l_組 = p_カバレッジ
                .Where(x => p_ユニティグ長.ContainsKey(x.Key) && x.Value > 0D)
                .Select(x => ((long)p_ユニティグ長[x.Key], x.Value));
            return StatsUtil.Get_長さ加重中央値(l_組);
        }

        #endregion
    }
}
