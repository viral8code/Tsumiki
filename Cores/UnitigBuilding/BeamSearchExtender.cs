using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Evidence;
using Tsumiki.Models.Reporting;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// 相互一意性の判定で決めきれなかった分岐を、先読み (ビームサーチ) で解く
    /// </summary>
    internal static class BeamSearchExtender
    {
        #region 定数

        /// <summary>
        /// 先読みで進む塩基数の上限
        /// </summary>
        private const int C_先読み倍率 = 4;

        /// <summary>
        /// 1 つの分岐あたりに保持する部分経路の数
        /// </summary>
        private const int C_ビーム幅 = 8;

        /// <summary>
        /// 先読みの 1 経路あたりの最大ステップ数 (暴走防止)
        /// </summary>
        private const int C_経路あたりの最大ステップ数 = 40;

        /// <summary>
        /// 並びの起点にする単一コピーの頂点を探して contig の末尾を遡る最大頂点数
        /// </summary>
        private const int C_経路の起点を遡る最大頂点数 = 8;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 結合が未確定 (-1) の頂点について、先読みで続きを決められるものを決める
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_結合"></param>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_コピー数"></param>
        /// <param name="p_インサートサイズ"></param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <param name="p_較正器">
        /// 支持を生カウントではなく期待本数との比で測るための較正器<br/>
        /// 渡さない (あるいは使えない) 場合は生カウントのままスコアリングする
        /// </param>
        /// <param name="p_コピー数区間">
        /// 観測された分散を踏まえたコピー数の妥当な範囲 (P1c)<br/>
        /// 反復を何度まで通ってよいかの予算には、点推定ではなくこの上限を使う (誤推定で真の経路を消さないため)<br/>
        /// 渡さない場合は従来どおり点推定 (p_コピー数) を予算にする
        /// </param>
        /// <param name="p_経路索引">
        /// 渡すと、contig の末尾から分岐元を通り抜けたリードの並びでも最初の 1 歩を決める<br/>
        /// ペアの先読みと別の 1 歩を示した分岐は繋がない
        /// </param>
        /// <returns>新たに確定した結合の数</returns>
        public static int V_延長_先読み(UnitigGraph p_グラフ, List<string> p_unitig配列, int[] p_結合, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, IReadOnlyDictionary<int, int> p_コピー数, int p_インサートサイズ, decimal p_優勢閾値, ulong p_最小証拠数, 証拠較正器? p_較正器 = null, IReadOnlyDictionary<int, コピー数区間>? p_コピー数区間 = null, ReadPathIndex? p_経路索引 = null)
        {
            var l_先読み塩基数 = Math.Max(p_インサートサイズ, 1) * C_先読み倍率;
            var l_確定数 = 0;

            for (var v = 2; v < p_グラフ.A_出辺.Count; v++)
            {
                if (p_結合[v] != -1)
                {
                    continue;
                }

                if (p_グラフ.A_出辺[v].Count == 0)
                {
                    continue;
                }

                var l_足場 = Get_足場(v, p_unitig配列, p_結合, p_インサートサイズ, p_コピー数);
                if (l_足場.Count == 0)
                {
                    continue;
                }

                var l_経路先 = Get_経路で優勢な1歩(p_経路索引, p_グラフ, p_unitig配列, v, p_結合, p_コピー数, p_インサートサイズ, p_優勢閾値, p_最小証拠数);
                var l_最良 = Get_最良1歩(p_グラフ, p_unitig配列, v, l_足場, p_ペア連結, p_コピー数, l_先読み塩基数, p_優勢閾値, p_最小証拠数, p_較正器, p_コピー数区間);
                if (l_経路先 is not null && l_最良 is not null && l_経路先 != l_最良)
                {
                    continue;
                }

                if ((l_経路先 ?? l_最良) is not { } l_選択)
                {
                    continue;
                }

                if (p_結合[l_選択 ^ 1] != -1 || p_結合[l_選択] == (v ^ 1))
                {
                    continue;
                }

                if (!p_グラフ.Is構造上一意な辺(v, l_選択) && (!p_グラフ.Is通過可能(p_コピー数, v, p_unitig配列) || !p_グラフ.Is通過可能(p_コピー数, l_選択 ^ 1, p_unitig配列)))
                {
                    continue;
                }

                p_結合[v] = l_選択;
                p_結合[l_選択 ^ 1] = v ^ 1;
                l_確定数 += 2;
            }

            return l_確定数;
        }

        /// <summary>
        /// 足場群から候補頂点への支持を集計する
        /// </summary>
        /// <param name="p_足場"></param>
        /// <param name="p_候補"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_較正器"></param>
        /// <param name="p_前進距離">分岐点から候補直前までの固有配列長</param>
        /// <returns></returns>
        internal static (long A_生, double A_正規化) Get_スコア(List<(int A_頂点, int A_距離)> p_足場, int p_候補, List<string> p_unitig配列, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, 証拠較正器? p_較正器, int p_前進距離)
        {
            var l_生スコア = 0L;
            var l_正規化スコア = 0D;
            var l_候補長 = p_unitig配列[p_候補].Length;
            foreach (var (l_足場頂点, l_足場距離) in p_足場)
            {
                var l_件数 = p_ペア連結.GetValueOrDefault((l_足場頂点, p_候補));
                l_生スコア += (long)l_件数;
                l_正規化スコア += p_較正器 is { A_Is使用可能: true } l_較正器
                    ? l_較正器.Get_正規化済み支持(l_件数, p_unitig配列[l_足場頂点].Length, l_候補長, p_ギャップ長: l_足場距離 + p_前進距離 - Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1))
                    : l_件数;
            }

            return (l_生スコア, l_正規化スコア);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 結合済みの contig 末尾を遡った最も近い単一コピーの頂点から分岐元を通り抜けた並びで、最初の 1 歩を決める
        /// </summary>
        /// <param name="p_経路索引"></param>
        /// <param name="p_グラフ"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_分岐元"></param>
        /// <param name="p_結合"></param>
        /// <param name="p_コピー数"></param>
        /// <param name="p_起点の最短長">起点にする単一コピーの頂点に要る配列長</param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <returns>決まらなければ null</returns>
        private static int? Get_経路で優勢な1歩(ReadPathIndex? p_経路索引, UnitigGraph p_グラフ, List<string> p_unitig配列, int p_分岐元, int[] p_結合, IReadOnlyDictionary<int, int> p_コピー数, int p_起点の最短長, decimal p_優勢閾値, ulong p_最小証拠数)
        {
            if (p_経路索引 is null)
            {
                return null;
            }

            List<int> l_上流 = [p_分岐元];
            var l_現在 = p_分岐元;
            while (l_上流.Count <= C_経路の起点を遡る最大頂点数)
            {
                var l_双子の結合 = p_結合[l_現在 ^ 1];
                if (l_双子の結合 == -1)
                {
                    return null;
                }

                var l_直前 = l_双子の結合 ^ 1;
                if (l_上流.Exists(x => (x >> 1) == (l_直前 >> 1)))
                {
                    return null;
                }

                l_上流.Insert(0, l_直前);
                if (p_コピー数.GetValueOrDefault(l_直前 >> 1, 1) <= 1 && p_unitig配列[l_直前].Length >= p_起点の最短長)
                {
                    return p_経路索引.Get_優勢な行き先(l_上流, p_グラフ.A_出辺[p_分岐元], p_優勢閾値, p_最小証拠数);
                }

                l_現在 = l_直前;
            }

            return null;
        }

        /// <summary>
        /// この unitig を先読み探索で何回まで通ってよいかの予算を求める
        /// </summary>
        /// <param name="p_unitigID"></param>
        /// <param name="p_コピー数">点推定</param>
        /// <param name="p_コピー数区間">観測された分散を踏まえた区間、無ければ点推定をそのまま予算にする</param>
        /// <returns></returns>
        private static int Get_通行予算(int p_unitigID, IReadOnlyDictionary<int, int> p_コピー数, IReadOnlyDictionary<int, コピー数区間>? p_コピー数区間)
        {
            return p_コピー数区間 is { } l_区間辞書 && l_区間辞書.TryGetValue(p_unitigID, out var l_区間) ? l_区間.A_上限 : p_コピー数.GetValueOrDefault(p_unitigID, 1);
        }

        /// <summary>
        /// contig 末尾のインサートサイズぶんの頂点のうち、単一コピーのものだけを集める
        /// </summary>
        /// <param name="p_頂点"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_結合"></param>
        /// <param name="p_インサートサイズ"></param>
        /// <param name="p_コピー数"></param>
        /// <returns>単一コピーとみなせる足場頂点の一覧</returns>
        private static List<(int A_頂点, int A_距離)> Get_足場(int p_頂点, List<string> p_unitig配列, int[] p_結合, int p_インサートサイズ, IReadOnlyDictionary<int, int> p_コピー数)
        {
            List<(int A_頂点, int A_距離)> l_足場 = [];
            List<int> l_通過済み = [p_頂点];
            if (p_コピー数.GetValueOrDefault(p_頂点 >> 1, 1) <= 1)
            {
                l_足場.Add((p_頂点, 0));
            }

            var l_重なり長 = Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1);
            var l_累積長 = Math.Max(0, p_unitig配列[p_頂点].Length - l_重なり長);
            var l_現在 = p_頂点;
            while (l_累積長 < p_インサートサイズ)
            {
                var l_双子の結合 = p_結合[l_現在 ^ 1];
                if (l_双子の結合 == -1)
                {
                    break;
                }

                var l_直前 = l_双子の結合 ^ 1;
                if (l_直前 == l_現在 || l_通過済み.Contains(l_直前))
                {
                    break;
                }

                l_通過済み.Add(l_直前);

                if (p_コピー数.GetValueOrDefault(l_直前 >> 1, 1) <= 1)
                {
                    l_足場.Add((l_直前, l_累積長));
                }

                l_累積長 += Math.Max(0, p_unitig配列[l_直前].Length - l_重なり長);
                l_現在 = l_直前;
            }

            return l_足場;
        }

        /// <summary>
        /// 分岐元からの各候補について先読みし、最初の 1 歩として最も支持される頂点を返す
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_分岐元"></param>
        /// <param name="p_足場"></param>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_コピー数"></param>
        /// <param name="p_先読み塩基数"></param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <param name="p_較正器"></param>
        /// <param name="p_コピー数区間">観測された分散を踏まえたコピー数の妥当な範囲 (P1c)、渡さない場合は点推定を予算にする</param>
        /// <returns>最初の 1 歩として最も支持される頂点、決めきれない場合は null</returns>
        private static int? Get_最良1歩(UnitigGraph p_グラフ, List<string> p_unitig配列, int p_分岐元, List<(int A_頂点, int A_距離)> p_足場, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, IReadOnlyDictionary<int, int> p_コピー数, int p_先読み塩基数, decimal p_優勢閾値, ulong p_最小証拠数, 証拠較正器? p_較正器, IReadOnlyDictionary<int, コピー数区間>? p_コピー数区間 = null)
        {
            List<先読み探索状態> l_ビーム = [];
            foreach (var l_候補 in p_グラフ.A_出辺[p_分岐元])
            {
                var l_予算 = Get_通行予算(l_候補 >> 1, p_コピー数, p_コピー数区間);
                if (l_予算 <= 0)
                {
                    continue;
                }

                var (l_生, l_正規化) = Get_スコア(p_足場, l_候補, p_unitig配列, p_ペア連結, p_較正器, 0);
                l_ビーム.Add(new 先読み探索状態 { A_現在の頂点 = l_候補, A_最初の1歩 = l_候補, A_スコア = l_正規化, A_生スコア = l_生, A_進んだ長さ = Math.Max(0, p_unitig配列[l_候補].Length - Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1)), A_使用回数 = new Dictionary<int, int> { [l_候補 >> 1] = 1 }, });
            }

            if (l_ビーム.Count == 0)
            {
                return null;
            }

            Dictionary<int, (double A_正規化, long A_生)> l_1歩ごとの最良 = [];
            foreach (var l_状態 in l_ビーム)
            {
                V_更新_各歩最良(l_1歩ごとの最良, l_状態);
            }

            var l_ステップ = 0;
            for (; l_ステップ < C_経路あたりの最大ステップ数 && l_ビーム.Count > 0; l_ステップ++)
            {
                List<先読み探索状態> l_次のビーム = [];
                foreach (var l_状態 in l_ビーム)
                {
                    if (l_状態.A_進んだ長さ >= p_先読み塩基数)
                    {
                        continue;
                    }

                    foreach (var l_候補 in p_グラフ.A_出辺[l_状態.A_現在の頂点])
                    {
                        var l_unitigID = l_候補 >> 1;
                        var l_予算 = Get_通行予算(l_unitigID, p_コピー数, p_コピー数区間);
                        if (l_状態.A_使用回数.GetValueOrDefault(l_unitigID) >= l_予算)
                        {
                            continue;
                        }

                        var l_使用回数 = new Dictionary<int, int>(l_状態.A_使用回数);
                        l_使用回数[l_unitigID] = l_使用回数.GetValueOrDefault(l_unitigID) + 1;
                        var (l_生, l_正規化) = l_状態.A_使用回数.ContainsKey(l_unitigID) ? (0L, 0D) : Get_スコア(p_足場, l_候補, p_unitig配列, p_ペア連結, p_較正器, l_状態.A_進んだ長さ);
                        l_次のビーム.Add(new 先読み探索状態 { A_現在の頂点 = l_候補, A_最初の1歩 = l_状態.A_最初の1歩, A_スコア = l_状態.A_スコア + l_正規化, A_生スコア = l_状態.A_生スコア + l_生, A_進んだ長さ = l_状態.A_進んだ長さ + Math.Max(0, p_unitig配列[l_候補].Length - Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1)), A_使用回数 = l_使用回数, });
                    }
                }

                if (l_次のビーム.Count == 0)
                {
                    break;
                }

                l_次のビーム.Sort((x, y) => y.A_スコア.CompareTo(x.A_スコア));
                if (l_次のビーム.Count > C_ビーム幅)
                {
                    l_次のビーム = l_次のビーム[..C_ビーム幅];
                }

                l_ビーム = l_次のビーム;

                foreach (var l_状態 in l_ビーム)
                {
                    V_更新_各歩最良(l_1歩ごとの最良, l_状態);
                }
            }

            var l_打ち切りにより終了 = l_ステップ >= C_経路あたりの最大ステップ数 && l_ビーム.Count > 0;

            var l_順位 = l_1歩ごとの最良.OrderByDescending(x => x.Value.A_正規化).ToList();
            var l_首位 = l_順位[0];
            var l_次点 = l_順位.Count > 1 ? l_順位[1].Value.A_正規化 : 0D;
            if ((ulong)Math.Max(0L, l_首位.Value.A_生) < p_最小証拠数)
            {
                AmbiguityRecorder.V_記録(l_打ち切りにより終了 ? 曖昧箇所の種別.探索打切り : 曖昧箇所の種別.支持なし, AmbiguityRecorder.Get_場所名(p_分岐元), l_首位.Value.A_正規化, l_次点, l_首位.Value.A_生);
                return null;
            }

            var l_合計 = l_順位.Sum(x => Math.Max(0D, x.Value.A_正規化));
            if (l_合計 <= 0D || (decimal)(l_首位.Value.A_正規化 / l_合計) < p_優勢閾値)
            {
                AmbiguityRecorder.V_記録(l_打ち切りにより終了 ? 曖昧箇所の種別.探索打切り : 曖昧箇所の種別.僅差, AmbiguityRecorder.Get_場所名(p_分岐元), l_首位.Value.A_正規化, l_次点, l_首位.Value.A_生);
                return null;
            }

            return l_首位.Key;
        }

        /// <summary>
        /// 最初の 1 歩ごとの最良スコアを更新する
        /// </summary>
        /// <param name="p_1歩ごとの最良"></param>
        /// <param name="p_状態"></param>
        private static void V_更新_各歩最良(Dictionary<int, (double A_正規化, long A_生)> p_1歩ごとの最良, 先読み探索状態 p_状態)
        {
            if (!p_1歩ごとの最良.TryGetValue(p_状態.A_最初の1歩, out var l_既存) || p_状態.A_スコア > l_既存.A_正規化)
            {
                p_1歩ごとの最良[p_状態.A_最初の1歩] = (p_状態.A_スコア, p_状態.A_生スコア);
            }
        }

        #endregion
    }
}
