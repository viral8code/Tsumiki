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
    /// <remarks>
    /// 相互一意性は 1 歩だけを見るため、分岐の直後は五分五分でも数歩先で片方だけがペアエンドの証拠と整合する状況を取りこぼす<br/>
    /// 各候補から複数の経路を並行して伸ばし、その間のペアエンドの支持を積算して比べる<br/>
    /// 安全側に倒す設計:- 上位の経路群が最初の 1 歩から割れていれば何もしない<br/>
    /// ビームサーチの利点は有力な仮説が一致する部分にだけコミットすることで、僅差で 1 本を選ぶことではない<br/>
    /// - どの候補にも支持が無ければ何もしない<br/>
    /// - コピー数を予算とし、反復配列を予算以上に通らない<br/>
    /// 予算が無いと同じ反復を何度でも通れてしまい、ありもしない長い経路ができる<br/>
    /// ゲノム全体を 1 本のオイラー路として探すと正解以外の経路も同数だけ存在し誤アセンブリを量産するため、あくまで局所的な近似に留める
    /// </remarks>
    internal static class BeamSearchExtender
    {
        #region 定数

        /// <summary>
        /// 先読みで進む塩基数の上限
        /// </summary>
        /// <remarks>
        /// 長くするほど遠くの証拠を使えるが、探索が広がるうえ、遠いほどペアエンドの証拠は届かなくなる<br/>
        /// インサートサイズの数倍あれば、跨げる範囲は使い切れる
        /// </remarks>
        private const int 先読み倍率 = 4;

        /// <summary>
        /// 1 つの分岐あたりに保持する部分経路の数
        /// </summary>
        private const int ビーム幅 = 8;

        /// <summary>
        /// 先読みの 1 経路あたりの最大ステップ数 (暴走防止)
        /// </summary>
        private const int 経路あたりの最大ステップ数 = 40;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 結合が未確定 (-1) の頂点について、先読みで続きを決められるものを決める
        /// </summary>
        /// <remarks>
        /// 結合の配列を直接書き換える<br/>
        /// 戻り値は新たに確定した結合の数 (有向、双子ぶんを含む)
        /// </remarks>
        /// <param name="p_グラフ"></param>
        /// <param name="p_ユニティグ配列"></param>
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
        /// <returns>新たに確定した結合の数</returns>
        public static int V_延長_先読み(UnitigGraph p_グラフ, List<string> p_ユニティグ配列, int[] p_結合, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, IReadOnlyDictionary<int, int> p_コピー数, int p_インサートサイズ, decimal p_優勢閾値, ulong p_最小証拠数, 証拠較正器? p_較正器 = null)
        {
            var l_先読み塩基数 = Math.Max(p_インサートサイズ, 1) * 先読み倍率;
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

                var l_足場 = Get_足場(v, p_ユニティグ配列, p_結合, p_インサートサイズ, p_コピー数);
                if (l_足場.Count == 0)
                {
                    // 単一コピーの足場が 1 つも取れない = いま反復配列の上にいて、
                    // どのコピーにいるのか分からない
                    // この状態で進む方向を選ぶ
                    // 根拠は原理的に存在しない
                    continue;
                }

                var l_最良 = Get_最良の1歩(p_グラフ, p_ユニティグ配列, v, l_足場, p_ペア連結, p_コピー数, l_先読み塩基数, p_優勢閾値, p_最小証拠数, p_較正器);
                if (l_最良 is not { } l_選択)
                {
                    continue;
                }

                // 相互一意性は保ったままにする
                // 行き先に既に別の結合が
                // 入っている場合は、そちらを壊してまで繋がない
                if (p_結合[l_選択 ^ 1] != -1 || p_結合[l_選択] == (v ^ 1))
                {
                    continue;
                }

                // 解きほぐされていない多コピーの反復を通り抜ける結合は作らない
                // A-R-B-R-C という構造で A→R と R→C はどちらも本物の隣接だが、
                // R を 1 回しか使えない walk でこれを連鎖させると中間の B が
                // 飛ばされる (詳細は ContigMaker 側の同名の判定を参照)
                if (!p_グラフ.Get_通り抜けてよいか(p_コピー数, v) || !p_グラフ.Get_通り抜けてよいか(p_コピー数, l_選択 ^ 1))
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
        /// <remarks>
        /// 生カウントの合計 (足切り判定用) と、較正器が使える場合は期待本数との比の合計 (ランキング・優勢判定用、較正器が使えない場合は生カウントと同じ値) を両方返す
        /// </remarks>
        /// <param name="p_足場"></param>
        /// <param name="p_候補"></param>
        /// <param name="p_ユニティグ配列"></param>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_較正器"></param>
        /// <param name="p_前進距離">分岐点から候補直前までの固有配列長</param>
        /// <returns></returns>
        internal static (long A_生, double A_正規化) Get_スコア(List<(int A_頂点, int A_距離)> p_足場, int p_候補, List<string> p_ユニティグ配列, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, 証拠較正器? p_較正器, int p_前進距離)
        {
            var l_生スコア = 0L;
            var l_正規化スコア = 0D;
            var l_候補長 = p_ユニティグ配列[p_候補].Length;
            foreach (var (l_足場頂点, l_足場距離) in p_足場)
            {
                var l_件数 = p_ペア連結.GetValueOrDefault((l_足場頂点, p_候補));
                l_生スコア += (long)l_件数;
                l_正規化スコア += p_較正器 is { A_使えるか: true } l_較正器
                    ? l_較正器.Get_正規化済み支持(l_件数, p_ユニティグ配列[l_足場頂点].Length, l_候補長, p_ギャップ長: l_足場距離 + p_前進距離 - Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1))
                    : l_件数;
            }
            return (l_生スコア, l_正規化スコア);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// contig 末尾のインサートサイズぶんの頂点のうち、単一コピーのものだけを集める
        /// </summary>
        /// <remarks>
        /// ここに載ったリードの相方が続きの証拠になる<br/>
        /// 直前の頂点は逆鎖対称性より 結合[v^1] の双子で辿れる<br/>
        /// 多コピーを足場から外すのが要点<br/>
        /// 反復内部から読まれたリードはどのコピー由来か区別できず、その証拠はどの行き先にも付くため、標本が少ないと偶然の偏りで誤った側を選ぶ<br/>
        /// 通過はするが足場には数えない (多コピー領域の向こう側にある単一コピー領域は証拠として有効なため)
        /// </remarks>
        /// <param name="p_頂点"></param>
        /// <param name="p_ユニティグ配列"></param>
        /// <param name="p_結合"></param>
        /// <param name="p_インサートサイズ"></param>
        /// <param name="p_コピー数"></param>
        /// <returns>単一コピーとみなせる足場頂点の一覧</returns>
        private static List<(int A_頂点, int A_距離)> Get_足場(int p_頂点, List<string> p_ユニティグ配列, int[] p_結合, int p_インサートサイズ, IReadOnlyDictionary<int, int> p_コピー数)
        {
            List<(int A_頂点, int A_距離)> l_足場 = [];
            List<int> l_通過済み = [p_頂点];
            if (p_コピー数.GetValueOrDefault(p_頂点 >> 1, 1) <= 1)
            {
                l_足場.Add((p_頂点, 0));
            }

            var l_重なり長 = Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1);
            var l_累積長 = Math.Max(0, p_ユニティグ配列[p_頂点].Length - l_重なり長);
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
                l_累積長 += Math.Max(0, p_ユニティグ配列[l_直前].Length - l_重なり長);
                l_現在 = l_直前;
            }
            return l_足場;
        }

        /// <summary>
        /// 分岐元からの各候補について先読みし、最初の 1 歩として最も支持される頂点を返す
        /// </summary>
        /// <remarks>
        /// 決めきれない場合は null
        /// </remarks>
        /// <param name="p_グラフ"></param>
        /// <param name="p_ユニティグ配列"></param>
        /// <param name="p_分岐元"></param>
        /// <param name="p_足場"></param>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_コピー数"></param>
        /// <param name="p_先読み塩基数"></param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <param name="p_較正器"></param>
        /// <returns>最初の 1 歩として最も支持される頂点、決めきれない場合は null</returns>
        private static int? Get_最良の1歩(UnitigGraph p_グラフ, List<string> p_ユニティグ配列, int p_分岐元, List<(int A_頂点, int A_距離)> p_足場, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, IReadOnlyDictionary<int, int> p_コピー数, int p_先読み塩基数, decimal p_優勢閾値, ulong p_最小証拠数, 証拠較正器? p_較正器)
        {
            List<先読み探索状態> l_ビーム = [];
            foreach (var l_候補 in p_グラフ.A_出辺[p_分岐元])
            {
                var l_予算 = p_コピー数.GetValueOrDefault(l_候補 >> 1, 1);
                if (l_予算 <= 0)
                {
                    continue;
                }
                var (l_生, l_正規化) = Get_スコア(p_足場, l_候補, p_ユニティグ配列, p_ペア連結, p_較正器, 0);
                l_ビーム.Add(new 先読み探索状態 { A_現在の頂点 = l_候補, A_最初の1歩 = l_候補, A_スコア = l_正規化, A_生スコア = l_生, A_進んだ長さ = Math.Max(0, p_ユニティグ配列[l_候補].Length - Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1)), A_使用回数 = new Dictionary<int, int> { [l_候補 >> 1] = 1 }, });
            }
            if (l_ビーム.Count == 0)
            {
                return null;
            }

            // 最初の 1 歩ごとの最良スコア (正規化値と、それに対応する生カウント) を追跡する
            Dictionary<int, (double A_正規化, long A_生)> l_1歩ごとの最良 = [];
            foreach (var l_状態 in l_ビーム)
            {
                V_更新_1歩ごとの最良(l_1歩ごとの最良, l_状態);
            }

            for (var l_ステップ = 0; l_ステップ < 経路あたりの最大ステップ数 && l_ビーム.Count > 0; l_ステップ++)
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
                        var l_ユニティグID = l_候補 >> 1;
                        var l_予算 = p_コピー数.GetValueOrDefault(l_ユニティグID, 1);
                        if (l_状態.A_使用回数.GetValueOrDefault(l_ユニティグID) >= l_予算)
                        {
                            // 予算切れ
                            // 反復を何度も通って架空の経路を作らないようにする
                            continue;
                        }
                        var l_使用回数 = new Dictionary<int, int>(l_状態.A_使用回数);
                        l_使用回数[l_ユニティグID] = l_使用回数.GetValueOrDefault(l_ユニティグID) + 1;
                        var (l_生, l_正規化) = l_状態.A_使用回数.ContainsKey(l_ユニティグID) ? (0L, 0D) : Get_スコア(p_足場, l_候補, p_ユニティグ配列, p_ペア連結, p_較正器, l_状態.A_進んだ長さ);
                        l_次のビーム.Add(new 先読み探索状態 { A_現在の頂点 = l_候補, A_最初の1歩 = l_状態.A_最初の1歩, A_スコア = l_状態.A_スコア + l_正規化, A_生スコア = l_状態.A_生スコア + l_生, A_進んだ長さ = l_状態.A_進んだ長さ + Math.Max(0, p_ユニティグ配列[l_候補].Length - Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1)), A_使用回数 = l_使用回数, });
                    }
                }

                if (l_次のビーム.Count == 0)
                {
                    break;
                }

                l_次のビーム.Sort((x, y) => y.A_スコア.CompareTo(x.A_スコア));
                if (l_次のビーム.Count > ビーム幅)
                {
                    l_次のビーム = l_次のビーム[..ビーム幅];
                }
                l_ビーム = l_次のビーム;

                foreach (var l_状態 in l_ビーム)
                {
                    V_更新_1歩ごとの最良(l_1歩ごとの最良, l_状態);
                }
            }

            var l_順位 = l_1歩ごとの最良.OrderByDescending(x => x.Value.A_正規化).ToList();
            var l_首位 = l_順位[0];
            var l_次点 = l_順位.Count > 1 ? l_順位[1].Value.A_正規化 : 0D;
            if ((ulong)Math.Max(0L, l_首位.Value.A_生) < p_最小証拠数)
            {
                // どの枝にもペアエンドの支持が無い
                // 根拠が無いので繋がない
                AmbiguityRecorder.V_記録(曖昧箇所の種別.支持なし, AmbiguityRecorder.Get_場所名(p_分岐元), l_首位.Value.A_正規化, l_次点, l_首位.Value.A_生);
                return null;
            }

            var l_合計 = l_順位.Sum(x => Math.Max(0D, x.Value.A_正規化));
            if (l_合計 <= 0D || (decimal)(l_首位.Value.A_正規化 / l_合計) < p_優勢閾値)
            {
                // 上位が割れている
                // 僅差で選ぶくらいなら繋がないほうがよい
                AmbiguityRecorder.V_記録(曖昧箇所の種別.僅差, AmbiguityRecorder.Get_場所名(p_分岐元), l_首位.Value.A_正規化, l_次点, l_首位.Value.A_生);
                return null;
            }

            return l_首位.Key;
        }

        /// <summary>
        /// 最初の 1 歩ごとの最良スコアを更新する
        /// </summary>
        /// <remarks>
        /// 正規化スコアが同点になりうる (較正器が無い場合は生カウントと一致する) ため、比較は正規化スコアで行い、対応する生カウントも一緒に持ち替える
        /// </remarks>
        /// <param name="p_1歩ごとの最良"></param>
        /// <param name="p_状態"></param>
        private static void V_更新_1歩ごとの最良(Dictionary<int, (double A_正規化, long A_生)> p_1歩ごとの最良, 先読み探索状態 p_状態)
        {
            if (!p_1歩ごとの最良.TryGetValue(p_状態.A_最初の1歩, out var l_既存) || p_状態.A_スコア > l_既存.A_正規化)
            {
                p_1歩ごとの最良[p_状態.A_最初の1歩] = (p_状態.A_スコア, p_状態.A_生スコア);
            }
        }

        #endregion
    }
}
