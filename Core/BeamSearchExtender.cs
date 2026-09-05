using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// 相互一意性の判定で決めきれなかった分岐を、先読み(ビームサーチ)で解く。
    ///
    /// 相互一意性は1歩だけを見るため、分岐の直後は五分五分でも数歩先で片方だけが
    /// ペアエンドの証拠と整合する状況を取りこぼす。各候補から複数の経路を並行して
    /// 伸ばし、その間のペアエンドの支持を積算して比べる。
    ///
    /// 安全側に倒す設計:
    /// - 上位の経路群が最初の1歩から割れていれば何もしない。ビームサーチの利点は
    ///   有力な仮説が一致する部分にだけコミットすることで、僅差で1本を選ぶことではない。
    /// - どの候補にも支持が無ければ何もしない。
    /// - コピー数を予算とし、反復配列を予算以上に通らない。予算が無いと同じ反復を
    ///   何度でも通れてしまい、ありもしない長い経路ができる。
    ///
    /// ゲノム全体を1本のオイラー路として探すと正解以外の経路も同数だけ存在し
    /// 誤アセンブリを量産するため、あくまで局所的な近似に留める。
    /// </summary>
    internal static class BeamSearchExtender
    {
        /// <summary>
        /// 先読みで進む塩基数の上限。長くするほど遠くの証拠を使えるが、
        /// 探索が広がるうえ、遠いほどペアエンドの証拠は届かなくなる。
        /// インサートサイズの数倍あれば、跨げる範囲は使い切れる。
        /// </summary>
        private const int 先読み倍率 = 4;

        /// <summary>1つの分岐あたりに保持する部分経路の数。</summary>
        private const int ビーム幅 = 8;

        /// <summary>先読みの1経路あたりの最大ステップ数(暴走防止)。</summary>
        private const int 経路あたりの最大ステップ数 = 40;

        /// <summary>
        /// 結合が未確定(-1)の頂点について、先読みで続きを決められるものを決める。
        /// 結合の配列を直接書き換える。戻り値は新たに確定した結合の数
        /// (有向、双子ぶんを含む)。
        /// </summary>
        public static int V_延長_先読み(
            UnitigGraph p_グラフ,
            List<string> p_ユニティグ配列,
            int[] p_結合,
            IReadOnlyDictionary<(int, int), ulong> p_ペア連結,
            IReadOnlyDictionary<int, int> p_コピー数,
            int p_インサートサイズ,
            decimal p_優勢閾値,
            ulong p_最小証拠数)
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
                    // 単一コピーの足場が1つも取れない = いま反復配列の上にいて、
                    // どのコピーにいるのか分からない。この状態で進む方向を選ぶ
                    // 根拠は原理的に存在しない。
                    continue;
                }

                var l_最良 = Get_最良の1歩(
                    p_グラフ, p_ユニティグ配列, v, l_足場, p_ペア連結, p_コピー数,
                    l_先読み塩基数, p_優勢閾値, p_最小証拠数);
                if (l_最良 is not { } l_選択)
                {
                    continue;
                }

                // 相互一意性は保ったままにする。行き先に既に別の結合が
                // 入っている場合は、そちらを壊してまで繋がない。
                if (p_結合[l_選択 ^ 1] != -1 || p_結合[l_選択] == (v ^ 1))
                {
                    continue;
                }

                // 解きほぐされていない多コピーの反復を通り抜ける結合は作らない。
                // A-R-B-R-C という構造で A→R と R→C はどちらも本物の隣接だが、
                // R を1回しか使えない walk でこれを連鎖させると中間の B が
                // 飛ばされる(詳細は ContigMaker 側の同名の判定を参照)。
                if (!Get_通り抜けてよいか(p_グラフ, p_コピー数, v) || !Get_通り抜けてよいか(p_グラフ, p_コピー数, l_選択 ^ 1))
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
        /// contig 末尾のインサートサイズぶんの頂点のうち、単一コピーのものだけを集める。
        /// ここに載ったリードの相方が続きの証拠になる。
        /// 直前の頂点は逆鎖対称性より 結合[v^1] の双子で辿れる。
        ///
        /// 多コピーを足場から外すのが要点。反復内部から読まれたリードはどのコピー
        /// 由来か区別できず、その証拠はどの行き先にも付くため、標本が少ないと
        /// 偶然の偏りで誤った側を選ぶ。
        ///
        /// 通過はするが足場には数えない(多コピー領域の向こう側にある単一コピー
        /// 領域は証拠として有効なため)。
        /// </summary>
        private static List<int> Get_足場(
            int p_頂点, List<string> p_ユニティグ配列, int[] p_結合,
            int p_インサートサイズ, IReadOnlyDictionary<int, int> p_コピー数)
        {
            List<int> l_足場 = [];
            List<int> l_通過済み = [p_頂点];
            if (p_コピー数.GetValueOrDefault(p_頂点 >> 1, 1) <= 1)
            {
                l_足場.Add(p_頂点);
            }

            var l_累積長 = p_ユニティグ配列[p_頂点].Length;
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
                    l_足場.Add(l_直前);
                }
                l_累積長 += p_ユニティグ配列[l_直前].Length;
                l_現在 = l_直前;
            }
            return l_足場;
        }

        /// <summary>
        /// 分岐元からの各候補について先読みし、最初の1歩として最も支持される頂点を返す。
        /// 決めきれない場合は null。
        /// </summary>
        private static int? Get_最良の1歩(
            UnitigGraph p_グラフ,
            List<string> p_ユニティグ配列,
            int p_分岐元,
            List<int> p_足場,
            IReadOnlyDictionary<(int, int), ulong> p_ペア連結,
            IReadOnlyDictionary<int, int> p_コピー数,
            int p_先読み塩基数,
            decimal p_優勢閾値,
            ulong p_最小証拠数)
        {
            List<先読み探索状態> l_ビーム = [];
            foreach (var l_候補 in p_グラフ.A_出辺[p_分岐元])
            {
                var l_予算 = p_コピー数.GetValueOrDefault(l_候補 >> 1, 1);
                if (l_予算 <= 0)
                {
                    continue;
                }
                l_ビーム.Add(new 先読み探索状態
                {
                    A_現在の頂点 = l_候補,
                    A_最初の1歩 = l_候補,
                    A_スコア = Get_スコア(p_足場, l_候補, p_ペア連結),
                    A_進んだ長さ = p_ユニティグ配列[l_候補].Length,
                    A_使用回数 = new Dictionary<int, int> { [l_候補 >> 1] = 1 },
                });
            }
            if (l_ビーム.Count == 0)
            {
                return null;
            }

            // 最初の1歩ごとの最良スコアを追跡する。
            Dictionary<int, long> l_1歩ごとの最良 = [];
            foreach (var l_状態 in l_ビーム)
            {
                l_1歩ごとの最良[l_状態.A_最初の1歩] =
                    Math.Max(l_1歩ごとの最良.GetValueOrDefault(l_状態.A_最初の1歩), l_状態.A_スコア);
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
                            // 予算切れ。反復を何度も通って架空の経路を作らないようにする。
                            continue;
                        }
                        var l_使用回数 = new Dictionary<int, int>(l_状態.A_使用回数);
                        l_使用回数[l_ユニティグID] = l_使用回数.GetValueOrDefault(l_ユニティグID) + 1;
                        l_次のビーム.Add(new 先読み探索状態
                        {
                            A_現在の頂点 = l_候補,
                            A_最初の1歩 = l_状態.A_最初の1歩,
                            A_スコア = l_状態.A_スコア + Get_スコア(p_足場, l_候補, p_ペア連結),
                            A_進んだ長さ = l_状態.A_進んだ長さ + p_ユニティグ配列[l_候補].Length,
                            A_使用回数 = l_使用回数,
                        });
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
                    l_1歩ごとの最良[l_状態.A_最初の1歩] =
                        Math.Max(l_1歩ごとの最良.GetValueOrDefault(l_状態.A_最初の1歩), l_状態.A_スコア);
                }
            }

            var l_順位 = l_1歩ごとの最良.OrderByDescending(x => x.Value).ToList();
            var l_首位 = l_順位[0];
            if ((ulong)Math.Max(0, l_首位.Value) < p_最小証拠数)
            {
                // どの枝にもペアエンドの支持が無い。根拠が無いので繋がない。
                return null;
            }

            var l_合計 = l_順位.Sum(x => Math.Max(0, x.Value));
            if (l_合計 <= 0 || (decimal)l_首位.Value / l_合計 < p_優勢閾値)
            {
                // 上位が割れている。僅差で選ぶくらいなら繋がないほうがよい。
                return null;
            }

            return l_首位.Key;
        }

        /// <summary>
        /// その頂点を「通り抜けて」よいか。多コピーの反復は、解きほぐされて
        /// 入次数・出次数がどちらも1になっている場合(=どのコピーにいるかが
        /// 確定している場合)にだけ通り抜けてよい。
        /// </summary>
        private static bool Get_通り抜けてよいか(
            UnitigGraph p_グラフ, IReadOnlyDictionary<int, int> p_コピー数, int p_頂点)
        {
            if (p_コピー数.GetValueOrDefault(p_頂点 >> 1, 1) <= 1)
            {
                return true;
            }
            return p_グラフ.A_出辺[p_頂点].Count == 1 && p_グラフ.Get_入次数(p_頂点) == 1;
        }

        private static long Get_スコア(
            List<int> p_足場, int p_候補, IReadOnlyDictionary<(int, int), ulong> p_ペア連結)
        {
            long l_スコア = 0;
            foreach (var l_足場頂点 in p_足場)
            {
                l_スコア += (long)p_ペア連結.GetValueOrDefault((l_足場頂点, p_候補));
            }
            return l_スコア;
        }
    }
}
