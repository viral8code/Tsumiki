using System.Text;
using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// unitig 間の隣接を、リードマッピングからの推測ではなく de Bruijn グラフそのものから構築する
    /// </summary>
    internal sealed class UnitigGraph
    {
        #region 定数

        /// <summary>
        /// 1 つの反復として扱う unitig 列の最大頂点数
        /// </summary>
        private const int 反復の鎖の最大頂点数 = 8;

        #endregion

        #region プロパティ

        /// <summary>
        /// 頂点ごとの出辺 (行き先の頂点インデックス)
        /// </summary>
        public List<List<int>> A_出辺 { get; }

        /// <summary>
        /// 末尾を 1 塩基伸ばすと自分の先頭 k-mer に戻る頂点
        /// </summary>
        public HashSet<int> A_自己ループ { get; }

        /// <summary>
        /// walk に通り抜けさせない頂点
        /// </summary>
        public HashSet<int> A_通り抜け禁止の頂点 { get; } = [];

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 構築済みの出辺と自己ループから頂点を作る
        /// </summary>
        /// <param name="p_出辺">頂点ごとの出辺</param>
        /// <param name="p_自己ループ">自己ループを持つ頂点</param>
        private UnitigGraph(List<List<int>> p_出辺, HashSet<int> p_自己ループ)
        {
            this.A_出辺 = p_出辺;
            this.A_自己ループ = p_自己ループ;
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 頂点の入次数
        /// </summary>
        /// <param name="p_頂点"></param>
        /// <returns></returns>
        public int Get_入次数(int p_頂点)
        {
            return this.A_出辺[p_頂点 ^ 1].Count;
        }

        /// <summary>
        /// その頂点を通り抜けてよいか
        /// </summary>
        /// <param name="p_コピー数"></param>
        /// <param name="p_頂点"></param>
        /// <param name="p_unitig配列">渡すと、分岐の片方が短い脇道だけの頂点も通り抜けを許す</param>
        /// <returns></returns>
        public bool Is通過可能(IReadOnlyDictionary<int, int>? p_コピー数, int p_頂点, IReadOnlyList<string>? p_unitig配列 = null)
        {
            if (this.A_通り抜け禁止の頂点.Contains(p_頂点))
            {
                return false;
            }

            var l_出次数 = this.A_出辺[p_頂点].Count;
            var l_入次数 = this.Get_入次数(p_頂点);
            if (l_出次数 == 1 && l_入次数 == 1)
            {
                return true;
            }

            var l_Isコピー数1以下 = (p_コピー数?.GetValueOrDefault(p_頂点 >> 1, 1) ?? 1) <= 1;
            if (l_出次数 < 2 || l_入次数 < 2)
            {
                return l_Isコピー数1以下;
            }

            return l_Isコピー数1以下 && p_unitig配列 is not null && this.Is短い脇道だけの分岐(p_頂点, p_unitig配列);
        }

        /// <summary>
        /// 始点の唯一の出辺が終点で、終点の唯一の入辺が始点かを返す
        /// </summary>
        /// <param name="p_始点"></param>
        /// <param name="p_終点"></param>
        /// <returns></returns>
        public bool Is構造上一意な辺(int p_始点, int p_終点)
        {
            return this.A_出辺[p_始点].Count == 1
                && this.Get_入次数(p_終点) == 1
                && !this.A_通り抜け禁止の頂点.Contains(p_始点)
                && !this.A_通り抜け禁止の頂点.Contains(p_終点 ^ 1);
        }

        /// <summary>
        /// 隣接グラフを構築する
        /// </summary>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_kmer辞書"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_曖昧kmerの番兵"></param>
        /// <returns></returns>
        public static UnitigGraph Get_グラフ(List<string> p_unitig配列, IReadOnlyDictionary<KmerKey, (int A_unitigID, int A_開始位置)> p_kmer辞書, int p_k長, int p_曖昧kmerの番兵)
        {
            List<List<int>> l_出辺 = [];
            HashSet<int> l_自己ループ = [];
            for (var i = 0; i < p_unitig配列.Count; i++)
            {
                l_出辺.Add([]);
            }

            var l_候補 = new byte[p_k長];

            for (var l_頂点 = 2; l_頂点 < p_unitig配列.Count; l_頂点++)
            {
                var l_配列 = p_unitig配列[l_頂点];
                if (l_配列.Length < p_k長)
                {
                    continue;
                }

                var l_末尾開始 = l_配列.Length - p_k長 + 1;
                var l_Has無効塩基 = false;
                for (var i = 0; i < p_k長 - 1; i++)
                {
                    var l_塩基ID = Util.Get_塩基ID(l_配列[l_末尾開始 + i]);
                    if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                    {
                        l_Has無効塩基 = true;
                        break;
                    }
                    l_候補[i] = l_塩基ID;
                }

                if (l_Has無効塩基)
                {
                    continue;
                }

                for (var l_末尾塩基 = Consts.塩基ID.A; l_末尾塩基 <= Consts.塩基ID.T; l_末尾塩基++)
                {
                    l_候補[p_k長 - 1] = l_末尾塩基;

                    if (!p_kmer辞書.TryGetValue(new KmerKey(l_候補.AsSpan()), out var l_ヒット))
                    {
                        continue;
                    }

                    if (l_ヒット.A_unitigID == p_曖昧kmerの番兵 || l_ヒット.A_開始位置 != 0)
                    {
                        continue;
                    }
                    var l_行き先 = ContigMaker.Get_頂点番号(l_ヒット.A_unitigID);

                    if (l_行き先 == l_頂点)
                    {
                        _ = l_自己ループ.Add(l_頂点);
                        continue;
                    }
                    l_出辺[l_頂点].Add(l_行き先);
                }
            }

            return new UnitigGraph(l_出辺, l_自己ループ);
        }

        /// <summary>
        /// 短い反復を、通り抜けたリードの並びとペアエンドの証拠に基づいて、入口と出口の組ごとに複製して解きほぐす
        /// </summary>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_支持"></param>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_反復長の上限"></param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <param name="p_r_mer検証器">
        /// 渡すと、対応付けが確定したあとに複製を確定する前段として、決着した入口・反復・出口の組それぞれを r-mer で検証する (ABySS RResolver 型の拒否権) <br/>
        /// どれか 1 組でも接合点を跨ぐ r-mer の支持が足りなければ複製を見送る<br/>
        /// 注意: head-repeat ・ repeat-tail はどちらの対応付けでも de Bruijn グラフ上の本物の辺なので、正しい対応付けのリードだけからも個々の接合点の存在は独立に確認できてしまう<br/>
        /// したがってこの検証は「対応付け A」と「対応付け B」のどちらが正しいかを区別する力は本質的に持たない (反復が反復である以上、局所的な文脈だけでは区別できないため) <br/>
        /// 実際に効くのは、ペア支持が示す対応付けについて個々の接合点すら生リードに一切裏付けられない (=そもそもその unitig 同士が隣接している根拠が生データに無い) 場合であり、この限定的だが無視できない安全網として使う
        /// </param>
        /// <param name="p_経路索引">渡すと、入口から反復を通って出口まで 1 本で読んだ並びを証拠に加える (複製に合わせて書き換える)</param>
        /// <param name="p_引き継ぎ経路索引">前段 k の確定経路の並び、この k の証拠が 1 件も無い反復に限って使う (複製に合わせて書き換える)</param>
        /// <returns>解きほぐした反復の数</returns>
        public int V_解決_短い反復(List<string> p_unitig配列, Dictionary<(int, int), ulong> p_支持, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, int p_反復長の上限, decimal p_優勢閾値, ulong p_最小証拠数, RepeatRMerVerifier? p_r_mer検証器 = null, ReadPathIndex? p_経路索引 = null, ReadPathIndex? p_引き継ぎ経路索引 = null)
        {
            var l_解決数 = 0;
            var l_r_mer検証で棄却した数 = 0;
            var l_長さで見送った数 = 0;
            var l_形で見送った数 = 0;
            var l_証拠不足の数 = 0;
            var l_僅差の数 = 0;
            var l_足場で解決した数 = 0;
            var l_経路で解決した数 = 0;
            var l_部分解決の数 = 0;
            var l_矛盾の数 = 0;
            var l_鎖で解決した数 = 0;
            var l_引き継ぎで解決した数 = 0;
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            var l_元の頂点数 = this.A_出辺.Count;

            for (var l_始点 = 2; l_始点 < l_元の頂点数; l_始点++)
            {
                if (this.Get_反復の鎖(l_始点) is not { } l_鎖 || l_鎖[0] > (l_鎖[^1] ^ 1))
                {
                    continue;
                }

                if (Get_鎖の配列長(p_unitig配列, l_鎖, l_k長) > p_反復長の上限)
                {
                    l_長さで見送った数++;
                    continue;
                }

                List<int> l_入口群 = [.. this.A_出辺[l_鎖[0] ^ 1].Select(x => x ^ 1).Distinct()];
                List<int> l_出口群 = [.. this.A_出辺[l_鎖[^1]].Distinct()];

                if (Is退化した形(l_鎖, l_入口群, l_出口群))
                {
                    l_形で見送った数++;
                    continue;
                }

                var l_経路行列 = Get_経路支持行列(p_経路索引, l_入口群, l_鎖, l_出口群);
                var l_ペア行列 = Get_ペア支持行列(p_ペア連結, l_入口群, l_出口群);
                var l_経路対応 = Get_決着した対応(l_経路行列, p_優勢閾値, p_最小証拠数);
                var l_ペア対応 = Get_決着した対応(l_ペア行列, p_優勢閾値, p_最小証拠数);
                var l_足場合計 = 0UL;
                var l_Is足場使用 = false;
                if (l_ペア対応.Contains(-1))
                {
                    var l_足場行列 = this.Get_足場ペア支持行列(l_入口群, l_出口群, p_unitig配列, p_ペア連結, p_反復長の上限, l_鎖);
                    l_足場合計 = Get_合計(l_足場行列);
                    var l_直接の決着数 = l_ペア対応.Count(x => x >= 0);
                    l_ペア対応 = Get_併合した対応(l_ペア対応, Get_決着した対応(l_足場行列, p_優勢閾値, p_最小証拠数));
                    l_Is足場使用 = l_ペア対応 is not null && l_ペア対応.Count(x => x >= 0) > l_直接の決着数;
                }

                var l_対応 = l_ペア対応 is null ? null : Get_併合した対応(l_経路対応, l_ペア対応);
                if (l_対応 is null)
                {
                    l_矛盾の数++;
                    continue;
                }

                var l_実測合計 = Get_合計(l_経路行列) + Get_合計(l_ペア行列) + l_足場合計;
                var l_Is引き継ぎ使用 = false;
                if (l_実測合計 == 0UL && p_引き継ぎ経路索引 is not null)
                {
                    l_対応 = Get_決着した対応(Get_経路支持行列(p_引き継ぎ経路索引, l_入口群, l_鎖, l_出口群), p_優勢閾値, 1UL);
                    l_Is引き継ぎ使用 = l_対応.Any(x => x >= 0);
                }

                List<(int A_入口, int A_出口)> l_決着 = [];
                for (var i = 0; i < l_入口群.Count; i++)
                {
                    if (l_対応[i] >= 0)
                    {
                        l_決着.Add((l_入口群[i], l_出口群[l_対応[i]]));
                    }
                }

                if (l_決着.Count == 0)
                {
                    _ = l_実測合計 < p_最小証拠数 ? l_証拠不足の数++ : l_僅差の数++;
                    continue;
                }

                List<int> l_残る入口 = [.. l_入口群.Where(x => !l_決着.Exists(y => y.A_入口 == x))];
                List<int> l_残る出口 = [.. l_出口群.Where(x => !l_決着.Exists(y => y.A_出口 == x))];

                List<(int A_入口, int A_出口)> l_一本道になる組 = [.. l_決着];
                if (l_残る入口.Count == 1 && l_残る出口.Count == 1)
                {
                    l_一本道になる組.Add((l_残る入口[0], l_残る出口[0]));
                }
                if (l_一本道になる組.Exists(x => !Is組が優勢(l_経路行列, l_入口群.IndexOf(x.A_入口), l_出口群.IndexOf(x.A_出口), p_優勢閾値) || !Is組が優勢(l_ペア行列, l_入口群.IndexOf(x.A_入口), l_出口群.IndexOf(x.A_出口), p_優勢閾値)))
                {
                    l_僅差の数++;
                    continue;
                }

                if (p_r_mer検証器 is not null)
                {
                    List<(int A_入口, int A_出口)> l_検証する組 = [.. l_決着];
                    if (l_残る入口.Count == 1 && l_残る出口.Count == 1)
                    {
                        l_検証する組.Add((l_残る入口[0], l_残る出口[0]));
                    }
                    var l_鎖配列 = Get_経路配列(p_unitig配列, l_鎖, l_k長);
                    if (l_検証する組.Exists(x => !p_r_mer検証器.Has接合点支持(p_unitig配列[x.A_入口], l_鎖配列, p_unitig配列[x.A_出口], Consts.r_mer接合点支持の閾値の既定値)))
                    {
                        l_r_mer検証で棄却した数++;
                        continue;
                    }
                }

                var l_移す数 = l_残る入口.Count == 0 || l_残る出口.Count == 0 ? l_決着.Count - 1 : l_決着.Count;
                HashSet<int> l_現在の入口 = [.. l_入口群];
                HashSet<int> l_現在の出口 = [.. l_出口群];
                for (var i = 0; i < l_移す数; i++)
                {
                    var (l_入口, l_出口) = l_決着[i];
                    var l_複製鎖 = this.V_複製_鎖(p_unitig配列, p_支持, l_鎖, l_入口, l_出口);
                    p_経路索引?.V_付け替え_複製(l_鎖, l_複製鎖, l_入口, l_出口, l_現在の入口, l_現在の出口);
                    p_引き継ぎ経路索引?.V_付け替え_複製(l_鎖, l_複製鎖, l_入口, l_出口, l_現在の入口, l_現在の出口);
                    _ = l_現在の入口.Remove(l_入口);
                    _ = l_現在の出口.Remove(l_出口);
                }

                l_解決数++;
                l_足場で解決した数 += l_Is足場使用 ? 1 : 0;
                l_経路で解決した数 += l_経路対応.Any(x => x >= 0) ? 1 : 0;
                l_部分解決の数 += l_現在の入口.Count > 1 || l_現在の出口.Count > 1 ? 1 : 0;
                l_鎖で解決した数 += l_鎖.Count > 1 ? 1 : 0;
                l_引き継ぎで解決した数 += l_Is引き継ぎ使用 ? 1 : 0;
            }

            if (p_r_mer検証器 is not null && l_r_mer検証で棄却した数 > 0)
            {
                Logger.V_出力(メッセージID.rMer検証による棄却, l_r_mer検証で棄却した数);
            }
            Logger.V_出力(メッセージID.短い反復の見送り内訳, l_長さで見送った数, l_形で見送った数, l_証拠不足の数, l_僅差の数, l_足場で解決した数);
            Logger.V_出力(メッセージID.反復解決の証拠内訳, l_経路で解決した数, l_部分解決の数, l_矛盾の数, l_鎖で解決した数, l_引き継ぎで解決した数);

            return l_解決数;
        }

        /// <summary>
        /// 単純バブルを検出し、リード支持が最も高い経路以外の辺を取り除く
        /// </summary>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_支持"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_敗者への引き継ぎ先"></param>
        /// <param name="p_長さ帯の割合"></param>
        /// <param name="p_長さ帯の下限"></param>
        /// <param name="p_類似度の下限"></param>
        /// <param name="p_経路長の上限"></param>
        /// <returns>取り除いた経路の数</returns>
        public int V_除去_単純バブル(List<string> p_unitig配列, IReadOnlyDictionary<(int, int), ulong> p_支持, int p_k長, List<string>? p_敗者への引き継ぎ先 = null, double p_長さ帯の割合 = 0.1D, int p_長さ帯の下限 = 3, double p_類似度の下限 = 0.7D, int p_経路長の上限 = 2_000)
        {
            var l_除去数 = 0;

            for (var l_分岐元 = 2; l_分岐元 < this.A_出辺.Count; l_分岐元++)
            {
                var l_出辺 = this.A_出辺[l_分岐元];
                if (l_出辺.Count < 2)
                {
                    continue;
                }

                Dictionary<int, List<List<int>>> l_再合流先ごと = [];
                foreach (var l_開始 in l_出辺)
                {
                    if (this.Get_単純経路(l_開始, p_unitig配列, p_k長, p_経路長の上限) is not { } l_結果)
                    {
                        continue;
                    }
                    var (l_経路, l_再合流先) = l_結果;
                    if (l_再合流先 == l_分岐元 || (l_再合流先 >> 1) == (l_経路[0] >> 1))
                    {
                        continue;
                    }
                    if (!l_再合流先ごと.TryGetValue(l_再合流先, out var l_経路一覧))
                    {
                        l_経路一覧 = [];
                        l_再合流先ごと[l_再合流先] = l_経路一覧;
                    }
                    l_経路一覧.Add(l_経路);
                }

                foreach (var (l_再合流先, l_経路群) in l_再合流先ごと)
                {
                    if (l_経路群.Count < 2)
                    {
                        continue;
                    }

                    var l_配列群 = l_経路群.Select(x => Get_経路配列(p_unitig配列, x, p_k長)).ToList();

                    var l_基準長 = l_配列群.Min(x => x.Length);
                    var l_差分 = Math.Max(p_長さ帯の下限, p_長さ帯の割合 * l_基準長);
                    if (l_配列群.Any(x => Math.Abs(x.Length - l_基準長) > l_差分))
                    {
                        continue;
                    }

                    var l_基準配列 = l_配列群[0];
                    if (l_配列群.Skip(1).Any(x => Get_類似度(l_基準配列, x) < p_類似度の下限))
                    {
                        continue;
                    }

                    var l_勝者index = Enumerable.Range(0, l_経路群.Count)
                        .OrderByDescending(i => p_支持.GetValueOrDefault((l_分岐元, l_経路群[i][0])))
                        .ThenByDescending(i => l_配列群[i].Length)
                        .First();

                    for (var i = 0; i < l_経路群.Count; i++)
                    {
                        if (i == l_勝者index)
                        {
                            continue;
                        }
                        var l_敗者経路 = l_経路群[i];
                        this.V_除去_双方向辺(l_分岐元, l_敗者経路[0]);
                        for (var j = 0; j + 1 < l_敗者経路.Count; j++)
                        {
                            this.V_除去_双方向辺(l_敗者経路[j], l_敗者経路[j + 1]);
                        }
                        this.V_除去_双方向辺(l_敗者経路[^1], l_再合流先);
                        l_除去数++;
                        p_敗者への引き継ぎ先?.Add(l_配列群[i]);
                    }
                }
            }

            return l_除去数;
        }

        /// <summary>
        /// 分岐のうち、行き止まりで終わる短い枝への辺を外す
        /// </summary>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_枝長の上限">外してよい枝の unitig 長の上限</param>
        /// <returns>外した辺の数</returns>
        public int V_除去_行き止まり枝(List<string> p_unitig配列, int p_枝長の上限)
        {
            var l_除去数 = 0;
            List<int> l_枝 = [];
            for (var v = 2; v < this.A_出辺.Count; v++)
            {
                var l_出辺 = this.A_出辺[v];
                if (l_出辺.Count < 2)
                {
                    continue;
                }

                l_枝.Clear();
                var l_続き = -1;
                var l_続き数 = 0;
                foreach (var w in l_出辺)
                {
                    if (this.Is行き止まり枝(w, p_unitig配列, p_枝長の上限))
                    {
                        l_枝.Add(w);
                    }
                    else
                    {
                        l_続き = w;
                        l_続き数++;
                    }
                }

                if (l_枝.Count == 0 || l_続き数 != 1 || this.Get_入次数(l_続き) != 1)
                {
                    continue;
                }

                if (this.Get_入次数(v) > 1)
                {
                    _ = this.A_通り抜け禁止の頂点.Add(v);
                    _ = this.A_通り抜け禁止の頂点.Add(v ^ 1);
                }

                foreach (var w in l_枝)
                {
                    this.V_除去_双方向辺(v, w);
                    l_除去数++;
                }
            }
            return l_除去数;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 入口と出口がともに複数ある反復の鎖を、入口側から順に取り出す
        /// </summary>
        /// <param name="p_始点"></param>
        /// <returns>該当しなければ null</returns>
        private List<int>? Get_反復の鎖(int p_始点)
        {
            if (this.Get_入次数(p_始点) < 2)
            {
                return null;
            }
            if (this.A_出辺[p_始点].Count >= 2)
            {
                return [p_始点];
            }

            List<int> l_鎖 = [p_始点];
            var l_現在 = p_始点;
            while (this.A_出辺[l_現在].Count == 1 && l_鎖.Count < 反復の鎖の最大頂点数)
            {
                var l_次 = this.A_出辺[l_現在][0];
                if (this.Get_入次数(l_次) != 1 || l_鎖.Exists(x => (x >> 1) == (l_次 >> 1)))
                {
                    return null;
                }
                l_鎖.Add(l_次);
                if (this.A_出辺[l_次].Count >= 2)
                {
                    return l_鎖;
                }
                l_現在 = l_次;
            }
            return null;
        }

        /// <summary>
        /// 鎖が表す配列の長さ
        /// </summary>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_鎖"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        private static int Get_鎖の配列長(List<string> p_unitig配列, List<int> p_鎖, int p_k長)
        {
            return p_鎖.Sum(x => p_unitig配列[x].Length) - ((p_鎖.Count - 1) * (p_k長 - 1));
        }

        /// <summary>
        /// 付け替えの意味が定まらない形か
        /// </summary>
        /// <param name="p_鎖"></param>
        /// <param name="p_入口群"></param>
        /// <param name="p_出口群"></param>
        /// <returns></returns>
        private static bool Is退化した形(List<int> p_鎖, List<int> p_入口群, List<int> p_出口群)
        {
            var l_鎖のID = p_鎖.Select(x => x >> 1).ToHashSet();
            return p_入口群.Concat(p_出口群).Any(x => l_鎖のID.Contains(x >> 1))
                || p_入口群.Select(x => x >> 1).Distinct().Count() != p_入口群.Count
                || p_出口群.Select(x => x >> 1).Distinct().Count() != p_出口群.Count;
        }

        /// <summary>
        /// 入口から鎖を通って出口まで続けて読んだ並びの件数を、入口×出口の行列にする
        /// </summary>
        /// <param name="p_経路索引"></param>
        /// <param name="p_入口群"></param>
        /// <param name="p_鎖"></param>
        /// <param name="p_出口群"></param>
        /// <returns></returns>
        private static ulong[,] Get_経路支持行列(ReadPathIndex? p_経路索引, List<int> p_入口群, List<int> p_鎖, List<int> p_出口群)
        {
            var l_行列 = new ulong[p_入口群.Count, p_出口群.Count];
            if (p_経路索引 is null)
            {
                return l_行列;
            }

            var l_部分列 = new int[p_鎖.Count + 2];
            p_鎖.CopyTo(l_部分列, 1);
            for (var i = 0; i < p_入口群.Count; i++)
            {
                l_部分列[0] = p_入口群[i];
                for (var j = 0; j < p_出口群.Count; j++)
                {
                    l_部分列[^1] = p_出口群[j];
                    l_行列[i, j] = p_経路索引.Get_出現数(l_部分列);
                }
            }
            return l_行列;
        }

        /// <summary>
        /// 入口と出口を直接結ぶペアの件数を、入口×出口の行列にする
        /// </summary>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_入口群"></param>
        /// <param name="p_出口群"></param>
        /// <returns></returns>
        private static ulong[,] Get_ペア支持行列(IReadOnlyDictionary<(int, int), ulong> p_ペア連結, List<int> p_入口群, List<int> p_出口群)
        {
            var l_行列 = new ulong[p_入口群.Count, p_出口群.Count];
            for (var i = 0; i < p_入口群.Count; i++)
            {
                for (var j = 0; j < p_出口群.Count; j++)
                {
                    l_行列[i, j] = p_ペア連結.GetValueOrDefault((p_入口群[i], p_出口群[j]));
                }
            }
            return l_行列;
        }

        /// <summary>
        /// 行列の要素の合計
        /// </summary>
        /// <param name="p_行列"></param>
        /// <returns></returns>
        private static ulong Get_合計(ulong[,] p_行列)
        {
            var l_合計 = 0UL;
            foreach (var l_値 in p_行列)
            {
                l_合計 += l_値;
            }
            return l_合計;
        }

        /// <summary>
        /// 支持行列から、入口ごとに決着した出口を求める
        /// </summary>
        /// <param name="p_行列">入口×出口の支持</param>
        /// <param name="p_優勢閾値"></param>
        /// <param name="p_最小証拠数"></param>
        /// <returns>入口の添字ごとの出口の添字、決着しなければ -1</returns>
        private static int[] Get_決着した対応(ulong[,] p_行列, decimal p_優勢閾値, ulong p_最小証拠数)
        {
            var l_入口数 = p_行列.GetLength(0);
            var l_出口数 = p_行列.GetLength(1);
            var l_対応 = new int[l_入口数];
            Array.Fill(l_対応, -1);
            if (Get_合計(p_行列) < p_最小証拠数)
            {
                return l_対応;
            }

            var l_組の最小証拠数 = Math.Max(1UL, (p_最小証拠数 + 1) / 2);
            var l_列和 = new ulong[l_出口数];
            for (var i = 0; i < l_入口数; i++)
            {
                for (var j = 0; j < l_出口数; j++)
                {
                    l_列和[j] += p_行列[i, j];
                }
            }

            for (var i = 0; i < l_入口数; i++)
            {
                var l_行和 = 0UL;
                var l_最良 = -1;
                var l_最良値 = 0UL;
                var l_Is同点 = false;
                for (var j = 0; j < l_出口数; j++)
                {
                    var l_値 = p_行列[i, j];
                    l_行和 += l_値;
                    if (l_値 > l_最良値)
                    {
                        (l_最良, l_最良値, l_Is同点) = (j, l_値, false);
                    }
                    else if (l_値 == l_最良値 && l_値 > 0UL)
                    {
                        l_Is同点 = true;
                    }
                }

                if (l_最良 < 0 || l_Is同点 || l_最良値 < l_組の最小証拠数 || (decimal)l_最良値 / l_行和 < p_優勢閾値 || (decimal)l_最良値 / l_列和[l_最良] < p_優勢閾値)
                {
                    continue;
                }
                l_対応[i] = l_最良;
            }

            for (var i = 0; i < l_入口数; i++)
            {
                if (l_対応[i] >= 0 && Array.FindAll(l_対応, x => x == l_対応[i]).Length > 1)
                {
                    var l_出口 = l_対応[i];
                    for (var j = 0; j < l_入口数; j++)
                    {
                        l_対応[j] = l_対応[j] == l_出口 ? -1 : l_対応[j];
                    }
                }
            }
            return l_対応;
        }

        /// <summary>
        /// 組が、その入口の行でもその出口の列でも優勢か
        /// </summary>
        /// <param name="p_行列">入口×出口の支持</param>
        /// <param name="p_行">入口の添字</param>
        /// <param name="p_列">出口の添字</param>
        /// <param name="p_優勢閾値"></param>
        /// <returns></returns>
        private static bool Is組が優勢(ulong[,] p_行列, int p_行, int p_列, decimal p_優勢閾値)
        {
            var l_行和 = 0UL;
            for (var j = 0; j < p_行列.GetLength(1); j++)
            {
                l_行和 += p_行列[p_行, j];
            }
            var l_列和 = 0UL;
            for (var i = 0; i < p_行列.GetLength(0); i++)
            {
                l_列和 += p_行列[i, p_列];
            }
            var l_値 = p_行列[p_行, p_列];
            return (l_行和 == 0UL || (decimal)l_値 / l_行和 >= p_優勢閾値) && (l_列和 == 0UL || (decimal)l_値 / l_列和 >= p_優勢閾値);
        }

        /// <summary>
        /// 2 つの証拠から得た対応を合わせる
        /// </summary>
        /// <param name="p_対応1"></param>
        /// <param name="p_対応2"></param>
        /// <returns>同じ入口に別の出口を示すか、合わせると 2 つの入口が同じ出口を指すなら null</returns>
        private static int[]? Get_併合した対応(int[] p_対応1, int[] p_対応2)
        {
            var l_対応 = (int[])p_対応1.Clone();
            for (var i = 0; i < l_対応.Length; i++)
            {
                if (p_対応2[i] < 0)
                {
                    continue;
                }
                if (l_対応[i] >= 0 && l_対応[i] != p_対応2[i])
                {
                    return null;
                }
                l_対応[i] = p_対応2[i];
            }
            var l_決着した出口 = l_対応.Where(x => x >= 0).ToList();
            return l_決着した出口.Distinct().Count() == l_決着した出口.Count ? l_対応 : null;
        }

        /// <summary>
        /// 鎖を複製し、入口 1 本と出口 1 本を複製側へ付け替える
        /// </summary>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_支持"></param>
        /// <param name="p_鎖"></param>
        /// <param name="p_入口"></param>
        /// <param name="p_出口"></param>
        /// <returns>鎖と同じ順の複製の頂点</returns>
        private List<int> V_複製_鎖(List<string> p_unitig配列, Dictionary<(int, int), ulong> p_支持, List<int> p_鎖, int p_入口, int p_出口)
        {
            List<int> l_複製鎖 = [];
            foreach (var l_頂点 in p_鎖)
            {
                l_複製鎖.Add(p_unitig配列.Count);
                p_unitig配列.Add(p_unitig配列[l_頂点]);
                p_unitig配列.Add(p_unitig配列[l_頂点 ^ 1]);
                this.A_出辺.Add([]);
                this.A_出辺.Add([]);
            }

            var l_入辺の支持 = p_支持.GetValueOrDefault((p_入口, p_鎖[0]));
            var l_出辺の支持 = p_支持.GetValueOrDefault((p_鎖[^1], p_出口));
            this.V_除去_双方向辺(p_入口, p_鎖[0]);
            this.V_除去_双方向辺(p_鎖[^1], p_出口);
            this.V_追加_双方向辺(p_入口, l_複製鎖[0]);
            this.V_追加_双方向辺(l_複製鎖[^1], p_出口);

            V_設定_双方向支持(p_支持, p_入口, l_複製鎖[0], l_入辺の支持);
            V_設定_双方向支持(p_支持, l_複製鎖[^1], p_出口, l_出辺の支持);
            for (var i = 0; i + 1 < p_鎖.Count; i++)
            {
                this.V_追加_双方向辺(l_複製鎖[i], l_複製鎖[i + 1]);
                V_設定_双方向支持(p_支持, l_複製鎖[i], l_複製鎖[i + 1], p_支持.GetValueOrDefault((p_鎖[i], p_鎖[i + 1])));
            }
            return l_複製鎖;
        }

        /// <summary>
        /// 辺とその逆鎖側の双子に同じ支持を設定する
        /// </summary>
        /// <param name="p_支持"></param>
        /// <param name="p_始点"></param>
        /// <param name="p_終点"></param>
        /// <param name="p_値"></param>
        private static void V_設定_双方向支持(Dictionary<(int, int), ulong> p_支持, int p_始点, int p_終点, ulong p_値)
        {
            p_支持[(p_始点, p_終点)] = p_値;
            p_支持[(p_終点 ^ 1, p_始点 ^ 1)] = p_値;
        }

        /// <summary>
        /// 反復の入口・出口から一意な鎖を断片長ぶん広げ、その鎖どうしを結ぶペアの件数を入口×出口の行列にする
        /// </summary>
        /// <param name="p_入口群"></param>
        /// <param name="p_出口群"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_ペア連結"></param>
        /// <param name="p_足場長">鎖を広げる長さ</param>
        /// <param name="p_反復の鎖"></param>
        /// <returns></returns>
        private ulong[,] Get_足場ペア支持行列(List<int> p_入口群, List<int> p_出口群, List<string> p_unitig配列, IReadOnlyDictionary<(int, int), ulong> p_ペア連結, int p_足場長, List<int> p_反復の鎖)
        {
            var l_反復のID群 = p_反復の鎖.Select(x => x >> 1).ToHashSet();
            var l_上流群 = p_入口群.Select(x => this.Get_一意な鎖(x, p_unitig配列, p_足場長, l_反復のID群, p_Is上流: true)).ToList();
            var l_下流群 = p_出口群.Select(x => this.Get_一意な鎖(x, p_unitig配列, p_足場長, l_反復のID群, p_Is上流: false)).ToList();
            var l_行列 = new ulong[p_入口群.Count, p_出口群.Count];
            for (var i = 0; i < p_入口群.Count; i++)
            {
                for (var j = 0; j < p_出口群.Count; j++)
                {
                    l_行列[i, j] = Get_鎖間ペア数(l_上流群[i], l_下流群[j], p_ペア連結);
                }
            }
            return l_行列;
        }

        /// <summary>
        /// 頂点から、他へ分かれも合流もしない辺だけを辿れる範囲の頂点を集める
        /// </summary>
        /// <param name="p_頂点">起点 (結果に含む)</param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_足場長">集める配列長の上限</param>
        /// <param name="p_反復のID群">辿らない反復の unitig</param>
        /// <param name="p_Is上流">上流へ辿るか</param>
        /// <returns></returns>
        private List<int> Get_一意な鎖(int p_頂点, List<string> p_unitig配列, int p_足場長, HashSet<int> p_反復のID群, bool p_Is上流)
        {
            List<int> l_鎖 = [p_頂点];
            var l_重なり長 = Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1);
            var l_累積長 = Math.Max(0, p_unitig配列[p_頂点].Length - l_重なり長);
            var l_現在 = p_頂点;
            while (l_累積長 < p_足場長)
            {
                var l_辿る辺 = p_Is上流 ? this.A_出辺[l_現在 ^ 1] : this.A_出辺[l_現在];
                if (l_辿る辺.Count != 1)
                {
                    break;
                }

                var l_次 = p_Is上流 ? l_辿る辺[0] ^ 1 : l_辿る辺[0];
                var l_次の反対側の次数 = p_Is上流 ? this.A_出辺[l_次].Count : this.Get_入次数(l_次);
                if (l_次の反対側の次数 != 1 || p_反復のID群.Contains(l_次 >> 1) || l_鎖.Contains(l_次))
                {
                    break;
                }
                l_鎖.Add(l_次);
                l_累積長 += Math.Max(0, p_unitig配列[l_次].Length - l_重なり長);
                l_現在 = l_次;
            }
            return l_鎖;
        }

        /// <summary>
        /// 上流の鎖から下流の鎖へ向かうペアの総数
        /// </summary>
        /// <param name="p_上流"></param>
        /// <param name="p_下流"></param>
        /// <param name="p_ペア連結"></param>
        /// <returns></returns>
        private static ulong Get_鎖間ペア数(List<int> p_上流, List<int> p_下流, IReadOnlyDictionary<(int, int), ulong> p_ペア連結)
        {
            var l_合計 = 0UL;
            foreach (var u in p_上流)
            {
                foreach (var w in p_下流)
                {
                    l_合計 += p_ペア連結.GetValueOrDefault((u, w));
                }
            }
            return l_合計;
        }

        /// <summary>
        /// 入口 2・出口 2 の分岐が、この頂点から出てこの頂点へ戻る短い脇道 1 本によるものか
        /// </summary>
        /// <param name="p_頂点"></param>
        /// <param name="p_unitig配列"></param>
        /// <returns></returns>
        private bool Is短い脇道だけの分岐(int p_頂点, IReadOnlyList<string> p_unitig配列)
        {
            var l_出辺 = this.A_出辺[p_頂点];
            if (l_出辺.Count != 2 || this.Get_入次数(p_頂点) != 2)
            {
                return false;
            }

            var l_脇道の長さ上限 = ConfigurationManager.A_実行時引数.A_k長 * 2;
            var l_脇道数 = 0;
            foreach (var w in l_出辺)
            {
                if ((w >> 1) == (p_頂点 >> 1) || !this.A_出辺[w].Contains(p_頂点))
                {
                    continue;
                }
                if (p_unitig配列[w].Length > l_脇道の長さ上限)
                {
                    return false;
                }
                l_脇道数++;
            }
            return l_脇道数 == 1;
        }

        /// <summary>
        /// 頂点が、この先に続きを持たず入口も 1 つしかない短い枝か
        /// </summary>
        /// <param name="p_頂点"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_枝長の上限"></param>
        /// <returns></returns>
        private bool Is行き止まり枝(int p_頂点, List<string> p_unitig配列, int p_枝長の上限)
        {
            return this.A_出辺[p_頂点].Count == 0 && this.Get_入次数(p_頂点) == 1 && p_unitig配列[p_頂点].Length <= p_枝長の上限;
        }

        /// <summary>
        /// 辺 v→w を、その逆鎖側の双子 w^1→v^1 と対にして取り除く
        /// </summary>
        /// <param name="p_始点"></param>
        /// <param name="p_終点"></param>
        private void V_除去_双方向辺(int p_始点, int p_終点)
        {
            _ = this.A_出辺[p_始点].Remove(p_終点);
            _ = this.A_出辺[p_終点 ^ 1].Remove(p_始点 ^ 1);
        }

        /// <summary>
        /// 辺 v→w を、その逆鎖側の双子 w^1→v^1 と対にして追加する
        /// </summary>
        /// <param name="p_始点"></param>
        /// <param name="p_終点"></param>
        private void V_追加_双方向辺(int p_始点, int p_終点)
        {
            this.A_出辺[p_始点].Add(p_終点);
            this.A_出辺[p_終点 ^ 1].Add(p_始点 ^ 1);
        }

        /// <summary>
        /// p_開始 から、途中に本物の分岐が無い限り辿れるだけ辿った経路と、その先の再合流先 (=最初に他からも入ってくる頂点) を返す
        /// </summary>
        /// <param name="p_開始"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_長さ上限"></param>
        /// <returns></returns>
        private (List<int> A_経路, int A_再合流先)? Get_単純経路(int p_開始, List<string> p_unitig配列, int p_k長, int p_長さ上限)
        {
            List<int> l_経路 = [];
            HashSet<int> l_訪問済み = [];
            var l_現在 = p_開始;
            var l_累積長 = 0;
            var l_重なり長 = p_k長 - 1;

            while (true)
            {
                if (this.Get_入次数(l_現在) != 1 || !l_訪問済み.Add(l_現在))
                {
                    return l_経路.Count > 0 ? (l_経路, l_現在) : null;
                }

                var l_この頂点の長さ = p_unitig配列[l_現在].Length;
                l_累積長 += l_経路.Count == 0 ? l_この頂点の長さ : Math.Max(0, l_この頂点の長さ - l_重なり長);

                if (l_累積長 > p_長さ上限)
                {
                    return null;
                }

                l_経路.Add(l_現在);

                var l_出辺 = this.A_出辺[l_現在];
                if (l_出辺.Count != 1)
                {
                    return null;
                }
                l_現在 = l_出辺[0];
            }
        }

        /// <summary>
        /// 経路が表す 1 本の配列を組み立てて返す
        /// </summary>
        /// <param name="p_unitig配列">unitig ID 順の配列</param>
        /// <param name="p_経路">辿る頂点の並び</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>組み立てた配列</returns>
        private static string Get_経路配列(List<string> p_unitig配列, List<int> p_経路, int p_k長)
        {
            var l_重なり長 = p_k長 - 1;
            var l_出力 = new StringBuilder(p_unitig配列[p_経路[0]]);
            for (var i = 1; i < p_経路.Count; i++)
            {
                var l_配列 = p_unitig配列[p_経路[i]];
                _ = l_出力.Append(l_配列.Length > l_重なり長 ? l_配列[l_重なり長..] : string.Empty);
            }
            return l_出力.ToString();
        }

        /// <summary>
        /// 2 本の配列の類似度を返す
        /// </summary>
        /// <param name="p_先頭配列">比べる配列</param>
        /// <param name="p_中間配列">比べる配列</param>
        /// <returns>0 から 1 の類似度</returns>
        private static double Get_類似度(string p_先頭配列, string p_中間配列)
        {
            var l_最大長 = Math.Max(p_先頭配列.Length, p_中間配列.Length);
            return l_最大長 == 0 ? 1D : 1D - ((double)Get_編集距離(p_先頭配列, p_中間配列) / l_最大長);
        }

        /// <summary>
        /// 2 本の配列の編集距離を返す
        /// </summary>
        /// <param name="p_先頭配列">比べる配列</param>
        /// <param name="p_中間配列">比べる配列</param>
        /// <returns>編集距離</returns>
        private static int Get_編集距離(string p_先頭配列, string p_中間配列)
        {
            var l_前行 = new int[p_中間配列.Length + 1];
            var l_今行 = new int[p_中間配列.Length + 1];
            for (var j = 0; j <= p_中間配列.Length; j++)
            {
                l_前行[j] = j;
            }
            for (var i = 1; i <= p_先頭配列.Length; i++)
            {
                l_今行[0] = i;
                for (var j = 1; j <= p_中間配列.Length; j++)
                {
                    var l_コスト = p_先頭配列[i - 1] == p_中間配列[j - 1] ? 0 : 1;
                    l_今行[j] = Math.Min(Math.Min(l_今行[j - 1] + 1, l_前行[j] + 1), l_前行[j - 1] + l_コスト);
                }
                (l_前行, l_今行) = (l_今行, l_前行);
            }
            return l_前行[p_中間配列.Length];
        }

        #endregion
    }
}
