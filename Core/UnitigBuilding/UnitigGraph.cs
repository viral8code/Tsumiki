using System.Text;
using Tsumiki.Common;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Core.UnitigBuilding
{
    /// <summary>
    /// unitig 間の隣接を、リードマッピングからの推測ではなく de Bruijn グラフ
    /// そのものから構築する<br/>
    /// 隣接の正しい根拠は「unitig A の末尾 k-mer を 1 塩基伸ばした k-mer が
    /// unitig B の先頭 k-mer に一致する」ことだけであり、それを満たす辺は
    /// 定義上ちょうど k-1 のオーバーラップを持つ<br/>
    /// 結合時にオーバーラップ長を
    /// 探索する必要がそもそも無くなる<br/>
    /// 頂点は向き付き: unitig ID u に対し 2 u(順鎖) と 2 u+1(逆鎖)、
    /// v の双子は v^1<br/>
    /// 構築方法より、辺 v→w があれば必ず w^1→v^1 もある
    /// </summary>
    internal sealed class UnitigGraph
    {
        /// <summary>
        /// 頂点ごとの出辺 (行き先の頂点インデックス)
        /// </summary>
        public List<List<int>> A_出辺 { get; }

        /// <summary>
        /// 末尾を 1 塩基伸ばすと自分の先頭 k-mer に戻る頂点<br/>
        /// 辺としては持てない (辿ると伸び続ける) が、分岐を持たない環状の
        /// 複製単位はこの形でしか現れないため、事実だけは残しておく
        /// </summary>
        public HashSet<int> A_自己ループ { get; }

        private UnitigGraph(List<List<int>> p_出辺, HashSet<int> p_自己ループ)
        {
            this.A_出辺 = p_出辺;
            this.A_自己ループ = p_自己ループ;
        }

        /// <summary>
        /// 頂点の入次数<br/>
        /// 辺の逆鎖対称性より、v の入次数は v^1 の出次数に等しい
        /// </summary>
        public int Get_入次数(int p_頂点)
        {
            return this.A_出辺[p_頂点 ^ 1].Count;
        }

        /// <summary>
        /// その頂点を通り抜けてよいか<br/>
        /// A-R-B-R-C(R は 2 コピーの反復) で A→R と R→C はどちらも本物の隣接だが、
        /// walk は各 unitig を 1 回しか使えないため、連鎖させると中間の B を
        /// 飛ばした A-R-C ができてしまう<br/>
        /// 通り抜けてよいのは反復が解きほぐされ
        /// 入次数・出次数がどちらも 1 になった、どのコピーにいるか確定した状態だけ
        /// </summary>
        public bool Get_通り抜けてよいか(IReadOnlyDictionary<int, int>? p_コピー数, int p_頂点)
        {
            return (p_コピー数?.GetValueOrDefault(p_頂点 >> 1, 1) ?? 1) <= 1 || this.A_出辺[p_頂点].Count == 1 && this.Get_入次数(p_頂点) == 1;
        }

        /// <summary>
        /// 隣接グラフを構築する<br/>
        /// 行き先が「先頭 k-mer である (開始位置==0)」ことを要求するのが要点で、
        /// これにより結合が必ず k-1 オーバーラップの単純連結になる<br/>
        /// 曖昧 k-mer は行き先を一意に決められないため辺を張らない
        /// </summary>
        public static UnitigGraph Get_グラフ(
            List<string> p_ユニティグ配列,
            IReadOnlyDictionary<KmerKey, (int A_ユニティグID, int A_開始位置)> p_kmer辞書,
            int p_k長,
            int p_曖昧kmerの番兵)
        {
            List<List<int>> l_出辺 = [];
            HashSet<int> l_自己ループ = [];
            for (var i = 0; i < p_ユニティグ配列.Count; i++)
            {
                l_出辺.Add([]);
            }

            // 末尾 k-mer から 1 塩基伸ばした候補を組み立てるための作業バッファ
            var l_候補 = new byte[p_k長];

            for (var l_頂点 = 2; l_頂点 < p_ユニティグ配列.Count; l_頂点++)
            {
                var l_配列 = p_ユニティグ配列[l_頂点];
                if (l_配列.Length < p_k長)
                {
                    continue;
                }

                // 末尾 k-mer の 2 文字目以降 (k-1 塩基) を候補の先頭に置く
                var l_末尾開始 = l_配列.Length - p_k長 + 1;
                var l_無効な塩基があるか = false;
                for (var i = 0; i < p_k長 - 1; i++)
                {
                    var l_塩基ID = Util.Get_塩基ID(l_配列[l_末尾開始 + i]);
                    if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                    {
                        l_無効な塩基があるか = true;
                        break;
                    }
                    l_候補[i] = l_塩基ID;
                }
                if (l_無効な塩基があるか)
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
                    if (l_ヒット.A_ユニティグID == p_曖昧kmerの番兵 || l_ヒット.A_開始位置 != 0)
                    {
                        // 開始位置 != 0 は「その k-mer が unitig の途中に現れる」
                        // ことを意味し、そこへ k-1 オーバーラップで連結することは
                        // できない (unitig 分割が正しければ本来起きないが、
                        // グラフ簡略化で k-mer を削った結果として起こりうる)
                        continue;
                    }
                    var l_行き先 = ContigMaker.Get_頂点番号(l_ヒット.A_ユニティグID);
                    if (l_行き先 == l_頂点)
                    {
                        // 自己ループは辿ると無限に伸びるため辺として持たない
                        // ただし環状に閉じている根拠そのものなので、事実は残す
                        _ = l_自己ループ.Add(l_頂点);
                        continue;
                    }
                    l_出辺[l_頂点].Add(l_行き先);
                }
            }

            return new UnitigGraph(l_出辺, l_自己ループ);
        }

        /// <summary>
        /// 辺 v→w を、その逆鎖側の双子 w^1→v^1 と対にして取り除く<br/>
        /// 片方だけ消すとグラフの逆鎖対称性が崩れ、順鎖側と逆鎖側で
        /// 別々の経路が組まれてしまう
        /// </summary>
        private void V_除去_辺の対(int p_始点, int p_終点)
        {
            _ = this.A_出辺[p_始点].Remove(p_終点);
            _ = this.A_出辺[p_終点 ^ 1].Remove(p_始点 ^ 1);
        }

        /// <summary>
        /// 辺 v→w を、その逆鎖側の双子 w^1→v^1 と対にして追加する
        /// </summary>
        private void V_追加_辺の対(int p_始点, int p_終点)
        {
            this.A_出辺[p_始点].Add(p_終点);
            this.A_出辺[p_終点 ^ 1].Add(p_始点 ^ 1);
        }

        /// <summary>
        /// 短い反復配列を、ペアエンドの証拠に基づいて経路ごとに複製して解きほぐす<br/>
        /// 反復 R が a→R→c と b→R→d の文脈を持つとき、グラフ上では R が
        /// 1 頂点に潰れて入次数 2・出次数 2 になる<br/>
        /// R 内部のリードはどちらのコピー
        /// 由来か区別できないため、リード支持では原理的に解けない<br/>
        /// 解ける唯一の手がかりは R を丸ごと跨いだフラグメントで、a-c と b-d の
        /// ペアが多く a-d / b-c に乏しければ対応が決まる<br/>
        /// R を複製して片方の経路を付け替えると、どちらも入次数 1・出次数 1 の
        /// 一本道になり既存の walk がそのまま伸ばせる<br/>
        /// R の配列が 2 回出力されるのは
        /// 実際に 2 回現れることの反映であって水増しではない
        /// </summary>
        /// <param name="p_r_mer検証器">
        /// 渡すと、対応付けが確定したあとに複製を確定する前段として、
        /// 提案された 2 本の経路 (勝ったペアリングそれぞれ) を r-mer で検証する
        /// (ABySS RResolver 型の拒否権)<br/>
        /// どちらか一方でも接合点を跨ぐ r-mer の
        /// 支持が足りなければ複製を見送る<br/>
        /// 注意: head-repeat・repeat-tail はどちらの対応付けでも de Bruijn
        /// グラフ上の本物の辺なので、正しい対応付けのリードだけからも
        /// 個々の接合点の存在は独立に確認できてしまう<br/>
        /// したがってこの検証は
        /// 「対応付け A」と「対応付け B」のどちらが正しいかを区別する力は
        /// 本質的に持たない (反復が反復である以上、局所的な文脈だけでは
        /// 区別できないため)<br/>
        /// 実際に効くのは、ペア支持が示す対応付けについて
        /// 個々の接合点すら生リードに一切裏付けられない (=そもそもその
        /// unitig 同士が隣接している根拠が生データに無い) 場合であり、
        /// この限定的だが無視できない安全網として使う
        /// </param>
        /// <returns>解きほぐした反復の数<br/>
        /// </returns>
        public int V_解決_短い反復(
            List<string> p_ユニティグ配列,
            Dictionary<(int, int), ulong> p_支持,
            IReadOnlyDictionary<(int, int), ulong> p_ペア連結,
            int p_反復長の上限,
            decimal p_優勢閾値,
            ulong p_最小証拠数,
            RepeatRMerVerifier? p_r_mer検証器 = null)
        {
            var l_解決数 = 0;
            var l_r_mer検証で棄却した数 = 0;
            // 複製で頂点が増えるが、増えた分 (複製そのもの) は対象にしない
            var l_元の頂点数 = this.A_出辺.Count;

            for (var l_反復頂点 = 2; l_反復頂点 < l_元の頂点数; l_反復頂点 += 2)
            {
                if (p_ユニティグ配列[l_反復頂点].Length > p_反復長の上限)
                {
                    continue;
                }

                var l_出辺 = this.A_出辺[l_反復頂点];
                var l_入辺の双子 = this.A_出辺[l_反復頂点 ^ 1];
                if (l_出辺.Count != 2 || l_入辺の双子.Count != 2)
                {
                    continue;
                }

                // 反復へ入ってくる頂点は、双子の出辺の双子
                var l_入1 = l_入辺の双子[0] ^ 1;
                var l_入2 = l_入辺の双子[1] ^ 1;
                var l_出1 = l_出辺[0];
                var l_出2 = l_出辺[1];

                // 反復自身が周囲に現れる (タンデム反復・自己ループ) 場合と、
                // 入口どうし・出口どうしが同じ unitig の場合 (逆向き反復の
                // ヘアピンなど) は、付け替えの意味が定まらないため触らない
                //
                // 一方、入口と出口に同じ unitig が現れることは退化ではない
                // 環状の複製単位に同じ反復が 2 回現れると、その間に挟まれた
                // 2 つの領域が必ず両方の役回りに立つ
                // 細菌のゲノムで最も
                // ありふれた反復の形なので、ここまで弾くと、証拠が揃っていて
                // 解ける反復が解けないまま残る
                var l_反復のID = l_反復頂点 >> 1;
                int[] l_周囲 = [l_入1 >> 1, l_入2 >> 1, l_出1 >> 1, l_出2 >> 1];
                if (l_周囲.Contains(l_反復のID)
                    || (l_入1 >> 1) == (l_入2 >> 1)
                    || (l_出1 >> 1) == (l_出2 >> 1))
                {
                    continue;
                }

                var l_平行 = p_ペア連結.GetValueOrDefault((l_入1, l_出1)) + p_ペア連結.GetValueOrDefault((l_入2, l_出2));
                var l_交差 = p_ペア連結.GetValueOrDefault((l_入1, l_出2)) + p_ペア連結.GetValueOrDefault((l_入2, l_出1));
                var l_合計 = l_平行 + l_交差;
                if (l_合計 < p_最小証拠数)
                {
                    continue;
                }

                var l_最良 = Math.Max(l_平行, l_交差);
                if ((decimal)l_最良 / l_合計 < p_優勢閾値)
                {
                    // どちらの対応付けとも決めきれない
                    // 無理に繋がない
                    continue;
                }

                // 勝った対応付けのうち片方を元の反復頂点に残し、もう片方を複製へ移す
                var (l_移す入辺, l_移す出辺) = l_平行 >= l_交差 ? (l_入2, l_出2) : (l_入2, l_出1);
                var l_残る出辺 = l_平行 >= l_交差 ? l_出1 : l_出2;

                if (p_r_mer検証器 is not null)
                {
                    // 集計されたペア支持は「跨いだリードが実在するか」を直接
                    // 確かめていない
                    // r-mer で両方の経路を独立に検証し、
                    // どちらか一方でも接合点の支持が足りなければ、この対応付け
                    // 自体を疑って複製しない (誤った複製は取りこぼしではなく
                    // 実在しない配列を作る偽陽性になるため、疑わしきは見送る)
                    var l_残る側で支持あるか = p_r_mer検証器.Get_接合点に支持があるか(
                        p_ユニティグ配列[l_入1], p_ユニティグ配列[l_反復頂点], p_ユニティグ配列[l_残る出辺],
                        Consts.r_mer接合点支持の閾値の既定値);
                    var l_移す側で支持あるか = p_r_mer検証器.Get_接合点に支持があるか(
                        p_ユニティグ配列[l_移す入辺], p_ユニティグ配列[l_反復頂点], p_ユニティグ配列[l_移す出辺],
                        Consts.r_mer接合点支持の閾値の既定値);
                    if (!l_残る側で支持あるか || !l_移す側で支持あるか)
                    {
                        l_r_mer検証で棄却した数++;
                        continue;
                    }
                }

                var l_複製 = p_ユニティグ配列.Count; // 常に偶数 = 順鎖側の頂点
                p_ユニティグ配列.Add(p_ユニティグ配列[l_反復頂点]);
                p_ユニティグ配列.Add(p_ユニティグ配列[l_反復頂点 ^ 1]);
                this.A_出辺.Add([]);
                this.A_出辺.Add([]);

                var l_入辺の支持 = p_支持.GetValueOrDefault((l_移す入辺, l_反復頂点));
                var l_出辺の支持 = p_支持.GetValueOrDefault((l_反復頂点, l_移す出辺));

                this.V_除去_辺の対(l_移す入辺, l_反復頂点);
                this.V_除去_辺の対(l_反復頂点, l_移す出辺);
                this.V_追加_辺の対(l_移す入辺, l_複製);
                this.V_追加_辺の対(l_複製, l_移す出辺);

                // 付け替えた辺の支持を複製側へ引き継ぐ (逆鎖側も対称に)
                p_支持[(l_移す入辺, l_複製)] = l_入辺の支持;
                p_支持[(l_複製 ^ 1, l_移す入辺 ^ 1)] = l_入辺の支持;
                p_支持[(l_複製, l_移す出辺)] = l_出辺の支持;
                p_支持[(l_移す出辺 ^ 1, l_複製 ^ 1)] = l_出辺の支持;

                l_解決数++;
            }

            if (p_r_mer検証器 is not null && l_r_mer検証で棄却した数 > 0)
            {
                Logger.V_出力(メッセージID.rMer検証による棄却, l_r_mer検証で棄却した数);
            }

            return l_解決数;
        }

        /// <summary>
        /// 単純バブル (u から分かれた枝が、途中に本物の分岐の無い 1 本の経路
        /// (1 つ以上の unitig の連なり) を経て同じ w へ再合流する構造) を検出し、
        /// リード支持が最も高い経路以外の辺を取り除く (SPAdes の
        /// AlternativesAnalyzer・MEGAHIT の ComplexBubbleRemover に相当)<br/>
        /// 相互一意を結合の条件にしているため、バブルがあると再合流点の入次数が
        /// 2 以上のままになり、その経路全体が結合されなくなる<br/>
        /// 半数体である
        /// 細菌ゲノムにバブルは本来存在しない (エラーか株レベルの変異)<br/>
        /// 経路は固定の長さ比ではなく delta = max(p_長さ帯の下限,
        /// p_長さ帯の割合 * 最短経路長) の帯で比較する<br/>
        /// 中間に unitig を複数
        /// 挟む経路や、長さがぴったり揃わない経路も対象になる<br/>
        /// 長さが揃っていても
        /// 配列が大きく異なる経路 (たまたま長さが一致した別の反復など) は、
        /// 編集距離ベースの類似度 (p_類似度の下限)で弾く<br/>
        /// 敗者の経路自体は削除しない<br/>
        /// 誤りだった場合の損害が大きく、辺だけ外せば
        /// 単独 contig として出力されるので内容は失われない<br/>
        /// p_敗者への引き継ぎ先 を
        /// 渡すと、敗者の配列 (MEGAHIT の careful_bubble)をそこへ集める<br/>
        /// 「この k では敗者と判断したが、それは決定であって事実ではない<br/>
        /// 次の k は自分の証拠で判断し直せる」という KmerCarryOver と同じ思想
        /// </summary>
        /// <returns>取り除いた経路の数<br/>
        /// </returns>
        public int V_除去_単純バブル(
            List<string> p_ユニティグ配列,
            IReadOnlyDictionary<(int, int), ulong> p_支持,
            int p_k長,
            List<string>? p_敗者への引き継ぎ先 = null,
            double p_長さ帯の割合 = 0.1D,
            int p_長さ帯の下限 = 3,
            double p_類似度の下限 = 0.7D,
            int p_経路長の上限 = 2000)
        {
            var l_除去数 = 0;

            for (var l_分岐元 = 2; l_分岐元 < this.A_出辺.Count; l_分岐元++)
            {
                var l_出辺 = this.A_出辺[l_分岐元];
                if (l_出辺.Count < 2)
                {
                    continue;
                }

                // 「途中に本物の分岐の無い経路を経て同じ頂点へ再合流する」枝を、
                // その再合流先ごとにまとめる
                Dictionary<int, List<List<int>>> l_再合流先ごと = [];
                foreach (var l_開始 in l_出辺)
                {
                    if (this.Get_単純経路(l_開始, p_ユニティグ配列, p_k長, p_経路長の上限) is not { } l_結果)
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

                    var l_配列群 = l_経路群.Select(x => Get_経路配列(p_ユニティグ配列, x, p_k長)).ToList();

                    var l_基準長 = l_配列群.Min(x => x.Length);
                    var l_delta = Math.Max(p_長さ帯の下限, p_長さ帯の割合 * l_基準長);
                    if (l_配列群.Any(x => Math.Abs(x.Length - l_基準長) > l_delta))
                    {
                        // 長さが揃っていない = 同じ領域の別表現ではなく
                        // 本物の分岐の可能性が高い
                        // 触らない
                        continue;
                    }

                    var l_基準配列 = l_配列群[0];
                    if (l_配列群.Skip(1).Any(x => Get_類似度(l_基準配列, x) < p_類似度の下限))
                    {
                        // 長さは近いが配列がまるで違う
                        // 同じ領域の別表現とは言えない
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
                        this.V_除去_辺の対(l_分岐元, l_敗者経路[0]);
                        for (var j = 0; j + 1 < l_敗者経路.Count; j++)
                        {
                            this.V_除去_辺の対(l_敗者経路[j], l_敗者経路[j + 1]);
                        }
                        this.V_除去_辺の対(l_敗者経路[^1], l_再合流先);
                        l_除去数++;
                        p_敗者への引き継ぎ先?.Add(l_配列群[i]);
                    }
                }
            }

            return l_除去数;
        }

        /// <summary>
        /// p_開始 から、途中に本物の分岐が無い限り辿れるだけ辿った経路と、
        /// その先の再合流先 (=最初に他からも入ってくる頂点) を返す<br/>
        /// 判定できない (開始点が既に他からも入られている、途中で行き止まる/
        /// さらに分岐する、循環する、長さの上限を超える) 場合は null
        /// </summary>
        private (List<int> A_経路, int A_再合流先)? Get_単純経路(
            int p_開始, List<string> p_ユニティグ配列, int p_k長, int p_長さ上限)
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
                    // 開始点以外から見て他からも入ってくる (=本物の再合流点)、
                    // あるいは循環に突入した
                    // 前者かつ経路が空でなければ、
                    // この頂点そのものが再合流先
                    return l_経路.Count > 0 ? (l_経路, l_現在) : null;
                }

                var l_この頂点の長さ = p_ユニティグ配列[l_現在].Length;
                l_累積長 += l_経路.Count == 0 ? l_この頂点の長さ : Math.Max(0, l_この頂点の長さ - l_重なり長);
                if (l_累積長 > p_長さ上限)
                {
                    return null;
                }

                l_経路.Add(l_現在);

                var l_出辺 = this.A_出辺[l_現在];
                if (l_出辺.Count != 1)
                {
                    // 行き止まり、または途中でさらに分岐している
                    // 単純な経路ではない
                    return null;
                }
                l_現在 = l_出辺[0];
            }
        }

        private static string Get_経路配列(List<string> p_ユニティグ配列, List<int> p_経路, int p_k長)
        {
            var l_重なり長 = p_k長 - 1;
            var l_出力 = new StringBuilder(p_ユニティグ配列[p_経路[0]]);
            for (var i = 1; i < p_経路.Count; i++)
            {
                var l_配列 = p_ユニティグ配列[p_経路[i]];
                _ = l_出力.Append(l_配列.Length > l_重なり長 ? l_配列[l_重なり長..] : string.Empty);
            }
            return l_出力.ToString();
        }

        private static double Get_類似度(string p_a, string p_b)
        {
            var l_最大長 = Math.Max(p_a.Length, p_b.Length);
            return l_最大長 == 0 ? 1.0 : 1.0 - ((double)Get_編集距離(p_a, p_b) / l_最大長);
        }

        private static int Get_編集距離(string p_a, string p_b)
        {
            var l_前行 = new int[p_b.Length + 1];
            var l_今行 = new int[p_b.Length + 1];
            for (var j = 0; j <= p_b.Length; j++)
            {
                l_前行[j] = j;
            }
            for (var i = 1; i <= p_a.Length; i++)
            {
                l_今行[0] = i;
                for (var j = 1; j <= p_b.Length; j++)
                {
                    var l_コスト = p_a[i - 1] == p_b[j - 1] ? 0 : 1;
                    l_今行[j] = Math.Min(Math.Min(l_今行[j - 1] + 1, l_前行[j] + 1), l_前行[j - 1] + l_コスト);
                }
                (l_前行, l_今行) = (l_今行, l_前行);
            }
            return l_前行[p_b.Length];
        }
    }
}
