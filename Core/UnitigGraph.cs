using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// unitig 間の隣接を、リードマッピングからの推測ではなく de Bruijn グラフ
    /// そのものから構築する。
    ///
    /// 隣接の正しい根拠は「unitig A の末尾 k-mer を1塩基伸ばした k-mer が
    /// unitig B の先頭 k-mer に一致する」ことだけであり、それを満たす辺は
    /// 定義上ちょうど k-1 のオーバーラップを持つ。結合時にオーバーラップ長を
    /// 探索する必要がそもそも無くなる。
    ///
    /// 頂点は向き付き: unitig ID u に対し 2u(順鎖)と 2u+1(逆鎖)、
    /// v の双子は v^1。構築方法より、辺 v→w があれば必ず w^1→v^1 もある。
    /// </summary>
    internal sealed class UnitigGraph
    {
        /// <summary>頂点ごとの出辺(行き先の頂点インデックス)。</summary>
        public List<List<int>> A_出辺 { get; }

        private UnitigGraph(List<List<int>> p_出辺)
        {
            this.A_出辺 = p_出辺;
        }

        /// <summary>頂点の入次数。辺の逆鎖対称性より、v の入次数は v^1 の出次数に等しい。</summary>
        public int Get_入次数(int p_頂点)
        {
            return this.A_出辺[p_頂点 ^ 1].Count;
        }

        /// <summary>
        /// 隣接グラフを構築する。
        ///
        /// 行き先が「先頭 k-mer である(開始位置==0)」ことを要求するのが要点で、
        /// これにより結合が必ず k-1 オーバーラップの単純連結になる。
        /// 曖昧 k-mer は行き先を一意に決められないため辺を張らない。
        /// </summary>
        public static UnitigGraph Get_グラフ(
            List<string> p_ユニティグ配列,
            IReadOnlyDictionary<KmerKey, (int A_ユニティグID, int A_開始位置)> p_kmer辞書,
            int p_k長,
            int p_曖昧kmerの番兵)
        {
            List<List<int>> l_出辺 = [];
            for (var i = 0; i < p_ユニティグ配列.Count; i++)
            {
                l_出辺.Add([]);
            }

            // 末尾 k-mer から 1 塩基伸ばした候補を組み立てるための作業バッファ。
            var l_候補 = new byte[p_k長];

            for (var l_頂点 = 2; l_頂点 < p_ユニティグ配列.Count; l_頂点++)
            {
                var l_配列 = p_ユニティグ配列[l_頂点];
                if (l_配列.Length < p_k長)
                {
                    continue;
                }

                // 末尾 k-mer の 2 文字目以降(k-1 塩基)を候補の先頭に置く。
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

                for (byte l_末尾塩基 = Consts.塩基ID.A; l_末尾塩基 <= Consts.塩基ID.T; l_末尾塩基++)
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
                        // できない(unitig 分割が正しければ本来起きないが、
                        // グラフ簡略化で k-mer を削った結果として起こりうる)。
                        continue;
                    }
                    var l_行き先 = ContigMaker.Get_頂点番号(l_ヒット.A_ユニティグID);
                    if (l_行き先 == l_頂点)
                    {
                        // 自己ループは辿ると無限に伸びるため辺として持たない。
                        continue;
                    }
                    l_出辺[l_頂点].Add(l_行き先);
                }
            }

            return new UnitigGraph(l_出辺);
        }

        /// <summary>
        /// 辺 v→w を、その逆鎖側の双子 w^1→v^1 と対にして取り除く。
        /// 片方だけ消すとグラフの逆鎖対称性が崩れ、順鎖側と逆鎖側で
        /// 別々の経路が組まれてしまう。
        /// </summary>
        private void V_除去_辺の対(int p_始点, int p_終点)
        {
            _ = this.A_出辺[p_始点].Remove(p_終点);
            _ = this.A_出辺[p_終点 ^ 1].Remove(p_始点 ^ 1);
        }

        /// <summary>辺 v→w を、その逆鎖側の双子 w^1→v^1 と対にして追加する。</summary>
        private void V_追加_辺の対(int p_始点, int p_終点)
        {
            this.A_出辺[p_始点].Add(p_終点);
            this.A_出辺[p_終点 ^ 1].Add(p_始点 ^ 1);
        }

        /// <summary>
        /// 短い反復配列を、ペアエンドの証拠に基づいて経路ごとに複製して解きほぐす。
        ///
        /// 反復 R が a→R→c と b→R→d の文脈を持つとき、グラフ上では R が
        /// 1頂点に潰れて入次数2・出次数2になる。R 内部のリードはどちらのコピー
        /// 由来か区別できないため、リード支持では原理的に解けない。
        /// 解ける唯一の手がかりは R を丸ごと跨いだフラグメントで、a-c と b-d の
        /// ペアが多く a-d / b-c に乏しければ対応が決まる。
        ///
        /// R を複製して片方の経路を付け替えると、どちらも入次数1・出次数1の
        /// 一本道になり既存の walk がそのまま伸ばせる。R の配列が2回出力されるのは
        /// 実際に2回現れることの反映であって水増しではない。
        /// </summary>
        /// <returns>解きほぐした反復の数。</returns>
        public int V_解決_短い反復(
            List<string> p_ユニティグ配列,
            Dictionary<(int, int), ulong> p_支持,
            IReadOnlyDictionary<(int, int), ulong> p_ペア連結,
            int p_反復長の上限,
            decimal p_優勢閾値,
            ulong p_最小証拠数)
        {
            var l_解決数 = 0;
            // 複製で頂点が増えるが、増えた分(複製そのもの)は対象にしない。
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

                // 反復へ入ってくる頂点は、双子の出辺の双子。
                var l_入1 = l_入辺の双子[0] ^ 1;
                var l_入2 = l_入辺の双子[1] ^ 1;
                var l_出1 = l_出辺[0];
                var l_出2 = l_出辺[1];

                // 同じ unitig が複数の役回りで現れる退化したケース(自己反復など)は
                // 付け替えの意味が定まらないため触らない。
                int[] l_関係する頂点 = [l_入1 >> 1, l_入2 >> 1, l_出1 >> 1, l_出2 >> 1, l_反復頂点 >> 1];
                if (l_関係する頂点.Distinct().Count() != l_関係する頂点.Length)
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
                    // どちらの対応付けとも決めきれない。無理に繋がない。
                    continue;
                }

                // 勝った対応付けのうち片方を元の反復頂点に残し、もう片方を複製へ移す。
                var (l_移す入辺, l_移す出辺) = l_平行 >= l_交差 ? (l_入2, l_出2) : (l_入2, l_出1);

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

                // 付け替えた辺の支持を複製側へ引き継ぐ(逆鎖側も対称に)。
                p_支持[(l_移す入辺, l_複製)] = l_入辺の支持;
                p_支持[(l_複製 ^ 1, l_移す入辺 ^ 1)] = l_入辺の支持;
                p_支持[(l_複製, l_移す出辺)] = l_出辺の支持;
                p_支持[(l_移す出辺 ^ 1, l_複製 ^ 1)] = l_出辺の支持;

                l_解決数++;
            }

            return l_解決数;
        }

        /// <summary>
        /// 単純バブル(u から分かれた枝が1本の unitig を経て同じ w へ再合流する構造)を
        /// 検出し、リード支持が最も高い枝以外の辺を取り除く。
        ///
        /// 相互一意を結合の条件にしているため、バブルがあると再合流点の入次数が
        /// 2以上のままになり、その経路全体が結合されなくなる。半数体である
        /// 細菌ゲノムにバブルは本来存在しない(エラーか株レベルの変異)。
        ///
        /// 敗者の配列自体は削除しない。誤りだった場合の損害が大きく、辺だけ外せば
        /// 単独 contig として出力されるので内容は失われない。
        /// 長さが大きく異なる枝はバブルではなく本物の分岐の可能性が高いため除く。
        /// </summary>
        /// <returns>取り除いた枝の数。</returns>
        public int V_除去_単純バブル(
            List<string> p_ユニティグ配列,
            IReadOnlyDictionary<(int, int), ulong> p_支持,
            double p_長さ比の上限 = 1.5)
        {
            var l_除去数 = 0;

            for (var l_分岐元 = 2; l_分岐元 < this.A_出辺.Count; l_分岐元++)
            {
                var l_出辺 = this.A_出辺[l_分岐元];
                if (l_出辺.Count < 2)
                {
                    continue;
                }

                // 「1本の unitig を経て同じ頂点へ再合流する」枝を、
                // その再合流先ごとにまとめる。
                Dictionary<int, List<int>> l_再合流先ごと = [];
                foreach (var l_枝 in l_出辺)
                {
                    // 枝は分岐元からのみ入られ、1箇所へのみ出て行く単純な中継でなければならない。
                    if (this.A_出辺[l_枝].Count != 1 || this.Get_入次数(l_枝) != 1)
                    {
                        continue;
                    }
                    var l_再合流先 = this.A_出辺[l_枝][0];
                    if (l_再合流先 == l_分岐元 || (l_再合流先 >> 1) == (l_枝 >> 1))
                    {
                        continue;
                    }
                    if (!l_再合流先ごと.TryGetValue(l_再合流先, out var l_枝一覧))
                    {
                        l_枝一覧 = [];
                        l_再合流先ごと[l_再合流先] = l_枝一覧;
                    }
                    l_枝一覧.Add(l_枝);
                }

                foreach (var (l_再合流先, l_枝群) in l_再合流先ごと)
                {
                    if (l_枝群.Count < 2)
                    {
                        continue;
                    }

                    var l_最短 = l_枝群.Min(x => p_ユニティグ配列[x].Length);
                    var l_最長 = l_枝群.Max(x => p_ユニティグ配列[x].Length);
                    if (l_最短 <= 0 || (double)l_最長 / l_最短 > p_長さ比の上限)
                    {
                        // 長さが揃っていない = 同じ領域の別表現ではなく
                        // 本物の分岐の可能性が高い。触らない。
                        continue;
                    }

                    var l_勝者 = l_枝群
                        .OrderByDescending(x => p_支持.GetValueOrDefault((l_分岐元, x)))
                        .ThenByDescending(x => p_ユニティグ配列[x].Length)
                        .First();

                    foreach (var l_敗者 in l_枝群)
                    {
                        if (l_敗者 == l_勝者)
                        {
                            continue;
                        }
                        this.V_除去_辺の対(l_分岐元, l_敗者);
                        this.V_除去_辺の対(l_敗者, l_再合流先);
                        l_除去数++;
                    }
                }
            }

            return l_除去数;
        }
    }
}
