using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// unitig 間の隣接関係を、リードマッピングからの推測ではなく
    /// de Bruijn グラフそのものから厳密に構築する。
    ///
    /// これを導入した経緯(実データでの計測):
    /// 従来の contig 結合は、リード上で連続して観測された unitig ペアを
    /// 隣接候補とし、結合時に「k-1 塩基のオーバーラップ」を試し、一致しなければ
    /// 長い方から順に任意長のオーバーラップを探すフォールバックを持っていた。
    /// しかし実データでは k-1(=30)での一致がほぼ起こらず、フォールバックが
    /// 平均 2.96 塩基という偶然の一致で unitig を接着していた(2135 箇所)。
    /// つまり結合のほぼ全てが誤結合であり、かつ contig 総長は unitig 総長の
    /// ちょうど 2.009 倍に膨れていた(順鎖・逆鎖の両方が別々の contig として
    /// 出力されていた)。
    ///
    /// 隣接の唯一の正しい根拠は de Bruijn グラフの辺であり、それは
    /// 「unitig A の末尾 k-mer から 1 塩基伸ばした k-mer が unitig B の
    /// 先頭 k-mer に一致する」ことと同値である。この条件を満たす辺は
    /// 定義上ちょうど k-1 塩基のオーバーラップを持つため、結合時に
    /// オーバーラップ長を探索する必要がそもそも無くなる。
    ///
    /// 頂点は「符号付き向き」を持つ: unitig ID u に対し
    /// 頂点 2u(順鎖)と 2u+1(逆鎖)。ある頂点 v の逆鎖側の双子は v^1。
    /// 本クラスの構築方法により、辺 v→w が存在すれば必ず w^1→v^1 も
    /// 存在する(逆相補を取れば同じ重なりが成立するため)。
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
        /// unitig 配列一覧(添字 2u=順鎖, 2u+1=逆鎖)と、k-mer から
        /// (符号付き unitig ID, その向きでの開始位置) への辞書から、
        /// 厳密な隣接グラフを構築する。
        ///
        /// 「先頭 k-mer である(開始位置==0)」ことを要求するのが要点で、
        /// これにより結合が必ず k-1 オーバーラップの単純連結になる。
        /// 複数 unitig に跨る曖昧 k-mer は行き先を一意に決められないため辺を張らない。
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
        /// 短い反復配列を、ペアエンドの証拠に基づいて「通り抜ける経路ごとに複製」して
        /// 解きほぐす。
        ///
        /// なぜ必要か: 反復配列 R がゲノム中に2回現れ、それぞれ a→R→c と b→R→d と
        /// いう文脈を持つ場合、de Bruijn グラフ上では R は1個の頂点に潰れ、
        /// 入次数2・出次数2になる。R の内部から読まれたリードはどちらのコピー由来か
        /// 区別できないため、R→c と R→d はどちらも同程度の支持を得てしまい、
        /// リード支持による分岐解決では原理的に解けない(比率が5割前後になり、
        /// 「根拠がないので繋がない」と正しく判断されて経路が途切れる)。
        ///
        /// 解ける唯一の手がかりは「R を丸ごと跨いだフラグメント」である。
        /// 片端が a に、もう片端が c に載ったペアが多数あり、かつ a-d / b-c の
        /// 組み合わせには乏しいなら、対応は (a,c) と (b,d) だと判断できる。
        ///
        /// 判断がついたら R を複製し、片方の経路を複製側へ付け替える。こうすると
        /// どちらの経路も入次数1・出次数1の一本道になり、既存の相互一意性の検査と
        /// walk がそのまま両方を伸ばせる。複製によって R の配列は2回出力されるが、
        /// これは「実際にゲノム中に2回現れる」ことの正しい反映であって水増しではない。
        ///
        /// 実データ(k=63)ではこの形(入次数2かつ出次数2)の unitig が151本あり、
        /// うち143本がフラグメント長(中央値245bp)より短く、跨げる見込みがあった。
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
        /// 単純バブル(ある頂点 u から分かれた複数の枝が、それぞれ1本の
        /// unitig を経て同じ頂点 w へ再合流する構造)を検出し、
        /// リード支持が最も高い枝だけを残して他の枝の辺を取り除く。
        ///
        /// これが必要な理由: 結合の採用条件を相互一意(vの唯一の行き先がw、
        /// かつwの唯一の来訪元がv)にしたため、バブルがあると再合流点 w の
        /// 入次数が2以上のままになり、u から w へ至る経路全体が一切
        /// 結合されなくなる。バブルは半数体である細菌ゲノムでは本来
        /// 存在しないはず(シーケンスエラーか株レベルの変異)なので、
        /// 支持の低い枝を経路から外すのが妥当。
        ///
        /// ただし敗者の unitig 配列そのものは削除しない。実配列を消すのは
        /// 誤りだった場合の損害が大きく、辺だけ外せば孤立した単独 contig
        /// として出力されるため内容は失われない。
        ///
        /// 長さが大きく異なる枝は「バブル」ではなく本物の分岐(反復配列の
        /// 出入口など)である可能性が高いため、長さ比が
        /// p_長さ比の上限 を超える組は対象外とする。
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
