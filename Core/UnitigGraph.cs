using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// unitig 間の隣接関係を、リードマッピングからの推測ではなく
    /// de Bruijn グラフそのものから厳密に構築する。
    ///
    /// これを導入した経緯(実データでの計測):
    /// 従来の UniteContigs は、リード上で連続して観測された unitig ペア
    /// (kmerPath)を隣接候補とし、結合時に「k-1 塩基のオーバーラップ」を
    /// 試し、一致しなければ長い方から順に任意長のオーバーラップを探す
    /// フォールバックを持っていた。しかし実データでは k-1(=30)での一致が
    /// ほぼ起こらず、フォールバックが平均 2.96 塩基という偶然の一致で
    /// unitig を接着していた(2135 箇所)。つまり結合のほぼ全てが誤結合で
    /// あり、かつ contig 総長は unitig 総長のちょうど 2.009 倍に膨れていた
    /// (順鎖・逆鎖の両方が別々の contig として出力されていた)。
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
        /// <summary>頂点数。unitigList と同じ長さ(添字 0,1 は未使用のダミー)。</summary>
        public int VertexCount => this.OutEdges.Count;

        /// <summary>頂点ごとの出辺(行き先の頂点インデックス)。</summary>
        public List<List<int>> OutEdges { get; }

        private UnitigGraph(List<List<int>> outEdges)
        {
            this.OutEdges = outEdges;
        }

        /// <summary>頂点 v の入次数。辺の逆鎖対称性より、v の入次数は v^1 の出次数に等しい。</summary>
        public int InDegree(int vertex)
        {
            return this.OutEdges[vertex ^ 1].Count;
        }

        /// <summary>
        /// unitigList(添字 2u=順鎖, 2u+1=逆鎖の配列)と、k-mer から
        /// (符号付き unitig ID, その向きでの開始位置) への辞書から、
        /// 厳密な隣接グラフを構築する。
        ///
        /// 「先頭 k-mer である(position==0)」ことを要求するのが要点で、
        /// これにより結合が必ず k-1 オーバーラップの単純連結になる。
        /// 複数 unitig に跨る曖昧 k-mer(ambiguousKmerSentinel)は
        /// 行き先を一意に決められないため辺を張らない。
        /// </summary>
        public static UnitigGraph Build(
            List<string> unitigList,
            IReadOnlyDictionary<KmerKey, (int UnitigId, int Position)> kmerDict,
            int kmerLength,
            int ambiguousKmerSentinel)
        {
            List<List<int>> outEdges = [];
            for (var i = 0; i < unitigList.Count; i++)
            {
                outEdges.Add([]);
            }

            // 末尾 k-mer から 1 塩基伸ばした候補を組み立てるための作業バッファ。
            var candidate = new byte[kmerLength];

            for (var vertex = 2; vertex < unitigList.Count; vertex++)
            {
                var seq = unitigList[vertex];
                if (seq.Length < kmerLength)
                {
                    continue;
                }

                // 末尾 k-mer の 2 文字目以降(k-1 塩基)を候補の先頭に置く。
                var tailStart = seq.Length - kmerLength + 1;
                var hasInvalidBase = false;
                for (var i = 0; i < kmerLength - 1; i++)
                {
                    var id = Util.GetSimpleNucleotideID(seq[tailStart + i]);
                    if (id is < Consts.NucleotideID.A or > Consts.NucleotideID.T)
                    {
                        hasInvalidBase = true;
                        break;
                    }
                    candidate[i] = id;
                }
                if (hasInvalidBase)
                {
                    continue;
                }

                for (byte last = Consts.NucleotideID.A; last <= Consts.NucleotideID.T; last++)
                {
                    candidate[kmerLength - 1] = last;
                    if (!kmerDict.TryGetValue(new KmerKey(candidate.AsSpan()), out var hit))
                    {
                        continue;
                    }
                    if (hit.UnitigId == ambiguousKmerSentinel || hit.Position != 0)
                    {
                        // Position != 0 は「その k-mer が unitig の途中に現れる」
                        // ことを意味し、そこへ k-1 オーバーラップで連結することは
                        // できない(unitig 分割が正しければ本来起きないが、
                        // グラフ簡略化で k-mer を削った結果として起こりうる)。
                        continue;
                    }
                    var target = ContigMaker.VertexIndex(hit.UnitigId);
                    if (target == vertex)
                    {
                        // 自己ループは辿ると無限に伸びるため辺として持たない。
                        continue;
                    }
                    outEdges[vertex].Add(target);
                }
            }

            return new UnitigGraph(outEdges);
        }

        /// <summary>
        /// 辺 v→w を、その逆鎖側の双子 w^1→v^1 と対にして取り除く。
        /// 片方だけ消すとグラフの逆鎖対称性が崩れ、順鎖側と逆鎖側で
        /// 別々の経路が組まれてしまう。
        /// </summary>
        private void RemoveEdgePair(int from, int to)
        {
            _ = this.OutEdges[from].Remove(to);
            _ = this.OutEdges[to ^ 1].Remove(from ^ 1);
        }

        /// <summary>辺 v→w を、その逆鎖側の双子 w^1→v^1 と対にして追加する。</summary>
        private void AddEdgePair(int from, int to)
        {
            this.OutEdges[from].Add(to);
            this.OutEdges[to ^ 1].Add(from ^ 1);
        }

        /// <summary>
        /// 短い反復配列を、ペアエンドの証拠に基づいて「通り抜ける経路ごとに複製」して
        /// 解きほぐす(repeat resolution)。
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
        /// <param name="unitigList">複製した配列を追加するため書き換える。</param>
        /// <param name="support">複製した頂点にも元の辺の支持を引き継がせるため書き換える。</param>
        /// <param name="pairLink">頂点対 (v, w) を跨いだフラグメントの本数。</param>
        /// <param name="maxRepeatLength">これより長い unitig は跨げる見込みが無いので対象外。</param>
        /// <returns>解きほぐした反復の数。</returns>
        public int ResolveShortRepeats(
            List<string> unitigList,
            Dictionary<(int, int), ulong> support,
            IReadOnlyDictionary<(int, int), ulong> pairLink,
            int maxRepeatLength,
            decimal uniteThreshold,
            ulong countThreshold)
        {
            var resolved = 0;
            // 複製で頂点が増えるが、増えた分(複製そのもの)は対象にしない。
            var originalVertexCount = this.VertexCount;

            for (var repeat = 2; repeat < originalVertexCount; repeat += 2)
            {
                if (unitigList[repeat].Length > maxRepeatLength)
                {
                    continue;
                }

                var outs = this.OutEdges[repeat];
                var insTwins = this.OutEdges[repeat ^ 1];
                if (outs.Count != 2 || insTwins.Count != 2)
                {
                    continue;
                }

                // repeat へ入ってくる頂点は、双子の出辺の双子。
                var a = insTwins[0] ^ 1;
                var b = insTwins[1] ^ 1;
                var c = outs[0];
                var d = outs[1];

                // 同じ unitig が複数の役回りで現れる退化したケース(自己反復など)は
                // 付け替えの意味が定まらないため触らない。
                int[] involved = [a >> 1, b >> 1, c >> 1, d >> 1, repeat >> 1];
                if (involved.Distinct().Count() != involved.Length)
                {
                    continue;
                }

                var straight = pairLink.GetValueOrDefault((a, c)) + pairLink.GetValueOrDefault((b, d));
                var crossed = pairLink.GetValueOrDefault((a, d)) + pairLink.GetValueOrDefault((b, c));
                var total = straight + crossed;
                if (total < countThreshold)
                {
                    continue;
                }

                var best = Math.Max(straight, crossed);
                if ((decimal)best / total < uniteThreshold)
                {
                    // どちらの対応付けとも決めきれない。無理に繋がない。
                    continue;
                }

                // 勝った対応付けのうち片方を元の repeat に残し、もう片方を複製へ移す。
                var (moveIn, moveOut) = straight >= crossed ? (b, d) : (b, c);

                var duplicate = unitigList.Count; // 常に偶数 = 順鎖側の頂点
                unitigList.Add(unitigList[repeat]);
                unitigList.Add(unitigList[repeat ^ 1]);
                this.OutEdges.Add([]);
                this.OutEdges.Add([]);

                var inSupport = support.GetValueOrDefault((moveIn, repeat));
                var outSupport = support.GetValueOrDefault((repeat, moveOut));

                this.RemoveEdgePair(moveIn, repeat);
                this.RemoveEdgePair(repeat, moveOut);
                this.AddEdgePair(moveIn, duplicate);
                this.AddEdgePair(duplicate, moveOut);

                // 付け替えた辺の支持を複製側へ引き継ぐ(逆鎖側も対称に)。
                support[(moveIn, duplicate)] = inSupport;
                support[(duplicate ^ 1, moveIn ^ 1)] = inSupport;
                support[(duplicate, moveOut)] = outSupport;
                support[(moveOut ^ 1, duplicate ^ 1)] = outSupport;

                resolved++;
            }

            return resolved;
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
        /// maxLengthRatio を超える組は対象外とする。
        /// </summary>
        /// <returns>取り除いた枝の数。</returns>
        public int PopSimpleBubbles(
            List<string> unitigList,
            IReadOnlyDictionary<(int, int), ulong> support,
            double maxLengthRatio = 1.5)
        {
            var popped = 0;

            for (var u = 2; u < this.VertexCount; u++)
            {
                var outs = this.OutEdges[u];
                if (outs.Count < 2)
                {
                    continue;
                }

                // 「1本の unitig を経て同じ頂点へ再合流する」枝を、
                // その再合流先ごとにまとめる。
                Dictionary<int, List<int>> byMergePoint = [];
                foreach (var branch in outs)
                {
                    // 枝は u からのみ入られ、1箇所へのみ出て行く単純な中継でなければならない。
                    if (this.OutEdges[branch].Count != 1 || this.InDegree(branch) != 1)
                    {
                        continue;
                    }
                    var mergePoint = this.OutEdges[branch][0];
                    if (mergePoint == u || (mergePoint >> 1) == (branch >> 1))
                    {
                        continue;
                    }
                    if (!byMergePoint.TryGetValue(mergePoint, out var list))
                    {
                        list = [];
                        byMergePoint[mergePoint] = list;
                    }
                    list.Add(branch);
                }

                foreach (var (mergePoint, branches) in byMergePoint)
                {
                    if (branches.Count < 2)
                    {
                        continue;
                    }

                    var shortest = branches.Min(b => unitigList[b].Length);
                    var longest = branches.Max(b => unitigList[b].Length);
                    if (shortest <= 0 || (double)longest / shortest > maxLengthRatio)
                    {
                        // 長さが揃っていない = 同じ領域の別表現ではなく
                        // 本物の分岐の可能性が高い。触らない。
                        continue;
                    }

                    var winner = branches
                        .OrderByDescending(b => support.GetValueOrDefault((u, b)))
                        .ThenByDescending(b => unitigList[b].Length)
                        .First();

                    foreach (var loser in branches)
                    {
                        if (loser == winner)
                        {
                            continue;
                        }
                        this.RemoveEdgePair(u, loser);
                        this.RemoveEdgePair(loser, mergePoint);
                        popped++;
                    }
                }
            }

            return popped;
        }
    }
}
