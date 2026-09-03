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
