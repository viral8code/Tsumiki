using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker.UniteContigs で確定した contigs.fasta を読み直し、
    /// ペアエンド由来の隣接情報(ContigMaker.PairPath)を使って
    /// contig 同士をさらに N 埋めで連結する(スキャフォールディング)。
    ///
    /// 実行タイミング: UniteContigs 完了後、Program.cs から別処理として呼び出す。
    /// 入力: 確定済み contigs.fasta ファイル(読み直す)。
    /// 出力: 新規ファイル scaffolds.fasta。contigs.fasta 自体は変更しない。
    /// </summary>
    internal class Scaffolder(ContigMaker contigMaker, string contigFilePath)
    {

        // contig ID(FastaWriter が振った 1 始まりの ID) -> 配列本体。
        private readonly Dictionary<int, string> contigSequences = [];

        // contig ID -> ID 文字列(先頭 ">" の次に書かれていた文字列。"NODE1" 等)。
        // 出力時に元の命名をある程度踏襲するために保持する。
        private readonly Dictionary<int, string> contigNames = [];

        /// <summary>
        /// 自動推定された(あるいは CLI で明示指定された)インサートサイズ。
        /// CLI 指定がある場合はそれを、ない場合は ContigMaker.InsertSizeSamples から
        /// 推定した値を保持する。推定に失敗した場合は null のままとなり、
        /// その場合 Run はスキャフォールディングを行わずに終了する。
        /// </summary>
        public int? EffectiveInsertSize { get; private set; }

        /// <summary>
        /// スキャフォールディングを実行し、scaffoldPath に結果を書き出す。
        /// InsertSize が(指定・推定いずれの方法でも)確定できなかった場合は、
        /// その旨をログに出力して何もせずに戻る(scaffoldPath は作成されない)。
        /// </summary>
        public void Run(string scaffoldPath)
        {
            if (!this.TryResolveInsertSize(out var insertSize))
            {
                Console.WriteLine("[Info] Scaffolding skipped: insert size was not specified and could not be estimated from mapped pairs.");
                return;
            }
            this.EffectiveInsertSize = insertSize;
            Console.WriteLine($"[Info] Scaffolding with insert size = {insertSize}");

            this.LoadContigs();

            if (this.contigSequences.Count == 0)
            {
                Console.WriteLine("[Info] Scaffolding skipped: no contigs were found.");
                return;
            }

            var placements = contigMaker.UnitigPlacements;
            var pairPath = contigMaker.PairPath;

            // contig 単位の頂点空間を作る。unitig 同様、各 contig を
            // 「順方向」「逆方向」の2頂点として扱う。
            // vertexIndex = contigId << 1 (順方向) / contigId << 1 | 1 (逆方向)
            var contigCount = this.contigSequences.Keys.Count == 0 ? 0 : this.contigSequences.Keys.Max();
            var vertexCount = (contigCount + 1) << 1;

            // contig 単位の隣接候補: vertexIndex -> List<(vertexIndex, count, List<gapEstimates>)>
            var adjacency = new List<(int To, ulong Count, List<int> GapSamples)>[vertexCount];
            for (var i = 0; i < vertexCount; i++)
            {
                adjacency[i] = [];
            }

            // (fromVertex,toVertex) -> 集計済みのカウントとギャップ推定サンプル。
            var edgeMap = new Dictionary<(int, int), (ulong Count, List<int> GapSamples)>();

            var skippedInteriorEdges = 0;
            var skippedUnplacedEdges = 0;

            foreach (var (pathKey, gapSamples) in pairPath)
            {
                var (fromUnitig, toUnitig) = pathKey;

                if (!TryResolveContigEndpoint(placements, fromUnitig, isOutgoing: true, out var fromVertex))
                {
                    if (!placements.ContainsKey(Math.Abs(fromUnitig)))
                    {
                        skippedUnplacedEdges++;
                    }
                    else
                    {
                        skippedInteriorEdges++;
                    }
                    continue;
                }

                if (!TryResolveContigEndpoint(placements, toUnitig, isOutgoing: false, out var toVertex))
                {
                    if (!placements.ContainsKey(Math.Abs(toUnitig)))
                    {
                        skippedUnplacedEdges++;
                    }
                    else
                    {
                        skippedInteriorEdges++;
                    }
                    continue;
                }

                // 自己ループ(同一 contig の同一末端同士)は無視する。
                if (fromVertex >> 1 == toVertex >> 1)
                {
                    continue;
                }

                var edgeKey = (fromVertex, toVertex);
                if (edgeMap.TryGetValue(edgeKey, out var existing))
                {
                    existing.GapSamples.AddRange(gapSamples);
                    edgeMap[edgeKey] = (existing.Count + (ulong)gapSamples.Count, existing.GapSamples);
                }
                else
                {
                    edgeMap[edgeKey] = ((ulong)gapSamples.Count, [.. gapSamples]);
                }
            }

            if (skippedInteriorEdges > 0)
            {
                Console.WriteLine($"[Info] {skippedInteriorEdges} pair-end candidate(s) pointed at unitigs interior to an already-joined contig and were skipped (endpoint already resolved by UniteContigs).");
            }
            if (skippedUnplacedEdges > 0)
            {
                Console.WriteLine($"[Info] {skippedUnplacedEdges} pair-end candidate(s) referenced unitigs that were not placed into any contig (e.g. too short) and were skipped.");
            }

            // 辺 v→w と、その逆鎖側の双子 w^1→v^1 は同一の物理的な隣接を表す。
            // ペアエンドの観測はどちらか一方の向きでしか記録されないため、
            // 対称化しないと「順鎖側では十分な支持があるのに逆鎖側では
            // 支持ゼロ」という状態になり、下の相互一意性の検査が常に落ちる。
            // 双方に同じ支持(サンプルの和集合)を持たせる。
            // 二重計上にはならない: ある観測は edgeMap 上のどちらか一方の
            // キーにしか入っていないため、和を取ると各観測はちょうど1回ずつ数えられる。
            Dictionary<(int, int), (ulong Count, List<int> GapSamples)> symmetric = [];
            foreach (var ((from, to), (count, gapSamples)) in edgeMap)
            {
                foreach (var key in new[] { (from, to), (to ^ 1, from ^ 1) })
                {
                    if (symmetric.TryGetValue(key, out var acc))
                    {
                        acc.GapSamples.AddRange(gapSamples);
                        symmetric[key] = (acc.Count + count, acc.GapSamples);
                    }
                    else
                    {
                        symmetric[key] = (count, [.. gapSamples]);
                    }
                }
            }

            foreach (var ((from, to), (count, gapSamples)) in symmetric)
            {
                adjacency[from].Add((to, count, gapSamples));
            }

            var uniteThreshold = ConfigurationManager.Arguments.PairUniteThreshold;
            var countThreshold = ConfigurationManager.Arguments.PairCountThreshold;

            Console.WriteLine($"[Info] Scaffold candidate edges (contig-level, before thresholding): {edgeMap.Count}");

            // 各頂点について、最多支持のエッジ1本だけを残す(FixPath と同じロジック)。
            var resolvedEdge = new (int To, int GapLength)?[vertexCount];
            for (var v = 2; v < vertexCount; v++)
            {
                this.FixScaffoldEdge(adjacency, v, uniteThreshold, countThreshold, resolvedEdge);
            }

            var resolvedCount = 0;
            for (var v = 2; v < vertexCount; v++)
            {
                if (resolvedEdge[v] != null)
                {
                    resolvedCount++;
                }
            }

            // 相互一意(mutual unique)な辺だけを採用する。v→w を繋いでよいのは
            // 「v の唯一の行き先が w」であり、かつ「w の唯一の来訪元が v」で
            // あるときに限る。後者は逆鎖対称性より resolvedEdge[w^1] が v^1 を
            // 指すことと同値。これを課さないと、複数の contig が同じ次の contig を
            // 指した場合に先着1本だけが繋がれ、残りは黙って千切れる
            // (どれが正しいかの根拠がないまま1本を選ぶことになる)。
            var candidateEdge = (( int To, int GapLength)?[])resolvedEdge.Clone();
            var rejectedByReciprocity = 0;
            for (var v = 2; v < vertexCount; v++)
            {
                if (candidateEdge[v] is not { } edge)
                {
                    continue;
                }
                var twin = edge.To ^ 1;
                if (twin >= vertexCount || candidateEdge[twin] is not { } back || back.To != (v ^ 1))
                {
                    resolvedEdge[v] = null;
                    rejectedByReciprocity++;
                }
            }
            Console.WriteLine($"[Info] Scaffold edges resolved after thresholding: {resolvedCount}; {rejectedByReciprocity} rejected by the mutual-uniqueness check, {resolvedCount - rejectedByReciprocity} kept.");

            // 「入ってくる結合を持たない」頂点が経路の始点。v への結合が
            // 存在することは、逆鎖対称性より resolvedEdge[v^1] != null と同値。
            var startVertices = new List<int>();
            for (var v = 2; v < vertexCount; v++)
            {
                if (this.contigSequences.ContainsKey(v >> 1) && (v ^ 1) < vertexCount && resolvedEdge[v ^ 1] == null)
                {
                    startVertices.Add(v);
                }
            }

            List<string> scaffoldList = [];
            var visited = new bool[vertexCount];
            foreach (var start in startVertices)
            {
                // startVertices には同一 contig の fwd/rev 両方の頂点が
                // 独立に含まれうる(両方とも enterCount==0 の場合)。
                // 先に処理された方の WalkScaffold が MarkContigVisited で
                // 両方向を visited にするため、後から来た方はここで
                // スキップしないと、同じ contig を起点とする scaffold が
                // 二重に生成されてしまう(contig 数の水増し・配列の重複の原因)。
                if (visited[start])
                {
                    continue;
                }
                var scaffold = this.WalkScaffold(resolvedEdge, start, visited);
                if (scaffold != null)
                {
                    scaffoldList.Add(scaffold);
                }
            }

            // まだ訪問されていない(=孤立した、あるいは循環に巻き込まれた)contig を
            // 単独スキャフォールドとして出力する。
            for (var contigId = 1; contigId <= contigCount; contigId++)
            {
                var fwd = contigId << 1;
                var rev = (contigId << 1) | 1;
                if (fwd < vertexCount && !visited[fwd] && !visited[rev] && this.contigSequences.TryGetValue(contigId, out var value))
                {
                    scaffoldList.Add(value);
                    visited[fwd] = true;
                    visited[rev] = true;
                }
            }

            using var writer = new FastaWriter(scaffoldPath);
            var scaffoldId = 1;
            long totalLength = 0;
            foreach (var scaffold in scaffoldList)
            {
                writer.Write($"SCAFFOLD{scaffoldId}", scaffold);
                scaffoldId++;
                totalLength += scaffold.Length;
            }

            Console.WriteLine($"[Info] Wrote {scaffoldList.Count} scaffold(s), total length {totalLength}, to {scaffoldPath}");
        }

        /// <summary>
        /// 同一unitig内サンプルを信頼してよいと判断する、
        /// 「unitig長 / 推定フラグメント長」の下限比。
        ///
        /// 同一unitig内サンプルの唯一の弱点は、unitig がフラグメントより短いと
        /// 両端が収まるペアしか観測できず短いフラグメントに偏ること。
        /// unitig がフラグメント長よりこの倍率以上に長ければ、その打ち切りは
        /// 事実上起きないため偏りは無視できる。
        /// </summary>
        private const int UnbiasedSameUnitigLengthRatio = 10;

        /// <summary>
        /// InsertSize を確定する。CLI で明示指定されていればそれを使う。
        ///
        /// 未指定の場合、2種類のサンプル群から選ぶ。
        ///
        /// 同一unitig内サンプルは、unitig自体がフラグメント長より短いと
        /// 両端が収まるペアしか観測できず、より短いフラグメントに偏った標本に
        /// なる。ただしこの打ち切りは unitig が十分長ければ起きないため、
        /// unitig の N50 が推定値の UnbiasedSameUnitigLengthRatio 倍以上ある
        /// 場合は、桁違いに多い標本数(実データで76万件 対 600件)を活かして
        /// こちらを採用する。
        ///
        /// resolved-edge由来は unitig 長に制約されない代わりに、
        /// 「結合が確定した辺」だけを対象とするため標本数が非常に少なく、
        /// 誤結合や誤マッピングの影響を受けやすい(実データで同一unitig側が
        /// 中央値245・分布も素直だったのに対し、resolved-edge側は437と
        /// 8割ほど上振れしていた)。
        /// </summary>
        private bool TryResolveInsertSize(out int insertSize)
        {
            if (ConfigurationManager.Arguments.InsertSize is { } specified)
            {
                insertSize = specified;
                return true;
            }

            var sameUnitigSamples = contigMaker.SameUnitigInsertSizeSamples;
            if (sameUnitigSamples.Count >= Consts.MinInsertSizeSampleCount)
            {
                var sameUnitigEstimate = Median(sameUnitigSamples);
                var unitigN50 = UnitigN50(contigMaker.UnitigLengths);
                if (sameUnitigEstimate > 0 && unitigN50 >= (long)sameUnitigEstimate * UnbiasedSameUnitigLengthRatio)
                {
                    insertSize = sameUnitigEstimate;
                    Console.WriteLine(
                        $"[Info] Insert size auto-estimated as {insertSize} from {sameUnitigSamples.Count} same-unitig sampled pairs " +
                        $"(median; unitig N50 {unitigN50} is >= {UnbiasedSameUnitigLengthRatio}x the estimate, so the short-fragment truncation bias does not apply).");
                    return true;
                }
            }

            var resolvedEdgeSamples = contigMaker.ResolvedEdgeInsertSizeSamples;
            if (resolvedEdgeSamples.Count >= Consts.MinInsertSizeSampleCount)
            {
                insertSize = Median(resolvedEdgeSamples);
                Console.WriteLine($"[Info] Insert size auto-estimated as {insertSize} from {resolvedEdgeSamples.Count} resolved-edge sampled pairs (median, preferred over same-unitig samples because the unitigs are not long enough for same-unitig samples to be unbiased).");
                return true;
            }

            var allSamples = contigMaker.InsertSizeSamples;
            if (allSamples.Count < Consts.MinInsertSizeSampleCount)
            {
                Console.WriteLine($"[Info] Insert size auto-estimation requires at least {Consts.MinInsertSizeSampleCount} samples; only {resolvedEdgeSamples.Count} resolved-edge and {allSamples.Count} total samples were collected.");
                insertSize = 0;
                return false;
            }

            insertSize = Median(allSamples);
            Console.WriteLine($"[Info] Insert size auto-estimated as {insertSize} from {allSamples.Count} sampled pairs (median; resolved-edge samples were too few ({resolvedEdgeSamples.Count}), fell back to the full pool which may be biased short).");
            return true;
        }

        /// <summary>
        /// unitig の N50(長い順に並べて累積長が全長の半分に達した時点の長さ)。
        /// 「同一unitig内サンプルが打ち切りバイアスを受けていないか」の判断に使う。
        /// 平均ではなく N50 を使うのは、短い断片が本数として多くても
        /// 実際にペアが観測される場所はゲノムの大部分を占める長い unitig に
        /// 偏るため、そちらの長さ水準を見るべきだから。
        /// </summary>
        private static long UnitigN50(IReadOnlyDictionary<int, int> unitigLengths)
        {
            if (unitigLengths.Count == 0)
            {
                return 0;
            }
            var lengths = unitigLengths.Values.OrderByDescending(x => x).ToList();
            var half = lengths.Sum(x => (long)x) / 2.0;
            long cumulative = 0;
            foreach (var length in lengths)
            {
                cumulative += length;
                if (cumulative >= half)
                {
                    return length;
                }
            }
            return lengths[^1];
        }

        private static int Median(List<int> values)
        {
            var sorted = values.OrderBy(x => x).ToList();
            var mid = sorted.Count / 2;
            return sorted.Count % 2 == 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
        }

        private void LoadContigs()
        {
            using var reader = new FastaReader(contigFilePath);
            var id = 1;
            while (reader.HasNext())
            {
                var seq = reader.NextSequence();
                this.contigNames[id] = seq.ID.TrimStart('>');
                this.contigSequences[id] = seq.Seq;
                id++;
            }
        }

        /// <summary>
        /// pairPath のキーに現れる unitig ID(符号付き)を、それが実際に
        /// contig の末端(スキャフォールディング候補として使える位置)に
        /// 配置されているかどうか判定し、配置されていれば対応する
        /// contig 頂点(vertexIndex = contigId &lt;&lt; 1 | reverseFlag)を返す。
        ///
        /// isOutgoing = true の場合(from 側、つまり読み進める方向の始点):
        ///   unitig がその向きで見て contig の「出口側」の末端、すなわち
        ///   - 向きが順鎖(id&gt;0)かつ contig 内で末尾(IsAtContigEnd)、または
        ///   - 向きが逆鎖(id&lt;0)かつ contig 内で先頭(IsAtContigStart)
        ///   である場合のみ有効な候補とみなす。
        /// isOutgoing = false の場合(to 側、読み進めた先の入口):
        ///   - 向きが順鎖(id&gt;0)かつ contig 内で先頭(IsAtContigStart)、または
        ///   - 向きが逆鎖(id&lt;0)かつ contig 内で末尾(IsAtContigEnd)
        ///   である場合のみ有効。
        ///
        /// contig の最終配列が正規化のため逆相補化されている場合
        /// (IsContigReverseComplemented)、walk 順ベースの先頭/末尾の意味が
        /// 反転するため、その分も考慮して vertex の reverseFlag を決める。
        /// </summary>
        private static bool TryResolveContigEndpoint(
            IReadOnlyDictionary<int, UnitigPlacement> placements,
            int signedUnitigId,
            bool isOutgoing,
            out int vertexIndex)
        {
            vertexIndex = 0;
            var unitigId = Math.Abs(signedUnitigId);
            var isForwardUnitig = signedUnitigId > 0;

            if (!placements.TryGetValue(unitigId, out var placement))
            {
                return false;
            }

            // unitig 自身が walk 中に逆鎖として使われていた場合、pairPath 上の
            // 向き(isForwardUnitig)は「unitig 単体の元の向き」を基準にしているため、
            // walk 内での実効的な向きに変換する(XOR)。
            var effectiveForward = isForwardUnitig != placement.IsUnitigReverseInWalk;

            var isAtRelevantEnd = isOutgoing
                ? effectiveForward ? placement.IsAtContigEnd : placement.IsAtContigStart
                : effectiveForward ? placement.IsAtContigStart : placement.IsAtContigEnd;
            if (!isAtRelevantEnd)
            {
                return false;
            }

            // contig 全体が正規化のために逆相補化されている場合、
            // 「walk 順で見た先頭/末尾」と「実際の contigs.fasta 上の先頭/末尾」が
            // 入れ替わる。スキャフォールディングは contigs.fasta 上の配列
            // (=実際に出力された向き)を基準に扱うため、ここで反転させる。
            var isForwardInFinalSequence = placement.IsContigReverseComplemented ? !effectiveForward : effectiveForward;

            vertexIndex = (placement.ContigId << 1) | (isForwardInFinalSequence ? 0 : 1);
            return true;
        }

        private void FixScaffoldEdge(
            List<(int To, ulong Count, List<int> GapSamples)>[] adjacency,
            int vertex,
            decimal uniteThreshold,
            ulong countThreshold,
            (int To, int GapLength)?[] resolvedEdge)
        {
            var candidates = adjacency[vertex];
            var filtered = candidates.Where(c => c.Count >= countThreshold).ToList();
            if (filtered.Count == 0)
            {
                resolvedEdge[vertex] = null;
                return;
            }

            var sum = filtered.Aggregate(0UL, (acc, c) => acc + c.Count);
            var (To, Count, GapSamples) = filtered.OrderByDescending(c => c.Count).First();

            if (sum == 0 || (decimal)Count / sum < uniteThreshold)
            {
                resolvedEdge[vertex] = null;
                return;
            }

            var gapLength = this.EstimateGapLength(GapSamples);
            resolvedEdge[vertex] = (To, gapLength);
        }

        /// <summary>
        /// あるエッジについて観測された「既に見えている長さ」のサンプル群から
        /// 実際に挿入する N の数を決める。ContigMaker が記録する各サンプルは
        /// 「read1長 + contig1末端までの残り + contig2先頭からの残り + read2長」
        /// であり、フラグメント長 = サンプル + ギャップ長 という関係が成り立つ。
        /// したがってギャップ長 = InsertSize - サンプル として個々の候補を計算し、
        /// その中央値を採用する。
        /// 中央値が Consts.MinimumGapLength を下回る場合はその最小値に丸める
        /// (負の推定値や 0 になった場合でも、隣接している事実自体は
        /// 相応の証拠があるため、少なくとも1つの N でギャップを明示する)。
        /// </summary>
        private int EstimateGapLength(List<int> spannedLengthSamples)
        {
            var insertSize = this.EffectiveInsertSize ?? 0;
            if (spannedLengthSamples.Count == 0)
            {
                return Consts.MinimumGapLength;
            }

            var gapEstimates = spannedLengthSamples
                .Select(spanned => insertSize - spanned)
                .OrderBy(x => x)
                .ToList();

            var mid = gapEstimates.Count / 2;
            var median = gapEstimates.Count % 2 == 0
                ? (gapEstimates[mid - 1] + gapEstimates[mid]) / 2
                : gapEstimates[mid];

            return Math.Max(Consts.MinimumGapLength, median);
        }

        private string? WalkScaffold((int To, int GapLength)?[] resolvedEdge, int start, bool[] visited)
        {
            var contigId = start >> 1;
            var isReverse = (start & 1) == 1;
            if (!this.contigSequences.TryGetValue(contigId, out var seq))
            {
                return null;
            }

            var sb = new StringBuilder(isReverse ? Util.ReverseComprement(seq) : seq);
            var current = start;
            // 頂点を「消費」した(=いずれかの向きでスキャフォールドに組み込んだ)際は、
            // その contig の両方の向きの頂点(fwd/rev)を visited にする。
            // 片方の頂点だけを visited にすると、同じ contig の反対向きの頂点が
            // 別の開始点や「未訪問の孤立 contig」判定で再度使われてしまう
            // (同じ contig が2回出力される)おそれがあるため。
            MarkContigVisited(visited, current);
            while (resolvedEdge[current] is { } edge && !visited[edge.To])
            {
                var nextContigId = edge.To >> 1;
                var nextIsReverse = (edge.To & 1) == 1;
                if (!this.contigSequences.TryGetValue(nextContigId, out var nextSeq))
                {
                    break;
                }

                _ = sb.Append('N', edge.GapLength);
                _ = sb.Append(nextIsReverse ? Util.ReverseComprement(nextSeq) : nextSeq);

                current = edge.To;
                MarkContigVisited(visited, current);
            }

            return sb.ToString();
        }

        /// <summary>
        /// vertexIndex が指す contig の両方の向きの頂点を visited にする。
        /// </summary>
        private static void MarkContigVisited(bool[] visited, int vertexIndex)
        {
            var contigId = vertexIndex >> 1;
            var fwd = contigId << 1;
            var rev = fwd | 1;
            if (rev < visited.Length)
            {
                visited[fwd] = true;
                visited[rev] = true;
            }
        }
    }
}