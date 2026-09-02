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

            foreach (var ((from, to), (count, gapSamples)) in edgeMap)
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

            // enterCount: 各頂点を「唯一の行き先」として指している頂点の数。
            // 2以上の場合、競合により WalkScaffold の visited 管理で
            // 最初に到達した経路のみが実際に結合される。
            var enterCount = new int[vertexCount];
            for (var v = 2; v < vertexCount; v++)
            {
                if (resolvedEdge[v] is { } edge)
                {
                    enterCount[edge.To]++;
                }
            }

            var resolvedCount = 0;
            for (var v = 2; v < vertexCount; v++)
            {
                if (resolvedEdge[v] != null)
                {
                    resolvedCount++;
                }
            }
            Console.WriteLine($"[Info] Scaffold edges resolved after thresholding: {resolvedCount}");

            var startVertices = new List<int>();
            for (var v = 2; v < vertexCount; v++)
            {
                if (this.contigSequences.ContainsKey(v >> 1) && enterCount[v] == 0)
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
        /// InsertSize を確定する。CLI で明示指定されていればそれを使う。
        /// 未指定の場合、ContigMaker.InsertSizeSamples から中央値を推定する。
        /// サンプル数が Consts.MinInsertSizeSampleCount 未満の場合は推定を諦め、
        /// false を返す(呼び出し側はスキャフォールディングをスキップする)。
        /// </summary>
        private bool TryResolveInsertSize(out int insertSize)
        {
            if (ConfigurationManager.Arguments.InsertSize is { } specified)
            {
                insertSize = specified;
                return true;
            }

            var samples = contigMaker.InsertSizeSamples;
            if (samples.Count < Consts.MinInsertSizeSampleCount)
            {
                Console.WriteLine($"[Info] Insert size auto-estimation requires at least {Consts.MinInsertSizeSampleCount} samples (both reads mapping uniquely to the same unitig); only {samples.Count} were collected.");
                insertSize = 0;
                return false;
            }

            var sorted = samples.OrderBy(x => x).ToList();
            var mid = sorted.Count / 2;
            var median = sorted.Count % 2 == 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
            Console.WriteLine($"[Info] Insert size auto-estimated as {median} from {samples.Count} sampled pairs (median).");
            insertSize = median;
            return true;
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
        /// あるエッジについて観測された「未読了長の合計」サンプル群から
        /// 実際に挿入する N の数を決める。ギャップ長 = InsertSize - 未読了長合計、
        /// という関係で個々のサンプルからギャップ長候補を計算し、その中央値を採用する。
        /// 中央値が Consts.MinimumGapLength を下回る場合はその最小値に丸める
        /// (負の推定値や 0 になった場合でも、隣接している事実自体は
        /// 相応の証拠があるため、少なくとも1つの N でギャップを明示する)。
        /// </summary>
        private int EstimateGapLength(List<int> totalRemainingLengthSamples)
        {
            var insertSize = this.EffectiveInsertSize ?? 0;
            if (totalRemainingLengthSamples.Count == 0)
            {
                return Consts.MinimumGapLength;
            }

            var gapEstimates = totalRemainingLengthSamples
                .Select(remaining => insertSize - remaining)
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