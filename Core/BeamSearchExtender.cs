using Tsumiki.Common;

namespace Tsumiki.Core
{
    /// <summary>
    /// 貪欲な相互一意性の判定では決めきれなかった分岐を、先読み(ビームサーチ)で
    /// 解けるだけ解く。
    ///
    /// 考え方:
    /// 相互一意性の検査は「その1歩だけ」を見て判断するため、分岐の直後だけを見ると
    /// 五分五分に見えるが、2〜3本先まで進めると片方だけがペアエンドの証拠と
    /// 整合する、という状況を取りこぼす。そこで各候補から数kb先まで複数の経路を
    /// 並行して伸ばし(ビーム)、その間に得られるペアエンドの支持を積算して比べる。
    ///
    /// 証拠は「いま組み上がっている contig の末尾インサートサイズぶんに載った
    /// リードの相方が、候補側のどの unitig に載ったか」で数える。片方の枝にだけ
    /// 相方が集中していれば、そちらが本当の続きである。
    ///
    /// 安全側の設計:
    /// - 探索の結果、上位の経路群が最初の1歩から割れている場合は何もしない。
    ///   ビームサーチの利点は「広く探して、全ての有力な仮説が一致する部分にだけ
    ///   コミットする」ことにあり、僅差で1本を選ぶことではない。
    /// - どの候補にもペアエンドの支持が無ければ何もしない(根拠が無い)。
    /// - コピー数を予算として持ち、反復配列を予算以上に通らない。予算が無いと
    ///   探索は同じ反復を何度でも通れてしまい、ありもしない長い経路を作る。
    ///
    /// これは「ゲノム全体を1本のオイラー路として探す」問題の、局所的で保守的な
    /// 近似にあたる。全体を一度に探索すると、正解以外の経路も同じだけ存在するため
    /// 誤アセンブリを大量に生む。実際に効くのは「あと1〜数歩だけ確実に伸ばす」ところ。
    /// </summary>
    internal static class BeamSearchExtender
    {
        /// <summary>
        /// 先読みで進む塩基数の上限。長くするほど遠くの証拠を使えるが、
        /// 探索が広がるうえ、遠いほどペアエンドの証拠は届かなくなる。
        /// インサートサイズの数倍あれば、跨げる範囲は使い切れる。
        /// </summary>
        private const int LookaheadMultiplier = 4;

        /// <summary>1つの分岐あたりに保持する部分経路の数。</summary>
        private const int DefaultBeamWidth = 8;

        /// <summary>先読みの1経路あたりの最大ステップ数(暴走防止)。</summary>
        private const int MaxStepsPerPath = 40;

        private sealed class State
        {
            public required int Current { get; init; }
            public required int FirstStep { get; init; }
            public required long Score { get; init; }
            public required int Length { get; init; }
            public required Dictionary<int, int> Used { get; init; }
        }

        /// <summary>
        /// merge が未確定(-1)の頂点について、先読みで続きを決められるものを決める。
        /// merge を直接書き換える。戻り値は新たに確定した結合の数(有向、双子ぶんを含む)。
        /// </summary>
        public static int Extend(
            UnitigGraph graph,
            List<string> unitigList,
            int[] merge,
            IReadOnlyDictionary<(int, int), ulong> pairLink,
            IReadOnlyDictionary<int, int> copyNumber,
            int insertSize,
            decimal dominanceThreshold,
            ulong minimumEvidence)
        {
            var lookaheadBases = Math.Max(insertSize, 1) * LookaheadMultiplier;
            var committed = 0;

            for (var v = 2; v < graph.VertexCount; v++)
            {
                if (merge[v] != -1)
                {
                    continue;
                }
                var candidates = graph.OutEdges[v];
                if (candidates.Count == 0)
                {
                    continue;
                }

                var anchors = CollectAnchors(v, unitigList, merge, insertSize, copyNumber);
                if (anchors.Count == 0)
                {
                    // 単一コピーの足場が1つも取れない = いま反復配列の上にいて、
                    // どのコピーにいるのか分からない。この状態で進む方向を選ぶ
                    // 根拠は原理的に存在しない。
                    continue;
                }

                var best = SearchBestFirstStep(
                    graph, unitigList, v, anchors, pairLink, copyNumber,
                    lookaheadBases, dominanceThreshold, minimumEvidence);
                if (best is not { } chosen)
                {
                    continue;
                }

                // 相互一意性は保ったままにする。行き先に既に別の結合が
                // 入っている場合は、そちらを壊してまで繋がない。
                if (merge[chosen ^ 1] != -1 || merge[chosen] == (v ^ 1))
                {
                    continue;
                }

                // 解きほぐされていない多コピーの反復を通り抜ける結合は作らない。
                // A-R-B-R-C という構造で A→R と R→C はどちらも本物の隣接だが、
                // R を1回しか使えない walk でこれを連鎖させると中間の B が
                // 飛ばされる(詳細は ContigMaker 側の同名の判定を参照)。
                if (!CanChainThrough(graph, copyNumber, v) || !CanChainThrough(graph, copyNumber, chosen ^ 1))
                {
                    continue;
                }

                merge[v] = chosen;
                merge[chosen ^ 1] = v ^ 1;
                committed += 2;
            }

            return committed;
        }

        /// <summary>
        /// v から遡って、いま組み上がっている contig の末尾インサートサイズぶんに
        /// あたる頂点のうち、<b>単一コピーのものだけ</b>を集める。
        /// ここに載ったリードの相方が、続きの証拠になる。
        /// 直前の頂点は、逆鎖対称性より merge[v^1] の双子で辿れる。
        ///
        /// 多コピーの unitig を足場から外すのが要点。反復配列の内部から読まれた
        /// リードは、どのコピー由来か区別できない(それが反復が解けない理由そのもの)。
        /// そこを起点にしたペアの証拠はどの行き先にも付いてしまい、標本数が少ないと
        /// 偶然の偏りが閾値を超えて誤った側を選ぶ。
        ///
        /// 実際、反復入りの合成ゲノム(A-R-B-R-C、R は150bpの2コピー反復)で、
        /// R 自身を足場にしたために A-R-C という中間を飛ばした contig が
        /// 出力されていた(真値照合で発覚)。単一コピーに限れば、A を足場として
        /// A-B のペアだけが支持されるため正しく B が選ばれる。
        ///
        /// 通過はするが足場には数えない、という扱いにする(多コピー領域の
        /// 向こう側にある単一コピー領域は、証拠として有効なため)。
        /// </summary>
        private static List<int> CollectAnchors(
            int v, List<string> unitigList, int[] merge, int insertSize, IReadOnlyDictionary<int, int> copyNumber)
        {
            List<int> anchors = [];
            List<int> walked = [v];
            if (copyNumber.GetValueOrDefault(v >> 1, 1) <= 1)
            {
                anchors.Add(v);
            }

            var accumulated = unitigList[v].Length;
            var current = v;
            while (accumulated < insertSize)
            {
                var twinMerge = merge[current ^ 1];
                if (twinMerge == -1)
                {
                    break;
                }
                var predecessor = twinMerge ^ 1;
                if (predecessor == current || walked.Contains(predecessor))
                {
                    break;
                }
                walked.Add(predecessor);
                if (copyNumber.GetValueOrDefault(predecessor >> 1, 1) <= 1)
                {
                    anchors.Add(predecessor);
                }
                accumulated += unitigList[predecessor].Length;
                current = predecessor;
            }
            return anchors;
        }

        /// <summary>
        /// v からの各候補について先読みし、最初の1歩として最も支持される頂点を返す。
        /// 決めきれない場合は null。
        /// </summary>
        private static int? SearchBestFirstStep(
            UnitigGraph graph,
            List<string> unitigList,
            int v,
            List<int> anchors,
            IReadOnlyDictionary<(int, int), ulong> pairLink,
            IReadOnlyDictionary<int, int> copyNumber,
            int lookaheadBases,
            decimal dominanceThreshold,
            ulong minimumEvidence)
        {
            List<State> beam = [];
            foreach (var w in graph.OutEdges[v])
            {
                var budget = copyNumber.GetValueOrDefault(w >> 1, 1);
                if (budget <= 0)
                {
                    continue;
                }
                beam.Add(new State
                {
                    Current = w,
                    FirstStep = w,
                    Score = ScoreOf(anchors, w, pairLink),
                    Length = unitigList[w].Length,
                    Used = new Dictionary<int, int> { [w >> 1] = 1 },
                });
            }
            if (beam.Count == 0)
            {
                return null;
            }

            // 最初の1歩ごとの最良スコアを追跡する。
            Dictionary<int, long> bestByFirstStep = [];
            foreach (var state in beam)
            {
                bestByFirstStep[state.FirstStep] = Math.Max(bestByFirstStep.GetValueOrDefault(state.FirstStep), state.Score);
            }

            for (var step = 0; step < MaxStepsPerPath && beam.Count > 0; step++)
            {
                List<State> next = [];
                foreach (var state in beam)
                {
                    if (state.Length >= lookaheadBases)
                    {
                        continue;
                    }
                    foreach (var w in graph.OutEdges[state.Current])
                    {
                        var unitigId = w >> 1;
                        var budget = copyNumber.GetValueOrDefault(unitigId, 1);
                        if (state.Used.GetValueOrDefault(unitigId) >= budget)
                        {
                            // 予算切れ。反復を何度も通って架空の経路を作らないようにする。
                            continue;
                        }
                        var used = new Dictionary<int, int>(state.Used);
                        used[unitigId] = used.GetValueOrDefault(unitigId) + 1;
                        next.Add(new State
                        {
                            Current = w,
                            FirstStep = state.FirstStep,
                            Score = state.Score + ScoreOf(anchors, w, pairLink),
                            Length = state.Length + unitigList[w].Length,
                            Used = used,
                        });
                    }
                }

                if (next.Count == 0)
                {
                    break;
                }

                next.Sort((x, y) => y.Score.CompareTo(x.Score));
                if (next.Count > DefaultBeamWidth)
                {
                    next = next[..DefaultBeamWidth];
                }
                beam = next;

                foreach (var state in beam)
                {
                    bestByFirstStep[state.FirstStep] = Math.Max(bestByFirstStep.GetValueOrDefault(state.FirstStep), state.Score);
                }
            }

            var ranked = bestByFirstStep.OrderByDescending(kv => kv.Value).ToList();
            var top = ranked[0];
            if ((ulong)Math.Max(0, top.Value) < minimumEvidence)
            {
                // どの枝にもペアエンドの支持が無い。根拠が無いので繋がない。
                return null;
            }

            var total = ranked.Sum(kv => Math.Max(0, kv.Value));
            if (total <= 0 || (decimal)top.Value / total < dominanceThreshold)
            {
                // 上位が割れている。僅差で選ぶくらいなら繋がないほうがよい。
                return null;
            }

            return top.Key;
        }

        /// <summary>
        /// その頂点を「通り抜けて」よいか。多コピーの反復は、解きほぐされて
        /// 入次数・出次数がどちらも1になっている場合(=どのコピーにいるかが
        /// 確定している場合)にだけ通り抜けてよい。
        /// </summary>
        private static bool CanChainThrough(UnitigGraph graph, IReadOnlyDictionary<int, int> copyNumber, int vertex)
        {
            if (copyNumber.GetValueOrDefault(vertex >> 1, 1) <= 1)
            {
                return true;
            }
            return graph.OutEdges[vertex].Count == 1 && graph.InDegree(vertex) == 1;
        }

        private static long ScoreOf(List<int> anchors, int candidate, IReadOnlyDictionary<(int, int), ulong> pairLink)
        {
            long score = 0;
            foreach (var anchor in anchors)
            {
                score += (long)pairLink.GetValueOrDefault((anchor, candidate));
            }
            return score;
        }
    }
}
