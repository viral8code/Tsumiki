using System.Buffers;
using System.Text;
using Tsumiki.Common;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    internal class BeamDBG(CountingBloomFilter bloomFilter, int B)
    {
        public class State
        {
            public int BranchCounts { get; set; }
            public required string End { get; set; }
            public int Len { get; set; }
            public int Revisits { get; set; }
            public ulong Hash { get; set; }
            public State? BeforeState { get; set; }
        }

        public Unitig ExtendFrom(string start)
        {
            PriorityQueue<State, double> beam = new();
            beam.Enqueue(new()
            {
                End = start,
                Len = 1,
                Revisits = 0,
                Hash = 0L,
                BeforeState = null,
            }, 0.0);

            State? bestSeen = null;
            var bestScore = double.MinValue;

            var cycle = 0;
            while (beam.Count > 0)
            {
                Console.WriteLine("cycle " + cycle++);
                List<(State state, double score)> cand = [];
                while (beam.TryDequeue(out var nowState, out var nowScore))
                {
                    var expanded = false;
                    for (var i = 0; i < 4; i++)
                    {
                        var nextNode = NodeExtend(nowState.End, i);
                        if (!bloomFilter.Contains(nextNode))
                        {
                            continue;
                        }
                        State nextState = new()
                        {
                            End = nextNode,
                            Len = nowState.Len + 1,
                            Hash = HashExtend(nowState.Hash, i),
                            BeforeState = nowState,
                        };
                        nextState.Revisits = (nextState.Hash == nowState.Hash) ? nowState.Revisits + 1 : 0;
                        for (var j = 0; j < 4; j++)
                        {
                            var subNextNode = NodeExtend(nextState.End, j);
                            if (bloomFilter.Contains(subNextNode))
                            {
                                nextState.BranchCounts += 1;
                            }
                        }
                        var sc = 0.0;
                        sc += ScoreEntropy(nextState);
                        sc -= LoopPenalty(nextState);
                        var nextScore = nowScore + sc;
                        if (Prunable(nextState, nextScore, bestSeen, bestScore))
                        {
                            continue;
                        }
                        cand.Add((nextState, nextScore));
                        expanded = true;
                    }
                    if (!expanded)
                    {
                        cand.Add((nowState, nowScore));
                    }
                }
                if (cand.Count == 0)
                {
                    break;
                }
                cand.Sort((s1, s2) => s2.score.CompareTo(s1.score));
                List<(State state, double score)> next = [];
                for (var i = 0; i < cand.Count && i < B; i++)
                {
                    next.Add(cand[i]);
                }
                (var head, var score) = next[0];
                if (bestSeen == null || score > bestScore)
                {
                    bestSeen = head;
                    bestScore = score;
                }
                if (Converged(next))
                {
                    break;
                }
                beam.Clear();
                foreach ((var state, var sore) in next)
                {
                    beam.Enqueue(state, score);
                }
            }
            return materialize(bestSeen!);
        }

        private static string NodeExtend(string node, int nextBase)
        {
            return node[1..] + Util.ByteToBase((byte)nextBase);
        }

        private static ulong HashExtend(ulong h, int nextBase)
        {
            h ^= (ulong)nextBase + 0x9e3779b97f4a7c15L;
            h ^= h << 13;
            h ^= (h >>> 7);
            h ^= h << 17;
            return h;
        }

        private static double ScoreEntropy(State t)
        {
            return t.BranchCounts switch
            {
                0 => 0.0,
                1 => 0.0,
                _ => Math.Log(1.0 / t.BranchCounts),
            };
        }

        private static double LoopPenalty(State t)
        {
            return t.Revisits > 0 ? 0.5 * t.Revisits : 0.0;
        }

        private static bool Prunable(State t, double nowScore, State? bestSeen, double bestScore)
        {
            if (bestSeen != null && nowScore < bestScore - 5.0)
            {
                return true;
            }
            return t.Revisits > 8;
        }

        private static bool Converged(List<(State state, double score)> beamTop)
        {
            if (beamTop.Count == 1)
            {
                return true;
            }
            var d = beamTop[0].score - beamTop[1].score;
            return d > 3.0 && beamTop[0].state.Len > 100;
        }

        private static Unitig materialize(State best)
        {
            StringBuilder stringBuilder = new();
            while (best.BeforeState != null)
            {
                _ = stringBuilder.Append(best.End[^1]);
                best = best.BeforeState;
            }
            foreach (var c in best.End.Reverse())
            {
                stringBuilder.Append(c);
            }
            return new Unitig(string.Join(string.Empty, stringBuilder.ToString().Reverse()));
        }

        public class Unitig(string sequence)
        {
            public override string ToString()
            {
                return sequence;
            }
        }
    }
}
