using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// unitig 間の隣接を de Bruijn グラフから厳密に構築する UnitigGraph の検証。
    /// 旧実装(リードマッピング由来の隣接候補 + 任意長オーバーラップ探索)は
    /// 実データで平均 2.96 塩基という偶然の一致で unitig を接着していたため、
    /// 「辺が張られる条件」そのものをここで固定する。
    /// </summary>
    public class UnitigGraphTests
    {
        private const int AmbiguousKmer = int.MinValue;

        /// <summary>
        /// ContigMaker のコンストラクタと同じ規則で kmerDict を組み立てる。
        /// 添字 2u が unitig u の順鎖、2u+1 が逆鎖。
        /// </summary>
        private static (List<string> UnitigList, Dictionary<KmerKey, (int UnitigId, int Position)> KmerDict) Build(
            int kmerLength,
            params string[] unitigs)
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = kmerLength, ThreadCount = 1 };
            List<string> unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> kmerDict = [];

            var id = 1;
            foreach (var seq in unitigs)
            {
                unitigList.Add(seq);
                unitigList.Add(Util.ReverseComprement(seq));

                for (var i = kmerLength; i <= seq.Length; i++)
                {
                    var startPos = i - kmerLength;
                    var key = new KmerKey(seq.AsSpan(startPos, kmerLength));
                    var revKey = key.ReverseComprement();
                    var revStartPos = seq.Length - i;
                    Register(kmerDict, key, id, startPos);
                    Register(kmerDict, revKey, -id, revStartPos);
                }
                id++;
            }
            return (unitigList, kmerDict);
        }

        private static void Register(Dictionary<KmerKey, (int, int)> dict, KmerKey key, int id, int position)
        {
            if (dict.TryGetValue(key, out var existing))
            {
                if (existing.Item1 is AmbiguousKmer || existing.Item1 == id)
                {
                    return;
                }
                dict[key] = (AmbiguousKmer, 0);
                return;
            }
            dict[key] = (id, position);
        }

        [Fact]
        public void Build_CreatesEdge_WhenOneUnitigTailExtendsIntoAnotherUnitigHead()
        {
            // unitig A の末尾 k-1 塩基が unitig B の先頭 k-1 塩基と一致する構成。
            // A の末尾 k-mer から 1 塩基伸ばすと、ちょうど B の先頭 k-mer になる。
            const int k = 8;
            const string shared = "CGTTACA"; // k-1 = 7 塩基の重なり
            var a = "GCTAAAGACAATTAC" + shared;      // 末尾が shared
            var b = shared + "GGATCCTTAGGCAAT";      // 先頭が shared

            var (unitigList, kmerDict) = Build(k, a, b);
            var graph = UnitigGraph.Build(unitigList, kmerDict, k, AmbiguousKmer);

            var aForward = ContigMaker.VertexIndex(1);
            var bForward = ContigMaker.VertexIndex(2);

            Assert.Contains(bForward, graph.OutEdges[aForward]);

            // 逆鎖対称性: A→B があるなら B' →A' も存在しなければならない。
            // これが崩れると順鎖側と逆鎖側で別々の経路が組まれ、同じ領域が
            // 2 通りに組み立てられてしまう。
            Assert.Contains(aForward ^ 1, graph.OutEdges[bForward ^ 1]);

            // 入次数は双子の出次数で表せる。
            Assert.Equal(graph.OutEdges[bForward ^ 1].Count, graph.InDegree(bForward));
        }

        [Fact]
        public void Build_CreatesNoEdge_WhenUnitigsDoNotOverlapByKMinusOne()
        {
            // 互いに無関係な2本。偶然の短い一致があっても辺は張られてはならない
            // (旧実装はここで平均3塩基程度の一致による誤結合を作っていた)。
            const int k = 8;
            const string a = "GCTAAAGACAATTACATAA";
            const string b = "TTGACCTGAATCCGGTTCA";

            var (unitigList, kmerDict) = Build(k, a, b);
            var graph = UnitigGraph.Build(unitigList, kmerDict, k, AmbiguousKmer);

            for (var v = 2; v < graph.VertexCount; v++)
            {
                Assert.Empty(graph.OutEdges[v]);
            }
        }

        [Fact]
        public void Build_CreatesNoEdge_WhenTargetKmerIsInteriorRatherThanHead()
        {
            // B の「途中」に一致する k-mer があっても、k-1 オーバーラップでの
            // 連結はできないため辺を張ってはならない(Position != 0 を弾く条件)。
            // A の末尾 k-mer を 1 塩基伸ばした k-mer が B の 3 塩基目から始まるよう構成する。
            const int k = 8;
            const string junction = "ACGGATCA"; // A の末尾から伸ばして得られる k-mer
            var a = "GCTAAAGACAATTAC" + junction[..(k - 1)]; // 末尾 k-1 が junction の先頭 k-1
            var b = "TT" + junction + "CCTTAGGCAAT";          // junction は B の位置 2 に現れる

            var (unitigList, kmerDict) = Build(k, a, b);
            var graph = UnitigGraph.Build(unitigList, kmerDict, k, AmbiguousKmer);

            var aForward = ContigMaker.VertexIndex(1);
            var bForward = ContigMaker.VertexIndex(2);
            Assert.DoesNotContain(bForward, graph.OutEdges[aForward]);
        }

        [Fact]
        public void Build_RecordsBothBranches_WhenTailExtendsIntoTwoDifferentUnitigs()
        {
            // 分岐: A の末尾から B へも C へも伸びられる構成。
            // 辺は両方張られ、どちらを選ぶかはリード支持に委ねられる。
            const int k = 8;
            const string shared = "CGTTACA";
            var a = "GCTAAAGACAATTAC" + shared;
            var b = shared + "GGATCCTTAGGCAAT";
            var c = shared + "TGATCCTTAGGCAAT"; // 分岐点の 1 塩基だけ B と異なる

            var (unitigList, kmerDict) = Build(k, a, b, c);
            var graph = UnitigGraph.Build(unitigList, kmerDict, k, AmbiguousKmer);

            var aForward = ContigMaker.VertexIndex(1);
            Assert.Equal(2, graph.OutEdges[aForward].Count);
            Assert.Contains(ContigMaker.VertexIndex(2), graph.OutEdges[aForward]);
            Assert.Contains(ContigMaker.VertexIndex(3), graph.OutEdges[aForward]);
        }

        [Fact]
        public void Build_CreatesNoSelfLoopEdge()
        {
            // 自己ループを辺として持つと walk が同じ unitig を無限に伸ばしうるため、
            // 構築段階で除外していることを確認する。
            const int k = 8;
            const string repeatUnit = "ACGGATCT";
            var a = repeatUnit + "GCTAAAGA" + repeatUnit[..(k - 1)];

            var (unitigList, kmerDict) = Build(k, a);
            var graph = UnitigGraph.Build(unitigList, kmerDict, k, AmbiguousKmer);

            var aForward = ContigMaker.VertexIndex(1);
            Assert.DoesNotContain(aForward, graph.OutEdges[aForward]);
        }
    }
}
