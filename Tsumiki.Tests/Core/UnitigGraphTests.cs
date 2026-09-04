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
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = kmerLength, A_スレッド数 = 1 };
            List<string> unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> kmerDict = [];

            var id = 1;
            foreach (var seq in unitigs)
            {
                unitigList.Add(seq);
                unitigList.Add(Util.V_逆相補(seq));

                for (var i = kmerLength; i <= seq.Length; i++)
                {
                    var startPos = i - kmerLength;
                    var key = new KmerKey(seq.AsSpan(startPos, kmerLength));
                    var revKey = key.Get_逆相補();
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
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, k, AmbiguousKmer);

            var aForward = ContigMaker.Get_頂点番号(1);
            var bForward = ContigMaker.Get_頂点番号(2);

            Assert.Contains(bForward, graph.A_出辺[aForward]);

            // 逆鎖対称性: A→B があるなら B' →A' も存在しなければならない。
            // これが崩れると順鎖側と逆鎖側で別々の経路が組まれ、同じ領域が
            // 2 通りに組み立てられてしまう。
            Assert.Contains(aForward ^ 1, graph.A_出辺[bForward ^ 1]);

            // 入次数は双子の出次数で表せる。
            Assert.Equal(graph.A_出辺[bForward ^ 1].Count, graph.Get_入次数(bForward));
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
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, k, AmbiguousKmer);

            for (var v = 2; v < graph.A_出辺.Count; v++)
            {
                Assert.Empty(graph.A_出辺[v]);
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
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, k, AmbiguousKmer);

            var aForward = ContigMaker.Get_頂点番号(1);
            var bForward = ContigMaker.Get_頂点番号(2);
            Assert.DoesNotContain(bForward, graph.A_出辺[aForward]);
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
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, k, AmbiguousKmer);

            var aForward = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, graph.A_出辺[aForward].Count);
            Assert.Contains(ContigMaker.Get_頂点番号(2), graph.A_出辺[aForward]);
            Assert.Contains(ContigMaker.Get_頂点番号(3), graph.A_出辺[aForward]);
        }

        /// <summary>
        /// 単純バブル(u から2本に分かれ、それぞれ1本の unitig を経て
        /// 同じ w へ再合流する)で、リード支持の高い枝だけが経路として
        /// 残ることを確認する。
        ///
        /// 結合の採用条件を相互一意にした結果、再合流点 w の入次数が
        /// 2 のままだと u から w へ至る経路が一切結合されなくなるため、
        /// この処理が無いとバブルのたびに contig が千切れる。
        /// </summary>
        [Fact]
        public void PopSimpleBubbles_KeepsTheBestSupportedBranch_AndRemovesTheOtherSymmetrically()
        {
            const int k = 8;
            // k=8 で全 unitig を通じて重複する正規化 k-mer が無いことを確認済みの構成。
            const string u = "GCTAAAGACAATTACGCA";
            const string b1 = "TTACGCAAGGATCCTGCACGT"; // u の末尾7塩基 + 'A' で始まり、w の先頭7塩基で終わる
            const string b2 = "TTACGCACTTAGCATGCACGT"; // 分岐点の1塩基だけ b1 と異なる同長の枝
            const string w = "TGCACGTAAGGCTTACCA";

            var (unitigList, kmerDict) = Build(k, u, b1, b2, w);
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, k, AmbiguousKmer);

            var uV = ContigMaker.Get_頂点番号(1);
            var b1V = ContigMaker.Get_頂点番号(2);
            var b2V = ContigMaker.Get_頂点番号(3);
            var wV = ContigMaker.Get_頂点番号(4);

            // 前提: バブル構造が実際に構築されている。
            Assert.Equal(2, graph.A_出辺[uV].Count);
            Assert.Equal(2, graph.Get_入次数(wV));

            // b1 側にだけリード支持を与える。
            Dictionary<(int, int), ulong> support = new()
            {
                [(uV, b1V)] = 40,
                [(uV, b2V)] = 3,
            };

            var popped = graph.V_除去_単純バブル(unitigList, support);

            Assert.Equal(1, popped);
            Assert.Equal([b1V], graph.A_出辺[uV]);
            Assert.Equal([wV], graph.A_出辺[b1V]);
            Assert.Empty(graph.A_出辺[b2V]);

            // 逆鎖側も対称に取り除かれていること(片側だけ消すと順鎖と逆鎖で
            // 別々の経路が組まれてしまう)。
            Assert.Equal(1, graph.Get_入次数(wV));
            Assert.DoesNotContain(b2V ^ 1, graph.A_出辺[wV ^ 1]);
            Assert.DoesNotContain(uV ^ 1, graph.A_出辺[b2V ^ 1]);
        }

        /// <summary>
        /// 長さが大きく異なる分岐は、同じ領域の別表現(バブル)ではなく
        /// 本物の分岐(反復配列の出入口など)である可能性が高いため、
        /// 支持の低い側であっても勝手に経路から外してはならない。
        /// </summary>
        [Fact]
        public void PopSimpleBubbles_LeavesBranchesOfVeryDifferentLengthsAlone()
        {
            const int k = 8;
            const string u = "GCTAAAGACAATTACGCA";
            const string b1 = "TTACGCAAGGATCCTGCACGT";
            // b2 は b1 より大幅に長い(長さ比が既定の閾値1.5を超える)。
            const string b2 = "TTACGCACTTAGCAGGTCCAATTGGACCAATGCACGT";
            const string w = "TGCACGTAAGGCTTACCA";

            var (unitigList, kmerDict) = Build(k, u, b1, b2, w);
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, k, AmbiguousKmer);

            var uV = ContigMaker.Get_頂点番号(1);
            Assert.Equal(2, graph.A_出辺[uV].Count);

            Dictionary<(int, int), ulong> support = new()
            {
                [(uV, ContigMaker.Get_頂点番号(2))] = 40,
                [(uV, ContigMaker.Get_頂点番号(3))] = 3,
            };

            var popped = graph.V_除去_単純バブル(unitigList, support);

            Assert.Equal(0, popped);
            Assert.Equal(2, graph.A_出辺[uV].Count);
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
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, k, AmbiguousKmer);

            var aForward = ContigMaker.Get_頂点番号(1);
            Assert.DoesNotContain(aForward, graph.A_出辺[aForward]);
        }
    }
}
