using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 短い反復配列の解きほぐし(repeat resolution)の検証。
    ///
    /// 反復配列 R がゲノム中に2回現れ、それぞれ A→R→C と B→R→D という文脈を
    /// 持つ場合、de Bruijn グラフ上では R は1個の頂点に潰れて入次数2・出次数2に
    /// なる。R の内部から読まれたリードはどちらのコピー由来か区別できないため、
    /// 分岐でのリード支持は原理的に5割前後にしかならず解けない。R を丸ごと
    /// 跨いだフラグメントだけが手がかりになる。
    ///
    /// 実データ(k=63)ではこの形の unitig が151本あり、うち143本が
    /// フラグメント長の中央値(245bp)より短かった。
    /// </summary>
    public class RepeatResolutionTests
    {
        private const int AmbiguousKmer = int.MinValue;
        private const int K = 8;

        // k=8 で5本すべてを通じて重複する正規化 k-mer が無いことを確認済みの構成。
        // A と B はどちらも R の先頭 k-1 塩基で終わり、C と D はどちらも
        // R の末尾 k-1 塩基で始まる(= R が入次数2・出次数2の反復になる)。
        private const string UnitigA = "ACAGTTCGCGAGCCCTCCGTC";
        private const string UnitigB = "TGTATTGAGGTCGTCTCCGTC";
        private const string UnitigR = "CTCCGTCAGCTTGTTTGGAGCAGA";
        private const string UnitigC = "GAGCAGAGTCGTTCTGCGAGG";
        private const string UnitigD = "GAGCAGACCGTCTGTAACAGC";

        private static (List<string> UnitigList, Dictionary<KmerKey, (int UnitigId, int Position)> KmerDict) Build(
            params string[] unitigs)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };
            List<string> unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> kmerDict = [];

            var id = 1;
            foreach (var seq in unitigs)
            {
                unitigList.Add(seq);
                unitigList.Add(Util.V_逆相補(seq));
                for (var i = K; i <= seq.Length; i++)
                {
                    var startPos = i - K;
                    var key = new KmerKey(seq.AsSpan(startPos, K));
                    Register(kmerDict, key, id, startPos);
                    Register(kmerDict, key.Get_逆相補(), -id, seq.Length - i);
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
        public void ResolveShortRepeats_UsesSpanningPairs_ToSplitTheRepeatIntoTwoCleanPaths()
        {
            var (unitigList, kmerDict) = Build(UnitigA, UnitigB, UnitigR, UnitigC, UnitigD);
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);

            var a = ContigMaker.Get_頂点番号(1);
            var b = ContigMaker.Get_頂点番号(2);
            var r = ContigMaker.Get_頂点番号(3);
            var c = ContigMaker.Get_頂点番号(4);
            var d = ContigMaker.Get_頂点番号(5);

            // 前提: R が入次数2・出次数2の反復として構築されている。
            Assert.Equal(2, graph.A_出辺[r].Count);
            Assert.Equal(2, graph.Get_入次数(r));

            // A-C と B-D を跨いだペアだけが観測された、という証拠を与える。
            Dictionary<(int, int), ulong> pairLink = new()
            {
                [(a, c)] = 30,
                [(b, d)] = 28,
            };
            Dictionary<(int, int), ulong> support = [];
            var vertexCountBefore = graph.A_出辺.Count;

            var resolved = graph.V_解決_短い反復(
                unitigList, support, pairLink, p_反復長の上限: 500, p_優勢閾値: 0.8m, p_最小証拠数: 5);

            Assert.Equal(1, resolved);

            // 反復が複製され、頂点が2つ(順鎖・逆鎖)増えているはず。
            Assert.Equal(vertexCountBefore + 2, graph.A_出辺.Count);
            Assert.Equal(UnitigR, unitigList[vertexCountBefore]);

            // どちらの経路も「反復のコピーを1つだけ通る一本道」になっていること。
            // 元の頂点 r と複製のどちらが A 側に残るかは辺の格納順に依存する
            // (unitig の番号とは無関係)ので、頂点の同一性ではなく
            // 経路の構造と対応付けを検証する。
            AssertUntangled(graph, unitigList, from: a, to: c, otherFrom: b, otherTo: d);
        }

        /// <summary>
        /// from →(反復のコピー)→ to と otherFrom →(別のコピー)→ otherTo という
        /// 2本の独立した一本道になっていることを検証する。
        /// </summary>
        private static void AssertUntangled(UnitigGraph graph, List<string> unitigList, int from, int to, int otherFrom, int otherTo)
        {
            var viaFirst = Assert.Single(graph.A_出辺[from]);
            var viaSecond = Assert.Single(graph.A_出辺[otherFrom]);

            // それぞれ別のコピーを通ること(同じ頂点を共有していたら解けていない)。
            Assert.NotEqual(viaFirst, viaSecond);

            // 通る頂点はどちらも反復配列そのもの。
            Assert.Equal(UnitigR, unitigList[viaFirst]);
            Assert.Equal(UnitigR, unitigList[viaSecond]);

            // 各コピーは入次数1・出次数1の一本道。
            Assert.Equal(1, graph.Get_入次数(viaFirst));
            Assert.Equal(1, graph.Get_入次数(viaSecond));
            Assert.Equal([to], graph.A_出辺[viaFirst]);
            Assert.Equal([otherTo], graph.A_出辺[viaSecond]);

            // 逆鎖側も対称であること(片側だけ付け替えるとグラフが壊れ、
            // 順鎖と逆鎖で別々の経路が組まれてしまう)。
            Assert.Contains(viaFirst ^ 1, graph.A_出辺[to ^ 1]);
            Assert.Contains(viaSecond ^ 1, graph.A_出辺[otherTo ^ 1]);
            Assert.Contains(from ^ 1, graph.A_出辺[viaFirst ^ 1]);
            Assert.Contains(otherFrom ^ 1, graph.A_出辺[viaSecond ^ 1]);
        }

        [Fact]
        public void ResolveShortRepeats_HonoursTheCrossedPairing_WhenThatIsWhatThePairsShow()
        {
            var (unitigList, kmerDict) = Build(UnitigA, UnitigB, UnitigR, UnitigC, UnitigD);
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);

            var a = ContigMaker.Get_頂点番号(1);
            var b = ContigMaker.Get_頂点番号(2);
            var c = ContigMaker.Get_頂点番号(4);
            var d = ContigMaker.Get_頂点番号(5);

            // 今度は A-D と B-C の組み合わせが支持されている。
            Dictionary<(int, int), ulong> pairLink = new()
            {
                [(a, d)] = 25,
                [(b, c)] = 31,
            };
            Dictionary<(int, int), ulong> support = [];

            var resolved = graph.V_解決_短い反復(
                unitigList, support, pairLink, p_反復長の上限: 500, p_優勢閾値: 0.8m, p_最小証拠数: 5);

            Assert.Equal(1, resolved);
            // 交差した対応付け: A は D へ、B は C へ繋がる。
            AssertUntangled(graph, unitigList, from: a, to: d, otherFrom: b, otherTo: c);
        }

        /// <summary>
        /// 両方の対応付けが同程度に支持されている場合、どちらが正しいか判断できない。
        /// 誤った繋ぎ方は誤アセンブリを生むため、繋がずに残すのが正しい。
        /// </summary>
        [Fact]
        public void ResolveShortRepeats_LeavesTheRepeatAlone_WhenPairsDoNotFavourEitherPairing()
        {
            var (unitigList, kmerDict) = Build(UnitigA, UnitigB, UnitigR, UnitigC, UnitigD);
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);

            var a = ContigMaker.Get_頂点番号(1);
            var b = ContigMaker.Get_頂点番号(2);
            var r = ContigMaker.Get_頂点番号(3);
            var c = ContigMaker.Get_頂点番号(4);
            var d = ContigMaker.Get_頂点番号(5);

            Dictionary<(int, int), ulong> pairLink = new()
            {
                [(a, c)] = 15,
                [(b, d)] = 14,
                [(a, d)] = 13,
                [(b, c)] = 16,
            };
            Dictionary<(int, int), ulong> support = [];
            var vertexCountBefore = graph.A_出辺.Count;

            var resolved = graph.V_解決_短い反復(
                unitigList, support, pairLink, p_反復長の上限: 500, p_優勢閾値: 0.8m, p_最小証拠数: 5);

            Assert.Equal(0, resolved);
            Assert.Equal(vertexCountBefore, graph.A_出辺.Count);
            Assert.Equal(2, graph.A_出辺[r].Count);
            Assert.Equal(2, graph.Get_入次数(r));
        }

        /// <summary>
        /// フラグメントで跨げない長さの反復は、そもそも証拠が得られないので対象外。
        /// (跨げていないのに偶然の対応付けで繋ぐと誤アセンブリになる。)
        /// </summary>
        [Fact]
        public void ResolveShortRepeats_SkipsRepeatsLongerThanTheFragmentCanSpan()
        {
            var (unitigList, kmerDict) = Build(UnitigA, UnitigB, UnitigR, UnitigC, UnitigD);
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);

            var a = ContigMaker.Get_頂点番号(1);
            var b = ContigMaker.Get_頂点番号(2);
            var c = ContigMaker.Get_頂点番号(4);
            var d = ContigMaker.Get_頂点番号(5);

            Dictionary<(int, int), ulong> pairLink = new()
            {
                [(a, c)] = 30,
                [(b, d)] = 28,
            };
            Dictionary<(int, int), ulong> support = [];
            var vertexCountBefore = graph.A_出辺.Count;

            // R は24bp なので、上限を10bp にすれば対象外になる。
            var resolved = graph.V_解決_短い反復(
                unitigList, support, pairLink, p_反復長の上限: 10, p_優勢閾値: 0.8m, p_最小証拠数: 5);

            Assert.Equal(0, resolved);
            Assert.Equal(vertexCountBefore, graph.A_出辺.Count);
        }

        /// <summary>
        /// 跨いだペアが少なすぎる場合も、偶然の一致で繋いでしまわないよう見送る。
        /// </summary>
        [Fact]
        public void ResolveShortRepeats_SkipsWhenSpanningPairsAreTooFew()
        {
            var (unitigList, kmerDict) = Build(UnitigA, UnitigB, UnitigR, UnitigC, UnitigD);
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);

            var a = ContigMaker.Get_頂点番号(1);
            var c = ContigMaker.Get_頂点番号(4);

            Dictionary<(int, int), ulong> pairLink = new() { [(a, c)] = 2 };
            Dictionary<(int, int), ulong> support = [];
            var vertexCountBefore = graph.A_出辺.Count;

            var resolved = graph.V_解決_短い反復(
                unitigList, support, pairLink, p_反復長の上限: 500, p_優勢閾値: 0.8m, p_最小証拠数: 10);

            Assert.Equal(0, resolved);
            Assert.Equal(vertexCountBefore, graph.A_出辺.Count);
        }
    }
}
