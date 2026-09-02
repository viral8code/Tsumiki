using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    public class ContigMakerAdjacencyTests
    {
        /// <summary>
        /// 旧実装(unitigCount^2 の総当り)の挙動をそのまま再現したリファレンス実装。
        /// BuildAdjacencyList(kmerPath を1回だけ走査する新実装)の結果が
        /// これと(頂点ごとの多重集合として)一致することを確認する。
        /// </summary>
        private static List<List<(int To, ulong Count)>> BruteForceReference(
            IReadOnlyDictionary<(int, int), ulong> kmerPath, int unitigCount, int vertexCount)
        {
            var adjacencyList = new List<List<(int, ulong)>>();
            for (var i = 0; i < vertexCount; i++)
            {
                adjacencyList.Add([]);
            }

            for (var i = 1; i <= unitigCount; i++)
            {
                for (var j = 1; j <= unitigCount; j++)
                {
                    if (i == j)
                    {
                        continue;
                    }
                    if (kmerPath.TryGetValue((i, j), out var count))
                    {
                        adjacencyList[i << 1].Add((j << 1, count));
                    }
                    if (kmerPath.TryGetValue((i, -j), out count))
                    {
                        adjacencyList[i << 1].Add(((j << 1) | 1, count));
                    }
                    if (kmerPath.TryGetValue((-i, j), out count))
                    {
                        adjacencyList[(i << 1) | 1].Add((j << 1, count));
                    }
                    if (kmerPath.TryGetValue((-i, -j), out count))
                    {
                        adjacencyList[(i << 1) | 1].Add(((j << 1) | 1, count));
                    }
                }
            }
            return adjacencyList;
        }

        private static void AssertSameMultiset(List<(int To, ulong Count)> expected, List<(int, ulong)> actual)
        {
            Assert.Equal(
                expected.OrderBy(e => e.To).ThenBy(e => e.Count),
                actual.OrderBy(e => e.Item1).ThenBy(e => e.Item2));
        }

        [Fact]
        public void BuildAdjacencyList_MatchesBruteForceReference_OnSmallSyntheticGraph()
        {
            const int unitigCount = 8;
            var vertexCount = 2 + (2 * unitigCount);

            // 順鎖・逆鎖・自己ループを含む小規模なテストケース。
            // 注: (id, -id) のような「符号違いだが絶対値が同じ」ペアは、
            // 旧実装(i,jとも正の値で回すループでi==jをスキップする)では
            // そもそもi==jの時点でスキップされ絶対に検出できない既知の
            // 別バグがあるため、ここでは含めない(そちらは別テストで検証する)。
            Dictionary<(int, int), ulong> kmerPath = new()
            {
                [(1, 2)] = 5,
                [(2, 3)] = 3,
                [(1, -4)] = 7,
                [(-3, 5)] = 2,
                [(-5, -6)] = 9,
                [(6, 7)] = 1,
                [(7, 8)] = 4,
                [(4, 4)] = 100, // from==to は除外される
            };

            var expected = BruteForceReference(kmerPath, unitigCount, vertexCount);
            var actual = ContigMaker.BuildAdjacencyList(kmerPath, vertexCount);

            Assert.Equal(expected.Count, actual.Count);
            for (var v = 0; v < vertexCount; v++)
            {
                AssertSameMultiset(expected[v], actual[v]);
            }
        }

        [Fact]
        public void BuildAdjacencyList_EmptyKmerPath_ProducesAllEmptyVertices()
        {
            var actual = ContigMaker.BuildAdjacencyList(new Dictionary<(int, int), ulong>(), vertexCount: 12);

            Assert.Equal(12, actual.Count);
            Assert.All(actual, list => Assert.Empty(list));
        }

        [Fact]
        public void BuildAdjacencyList_SkipsSelfLoopsWithIdenticalSignedId()
        {
            Dictionary<(int, int), ulong> kmerPath = new()
            {
                [(3, 3)] = 42,
            };

            var actual = ContigMaker.BuildAdjacencyList(kmerPath, vertexCount: 12);

            Assert.All(actual, list => Assert.Empty(list));
        }

        [Theory]
        [InlineData(1, 2)]
        [InlineData(-1, 3)]
        [InlineData(5, 10)]
        [InlineData(-5, 11)]
        public void VertexIndex_EncodesMagnitudeAndSign(int signedId, int expected)
        {
            Assert.Equal(expected, ContigMaker.VertexIndex(signedId));
        }

        /// <summary>
        /// 旧実装(i,jとも1..unitigCountの正の値で回す二重ループでi==jをスキップする)は、
        /// (id, -id)のような「絶対値は同じだが符号違い」の隣接ペアを、i==jの時点で
        /// 丸ごとスキップしてしまい絶対に検出できなかった(from!=toではあるものの、
        /// ループがそこへ到達すらしない既存のバグ)。新実装はkmerPathのキーを
        /// そのまま1回走査するだけなので、このような正当なエッジも正しく拾える。
        /// </summary>
        [Fact]
        public void BuildAdjacencyList_IncludesSameMagnitudeOppositeSignEdges_FixingOldBlindSpot()
        {
            Dictionary<(int, int), ulong> kmerPath = new()
            {
                [(2, -2)] = 6,
            };

            var actual = ContigMaker.BuildAdjacencyList(kmerPath, vertexCount: 12);

            var fromVertex = ContigMaker.VertexIndex(2);
            var toVertex = ContigMaker.VertexIndex(-2);
            Assert.Contains((toVertex, 6UL), actual[fromVertex]);
        }

        /// <summary>
        /// 旧実装(unitigCount^2)なら未エラー訂正の実データ規模(unitig数が上限
        /// 10万に近い)で現実的な時間に終わらなかった規模を、新実装が
        /// 数秒以内に処理できることを確認する回帰テスト。
        /// </summary>
        [Fact]
        public void BuildAdjacencyList_HandlesLargeSparseGraph_QuicklyAndCorrectly()
        {
            const int unitigCount = 90_000;
            var vertexCount = 2 + (2 * unitigCount);
            var random = new Random(12345);

            Dictionary<(int, int), ulong> kmerPath = new();
            for (var i = 0; i < 150_000; i++)
            {
                var from = random.Next(1, unitigCount + 1) * (random.Next(2) == 0 ? 1 : -1);
                int to;
                do
                {
                    to = random.Next(1, unitigCount + 1) * (random.Next(2) == 0 ? 1 : -1);
                } while (to == from);
                kmerPath[(from, to)] = (ulong)random.Next(1, 100);
            }

            var sw = System.Diagnostics.Stopwatch.StartNew();
            var actual = ContigMaker.BuildAdjacencyList(kmerPath, vertexCount);
            sw.Stop();

            Assert.True(sw.Elapsed < TimeSpan.FromSeconds(5), $"BuildAdjacencyList took too long: {sw.Elapsed}");

            var totalEdges = actual.Sum(list => list.Count);
            Assert.Equal(kmerPath.Count, totalEdges);
        }
    }
}
