using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 提案F: 短い反復解決に対する r-mer 拒否権(ABySS RResolver型)の検証。
    ///
    /// 反復配列 R が A→R→C, B→R→D という文脈を持つ場合(RepeatResolutionTests
    /// と同じ構造)、A-R・R-D(逆に言えば B-R・R-C も)はどちらの対応付けを
    /// 検証する場合でも de Bruijn グラフ上の本物の辺であり、正しい方の
    /// 組み合わせ(A-R-C, B-R-D)のリードだけからも個々の接合点の存在は
    /// 独立に確認できてしまう。したがって「個々の接合点が存在するか」の
    /// 確認だけでは、A-R-C/B-R-D と A-R-D/B-R-C のどちらの対応付けが
    /// 正しいかを区別できない -- これは実装の欠陥ではなく、反復配列が
    /// まさに「局所的な文脈だけでは区別できない」ことの裏返しである
    /// (区別できるならそもそも反復として1頂点に潰れていない)。
    ///
    /// この拒否権が実際に効くのは、ペア支持が示す対応付けについて
    /// 個々の接合点すら生リードに一切裏付けられない(=そもそもその
    /// unitig 同士が隣接している根拠が生データに無い、破損したデータや
    /// 完全に的外れなペア支持を想定)場合である。
    /// </summary>
    public class RepeatRMerVetoTests
    {
        private const int AmbiguousKmer = int.MinValue;
        private const int K = 8;
        private const int R = K + 10; // AssemblyPipeline の既定の決め方(k + rMer長のk超過分の既定値)と同じ

        // RepeatResolutionTests と同じ構成: A/B が R の先頭 k-1 塩基を、
        // C/D が R の末尾 k-1 塩基を共有する。真の文脈は A→R→C, B→R→D。
        private const string UnitigA = "ACAGTTCGCGAGCCCTCCGTC";
        private const string UnitigB = "TGTATTGAGGTCGTCTCCGTC";
        private const string UnitigR = "CTCCGTCAGCTTGTTTGGAGCAGA";
        private const string UnitigC = "GAGCAGAGTCGTTCTGCGAGG";
        private const string UnitigD = "GAGCAGACCGTCTGTAACAGC";

        private static (List<string> UnitigList, Dictionary<KmerKey, (int UnitigId, int Position)> KmerDict) Build()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };
            List<string> unitigList = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int UnitigId, int Position)> kmerDict = [];

            var id = 1;
            foreach (var seq in new[] { UnitigA, UnitigB, UnitigR, UnitigC, UnitigD })
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

        /// <summary>
        /// 比較対象: 検証器を渡さなければ、RepeatResolutionTests と同じく
        /// ペア支持のみで交差した対応付けが解決される
        /// (=r-mer検証が無いと誤った複製を止められないことの確認)。
        /// </summary>
        [Fact]
        public void WithoutVerifier_TheMisleadingCrossedPairingIsResolvedByPairSupportAlone()
        {
            var (unitigList, kmerDict) = Build();
            var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);

            var a = ContigMaker.Get_頂点番号(1);
            var b = ContigMaker.Get_頂点番号(2);
            var c = ContigMaker.Get_頂点番号(4);
            var d = ContigMaker.Get_頂点番号(5);

            // 現実には A-C, B-D しか読まれていないのに、何らかの理由で
            // ペア支持は交差した組み合わせ(A-D, B-C)を優勢だと示している
            // (誤マッピング等を想定した意図的な誤情報)。
            Dictionary<(int, int), ulong> pairLink = new() { [(a, d)] = 25, [(b, c)] = 31 };
            Dictionary<(int, int), ulong> support = [];

            var resolved = graph.V_解決_短い反復(
                unitigList, support, pairLink, p_反復長の上限: 500, p_優勢閾値: 0.8m, p_最小証拠数: 5);

            Assert.Equal(1, resolved);
        }

        /// <summary>
        /// 拒否権の実際の効き目には限界がある: head-repeat・repeat-tail は
        /// どちらのペアリングでも de Bruijn グラフ上の本物の辺であるため、
        /// 「A-R が存在する」「R-D が存在する」という個々の接合点の確認は
        /// A-R-C / B-R-D という(こちらが正しい)組み合わせのリードだけからも
        /// 満たされてしまう(A-R は A-R-C 由来のリードで、R-D は B-R-D 由来の
        /// リードで、それぞれ独立に確認できてしまうため)。個々の接合点の
        /// 存在確認だけでは、どちらの対応付けが正しいかを区別する情報には
        /// ならない -- これは実装の欠陥ではなく、反復配列がまさに
        /// 「局所的な文脈だけでは区別できない」ことの裏返しである。
        ///
        /// この拒否権が実際に効くのは、ペア支持が示す対応付けについて
        /// 個々の接合点すら生リードに一切裏付けられない(=そもそも
        /// その unitig 同士が本当に隣接している根拠が生データに無い)場合。
        /// ここでは生リードを一切与えず、ペア支持だけで対応付けようとしても
        /// 拒否されることを確認する。
        /// </summary>
        [Fact]
        public void WithVerifier_VetoesTheWinningPairing_WhenNoRawReadDataConfirmsEitherJunctionAtAll()
        {
            var tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_rmer_veto_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(tempDir);
            try
            {
                var (unitigList, kmerDict) = Build();
                var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);

                var a = ContigMaker.Get_頂点番号(1);
                var b = ContigMaker.Get_頂点番号(2);
                var r = ContigMaker.Get_頂点番号(3);
                var c = ContigMaker.Get_頂点番号(4);
                var d = ContigMaker.Get_頂点番号(5);

                Dictionary<(int, int), ulong> pairLink = new() { [(a, d)] = 25, [(b, c)] = 31 };
                Dictionary<(int, int), ulong> support = [];
                var vertexCountBefore = graph.A_出辺.Count;

                // 生データを一切与えない(=マッピング元のリードが実在しない、
                // 破損したデータセット等を想定)。
                var emptyPath = Path.Combine(tempDir, "empty.fq");
                File.WriteAllText(emptyPath, string.Empty);
                var verifier = RepeatRMerVerifier.V_構築([emptyPath, string.Empty], R);

                var resolved = graph.V_解決_短い反復(
                    unitigList, support, pairLink, p_反復長の上限: 500, p_優勢閾値: 0.8m, p_最小証拠数: 5,
                    p_r_mer検証器: verifier);

                Assert.Equal(0, resolved);
                Assert.Equal(vertexCountBefore, graph.A_出辺.Count);
                // 反復は解決されず、入次数2・出次数2のまま残る。
                Assert.Equal(2, graph.A_出辺[r].Count);
                Assert.Equal(2, graph.Get_入次数(r));
            }
            finally
            {
                Directory.Delete(tempDir, recursive: true);
            }
        }

        /// <summary>
        /// 対照実験: 実際に A-D, B-C を跨いだリードがあれば
        /// (=これが真の文脈である場合)、r-mer 検証器を渡しても
        /// 複製は正しく実行される(拒否権は正しい対応付けまで妨げない)。
        /// </summary>
        [Fact]
        public void WithVerifier_StillResolves_WhenTheCrossedPairingIsActuallyReal()
        {
            var tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_rmer_veto_real_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(tempDir);
            try
            {
                var (unitigList, kmerDict) = Build();
                var graph = UnitigGraph.Get_グラフ(unitigList, kmerDict, K, AmbiguousKmer);

                var a = ContigMaker.Get_頂点番号(1);
                var b = ContigMaker.Get_頂点番号(2);
                var c = ContigMaker.Get_頂点番号(4);
                var d = ContigMaker.Get_頂点番号(5);

                Dictionary<(int, int), ulong> pairLink = new() { [(a, d)] = 25, [(b, c)] = 31 };
                Dictionary<(int, int), ulong> support = [];

                // 今度は実際に A→R→D, B→R→C を跨ぐリードを用意する
                // (=交差した対応付けが真の文脈であるケース)。
                var walkAD = UnitigA + UnitigR[(K - 1)..] + UnitigD[(K - 1)..];
                var walkBC = UnitigB + UnitigR[(K - 1)..] + UnitigC[(K - 1)..];
                List<string> reads = [];
                foreach (var walk in new[] { walkAD, walkBC })
                {
                    for (var i = 0; i + 25 <= walk.Length; i++)
                    {
                        reads.Add(walk.Substring(i, 25));
                    }
                }
                var path = Path.Combine(tempDir, "reads.fq");
                using (var writer = new StreamWriter(path))
                {
                    var idx = 0;
                    foreach (var seq in reads)
                    {
                        writer.WriteLine($"@r{idx}");
                        writer.WriteLine(seq);
                        writer.WriteLine("+");
                        writer.WriteLine(new string('I', seq.Length));
                        idx++;
                    }
                }
                var verifier = RepeatRMerVerifier.V_構築([path, string.Empty], R);

                var resolved = graph.V_解決_短い反復(
                    unitigList, support, pairLink, p_反復長の上限: 500, p_優勢閾値: 0.8m, p_最小証拠数: 5,
                    p_r_mer検証器: verifier);

                Assert.Equal(1, resolved);
            }
            finally
            {
                Directory.Delete(tempDir, recursive: true);
            }
        }
    }
}
