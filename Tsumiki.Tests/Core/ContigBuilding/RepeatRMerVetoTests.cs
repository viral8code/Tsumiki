using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 提案 F: 短い反復解決に対する r-mer 拒否権 (ABySS RResolver 型) の検証
    /// </summary>
    /// <remarks>
    /// 反復配列 R が A→R→C, B→R→D という文脈を持つ場合 (RepeatResolutionTests と同じ構造)、A-R ・ R-D (逆に言えば B-R ・ R-C も) はどちらの対応付けを検証する場合でも de Bruijn グラフ上の本物の辺であり、正しい方の組み合わせ (A-R-C, B-R-D) のリードだけからも個々の接合点の存在は独立に確認できてしまう<br/>
    /// したがって「個々の接合点が存在するか」の確認だけでは、A-R-C/B-R-D と A-R-D/B-R-C のどちらの対応付けが正しいかを区別できない -- これは実装の欠陥ではなく、反復配列がまさに「局所的な文脈だけでは区別できない」ことの裏返しである (区別できるならそもそも反復として 1 頂点に潰れていない) <br/>
    /// この拒否権が実際に効くのは、ペア支持が示す対応付けについて個々の接合点すら生リードに一切裏付けられない (=そもそもその unitig 同士が隣接している根拠が生データに無い、破損したデータや完全に的外れなペア支持を想定) 場合である
    /// </remarks>
    public class RepeatRMerVetoTests
    {
        #region 定数

        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int 曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int k長 = 8;

        /// <summary>
        /// この検証で使う r-mer 長
        /// </summary>
        private const int r長 = k長 + 10; // AssemblyPipeline の既定の決め方 (k + rMer 長の k 超過分の既定値) と同じ

        // RepeatResolutionTests と同じ構成: A/B が R の先頭 k-1 塩基を、
        // C/D が R の末尾 k-1 塩基を共有する
        // 真の文脈は A→R→C, B→R→D

        /// <summary>
        /// 反復の手前にある片方の入口
        /// </summary>
        private const string 入口ユニティグ = "ACAGTTCGCGAGCCCTCCGTC";

        /// <summary>
        /// 反復の手前にあるもう片方の入口
        /// </summary>
        private const string 代替入口ユニティグ = "TGTATTGAGGTCGTCTCCGTC";

        /// <summary>
        /// 入口と出口に挟まれた反復配列
        /// </summary>
        private const string 反復ユニティグ = "CTCCGTCAGCTTGTTTGGAGCAGA";

        /// <summary>
        /// 反復の先にある片方の出口
        /// </summary>
        private const string 出口ユニティグ = "GAGCAGAGTCGTTCTGCGAGG";

        /// <summary>
        /// 反復の先にあるもう片方の出口
        /// </summary>
        private const string 代替出口ユニティグ = "GAGCAGACCGTCTGTAACAGC";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 比較対象: 検証器を渡さなければ、RepeatResolutionTests と同じくペア支持のみで交差した対応付けが解決される (=r-mer 検証が無いと誤った複製を止められないことの確認)
        /// </summary>
        [Fact]
        public void V_WithoutVerifier_TheMisleadingCrossedPairingIsResolvedByPairSupportAlone()
        {
            var (l_ユニティグ一覧, l_kmer辞書) = V_構築();
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);

            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            var l_終端配列 = ContigMaker.Get_頂点番号(5);

            // 現実には A-C, B-D しか読まれていないのに、何らかの理由で
            // ペア支持は交差した組み合わせ (A-D, B-C) を優勢だと示している
            // (誤マッピング等を想定した意図的な誤情報)
            Dictionary<(int, int), ulong> l_ペア連結 = new() { [(l_先頭配列, l_終端配列)] = 25UL, [(l_中間配列, l_末尾配列)] = 31UL };
            Dictionary<(int, int), ulong> l_支持 = [];

            var l_解決数 = l_グラフ.V_解決_短い反復(l_ユニティグ一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);

            Assert.Equal(1, l_解決数);
        }

        /// <summary>
        /// 拒否権の効き目には限界がある
        /// </summary>
        /// <remarks>
        /// head-repeat と repeat-tail はどちらのペアリングでも de Bruijn グラフ上の本物の辺なので、A-R が存在する、R-D が存在するという個々の接合点の確認は、正しい組み合わせである A-R-C と B-R-D のリードだけからも満たされてしまう<br/>
        /// A-R は A-R-C 由来のリードで、R-D は B-R-D 由来のリードで、それぞれ独立に確認できるため
        /// </remarks>
        /// <remarks>
        /// 個々の接合点の存在確認だけでは、どちらの対応付けが正しいかを区別する情報にはならない -- これは実装の欠陥ではなく、反復配列がまさに「局所的な文脈だけでは区別できない」ことの裏返しである<br/>
        /// この拒否権が実際に効くのは、ペア支持が示す対応付けについて個々の接合点すら生リードに一切裏付けられない (=そもそもその unitig 同士が本当に隣接している根拠が生データに無い) 場合<br/>
        /// ここでは生リードを一切与えず、ペア支持だけで対応付けようとしても拒否されることを確認する
        /// </remarks>
        [Fact]
        public void V_WithVerifier_VetoesTheWinningPairing_WhenNoRawReadDataConfirmsEitherJunctionAtAll()
        {
            var l_作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_rmer_veto_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_作業ディレクトリ);
            try
            {
                var (l_ユニティグ一覧, l_kmer辞書) = V_構築();
                var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);

                var l_先頭配列 = ContigMaker.Get_頂点番号(1);
                var l_中間配列 = ContigMaker.Get_頂点番号(2);
                var l_反復配列 = ContigMaker.Get_頂点番号(3);
                var l_末尾配列 = ContigMaker.Get_頂点番号(4);
                var l_終端配列 = ContigMaker.Get_頂点番号(5);

                Dictionary<(int, int), ulong> l_ペア連結 = new() { [(l_先頭配列, l_終端配列)] = 25UL, [(l_中間配列, l_末尾配列)] = 31UL };
                Dictionary<(int, int), ulong> l_支持 = [];
                var l_変更前頂点数 = l_グラフ.A_出辺.Count;

                // 生データを一切与えない (=マッピング元のリードが実在しない、
                // 破損したデータセット等を想定)
                var l_空ファイルパス = Path.Combine(l_作業ディレクトリ, "empty.fq");
                File.WriteAllText(l_空ファイルパス, string.Empty);
                var l_検証器 = RepeatRMerVerifier.V_構築([l_空ファイルパス, string.Empty], r長);

                var l_解決数 = l_グラフ.V_解決_短い反復(l_ユニティグ一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL, p_r_mer検証器: l_検証器);

                Assert.Equal(0, l_解決数);
                Assert.Equal(l_変更前頂点数, l_グラフ.A_出辺.Count);
                // 反復は解決されず、入次数 2 ・出次数 2 のまま残る
                Assert.Equal(2, l_グラフ.A_出辺[l_反復配列].Count);
                Assert.Equal(2, l_グラフ.Get_入次数(l_反復配列));
            }
            finally
            {
                Directory.Delete(l_作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 対照実験: 実際に A-D, B-C を跨いだリードがあれば (=これが真の文脈である場合)、r-mer 検証器を渡しても複製は正しく実行される (拒否権は正しい対応付けまで妨げない)
        /// </summary>
        [Fact]
        public void V_WithVerifier_StillResolves_WhenTheCrossedPairingIsActuallyReal()
        {
            var l_作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_rmer_veto_real_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(l_作業ディレクトリ);
            try
            {
                var (l_ユニティグ一覧, l_kmer辞書) = V_構築();
                var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);

                var l_先頭配列 = ContigMaker.Get_頂点番号(1);
                var l_中間配列 = ContigMaker.Get_頂点番号(2);
                var l_末尾配列 = ContigMaker.Get_頂点番号(4);
                var l_終端配列 = ContigMaker.Get_頂点番号(5);

                Dictionary<(int, int), ulong> l_ペア連結 = new() { [(l_先頭配列, l_終端配列)] = 25UL, [(l_中間配列, l_末尾配列)] = 31UL };
                Dictionary<(int, int), ulong> l_支持 = [];

                // 今度は実際に A→R→D, B→R→C を跨ぐリードを用意する
                // (=交差した対応付けが真の文脈であるケース)
                var l_主経路 = 入口ユニティグ + 反復ユニティグ[(k長 - 1)..] + 代替出口ユニティグ[(k長 - 1)..];
                var l_代替経路 = 代替入口ユニティグ + 反復ユニティグ[(k長 - 1)..] + 出口ユニティグ[(k長 - 1)..];
                List<string> l_リード群 = [];
                foreach (var l_経路 in new[] { l_主経路, l_代替経路 })
                {
                    for (var i = 0; i + 25 <= l_経路.Length; i++)
                    {
                        l_リード群.Add(l_経路.Substring(i, 25));
                    }
                }
                var l_パス = Path.Combine(l_作業ディレクトリ, "reads.fq");
                using (var l_書き込み = new StreamWriter(l_パス))
                {
                    var l_位置 = 0;
                    foreach (var l_配列 in l_リード群)
                    {
                        l_書き込み.WriteLine($"@r{l_位置}");
                        l_書き込み.WriteLine(l_配列);
                        l_書き込み.WriteLine("+");
                        l_書き込み.WriteLine(new string('I', l_配列.Length));
                        l_位置++;
                    }
                }
                var l_検証器 = RepeatRMerVerifier.V_構築([l_パス, string.Empty], r長);

                var l_解決数 = l_グラフ.V_解決_短い反復(l_ユニティグ一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL, p_r_mer検証器: l_検証器);

                Assert.Equal(1, l_解決数);
            }
            finally
            {
                Directory.Delete(l_作業ディレクトリ, recursive: true);
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// A ・ B ・ R ・ C ・ D のユニティグ配列と、それらから作った kmer 辞書を組み立てる
        /// </summary>
        /// <returns></returns>
        private static (List<string> A_ユニティグ一覧, Dictionary<KmerKey, (int A_ユニティグID, int A_位置)> A_kmer辞書) V_構築()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            List<string> l_ユニティグ配列 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_ユニティグID, int A_位置)> l_kmer辞書 = [];

            var l_ID = 1;
            foreach (var l_配列 in new[] { 入口ユニティグ, 代替入口ユニティグ, 反復ユニティグ, 出口ユニティグ, 代替出口ユニティグ })
            {
                l_ユニティグ配列.Add(l_配列);
                l_ユニティグ配列.Add(Util.V_逆相補(l_配列));
                for (var i = k長; i <= l_配列.Length; i++)
                {
                    var l_開始位置 = i - k長;
                    var l_キー = new KmerKey(l_配列.AsSpan(l_開始位置, k長));
                    V_登録(l_kmer辞書, l_キー, l_ID, l_開始位置);
                    V_登録(l_kmer辞書, l_キー.Get_逆相補(), -l_ID, l_配列.Length - i);
                }
                l_ID++;
            }
            return (l_ユニティグ配列, l_kmer辞書);
        }

        /// <summary>
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_キー">登録する k-mer</param>
        /// <param name="p_ID">ユニティグ ID</param>
        /// <param name="p_位置">ユニティグ内の開始位置</param>
        private static void V_登録(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存))
            {
                if (l_既存.Item1 is 曖昧kmer番号 || l_既存.Item1 == p_ID)
                {
                    return;
                }
                p_辞書[p_キー] = (曖昧kmer番号, 0);
                return;
            }
            p_辞書[p_キー] = (p_ID, p_位置);
        }

        #endregion

    }
}
