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
    public class RepeatRMerVetoTests
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki rmer veto tests
        /// </summary>
        private const string C_項目_tsumiki_rmer_veto_tests = "tsumiki_rmer_veto_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 empty fq
        /// </summary>
        private const string C_ファイル名_empty_fq = "empty.fq";

        /// <summary>
        /// 項目 tsumiki rmer veto real tests
        /// </summary>
        private const string C_項目_tsumiki_rmer_veto_real_tests = "tsumiki_rmer_veto_real_tests_";

        /// <summary>
        /// ファイル名 reads fq
        /// </summary>
        private const string C_ファイル名_reads_fq = "reads.fq";

        /// <summary>
        /// FASTQ 品質区切り
        /// </summary>
        private const string C_FASTQ品質区切り = "+";

        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int C_曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k長 = 8;

        /// <summary>
        /// この検証で使う r-mer 長
        /// </summary>
        private const int C_r長 = C_k長 + 10;

        /// <summary>
        /// 反復の手前にある片方の入口
        /// </summary>
        private const string C_入口unitig = "ACAGTTCGCGAGCCCTCCGTC";

        /// <summary>
        /// 反復の手前にあるもう片方の入口
        /// </summary>
        private const string C_代替入口unitig = "TGTATTGAGGTCGTCTCCGTC";

        /// <summary>
        /// 入口と出口に挟まれた反復配列
        /// </summary>
        private const string C_反復unitig = "CTCCGTCAGCTTGTTTGGAGCAGA";

        /// <summary>
        /// 反復の先にある片方の出口
        /// </summary>
        private const string C_出口unitig = "GAGCAGAGTCGTTCTGCGAGG";

        /// <summary>
        /// 反復の先にあるもう片方の出口
        /// </summary>
        private const string C_代替出口unitig = "GAGCAGACCGTCTGTAACAGC";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 比較対象: 検証器を渡さなければ、RepeatResolutionTests と同じくペア支持のみで交差した対応付けが解決される (=r-mer 検証が無いと誤った複製を止められないことの確認)
        /// </summary>
        [Fact]
        public void V_WithoutVerifier_TheMisleadingCrossedPairingIsResolvedByPairSupportAlone()
        {
            var (l_unitig一覧, l_kmer辞書) = V_構築();
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号);
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            var l_終端配列 = ContigMaker.Get_頂点番号(5);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 25UL,
                [(l_中間配列, l_末尾配列)] = 31UL
            };
            Dictionary<(int, int), ulong> l_支持 = [];
            var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(1, l_解決数);
        }

        /// <summary>
        /// 拒否権の効き目には限界がある
        /// </summary>
        [Fact]
        public void V_WithVerifier_VetoesTheWinningPairing_WhenNoRawReadDataConfirmsEitherJunctionAtAll()
        {
            var l_作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_rmer_veto_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(l_作業ディレクトリ);
            try
            {
                var (l_unitig一覧, l_kmer辞書) = V_構築();
                var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号);
                var l_先頭配列 = ContigMaker.Get_頂点番号(1);
                var l_中間配列 = ContigMaker.Get_頂点番号(2);
                var l_反復配列 = ContigMaker.Get_頂点番号(3);
                var l_末尾配列 = ContigMaker.Get_頂点番号(4);
                var l_終端配列 = ContigMaker.Get_頂点番号(5);
                Dictionary<(int, int), ulong> l_ペア連結 = new()
                {
                    [(l_先頭配列, l_終端配列)] = 25UL,
                    [(l_中間配列, l_末尾配列)] = 31UL
                };
                Dictionary<(int, int), ulong> l_支持 = [];
                var l_変更前頂点数 = l_グラフ.A_出辺.Count;
                var l_空ファイルパス = Path.Combine(l_作業ディレクトリ, C_ファイル名_empty_fq);
                File.WriteAllText(l_空ファイルパス, string.Empty);
                var l_検証器 = RepeatRMerVerifier.V_構築([l_空ファイルパス, string.Empty], C_r長);
                var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL, p_r_mer検証器: l_検証器);
                Assert.Equal(0, l_解決数);
                Assert.Equal(l_変更前頂点数, l_グラフ.A_出辺.Count);
                Assert.Equal(2, l_グラフ.A_出辺[l_反復配列].Count);
                Assert.Equal(2, l_グラフ.Get_入次数(l_反復配列));
            }
            finally
            {
                Directory.Delete(l_作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 対照実験: 実際に A-D, B-C を跨いだリードがあれば (=これが真の文脈である場合) 、r-mer 検証器を渡しても複製は正しく実行される (拒否権は正しい対応付けまで妨げない)
        /// </summary>
        [Fact]
        public void V_WithVerifier_StillResolves_WhenTheCrossedPairingIsActuallyReal()
        {
            var l_作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_rmer_veto_real_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(l_作業ディレクトリ);
            try
            {
                var (l_unitig一覧, l_kmer辞書) = V_構築();
                var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号);
                var l_先頭配列 = ContigMaker.Get_頂点番号(1);
                var l_中間配列 = ContigMaker.Get_頂点番号(2);
                var l_末尾配列 = ContigMaker.Get_頂点番号(4);
                var l_終端配列 = ContigMaker.Get_頂点番号(5);
                Dictionary<(int, int), ulong> l_ペア連結 = new()
                {
                    [(l_先頭配列, l_終端配列)] = 25UL,
                    [(l_中間配列, l_末尾配列)] = 31UL
                };
                Dictionary<(int, int), ulong> l_支持 = [];
                var l_主経路 = C_入口unitig + C_反復unitig[(C_k長 - 1)..] + C_代替出口unitig[(C_k長 - 1)..];
                var l_代替経路 = C_代替入口unitig + C_反復unitig[(C_k長 - 1)..] + C_出口unitig[(C_k長 - 1)..];
                List<string> l_リード群 = [];
                foreach (var l_経路 in new[]
                {
                    l_主経路,
                    l_代替経路
                }

                )
                {
                    for (var i = 0; i + 25 <= l_経路.Length; i++)
                    {
                        l_リード群.Add(l_経路.Substring(i, 25));
                    }
                }

                var l_パス = Path.Combine(l_作業ディレクトリ, C_ファイル名_reads_fq);
                using (var l_書き込み = new StreamWriter(l_パス))
                {
                    var l_位置 = 0;
                    foreach (var l_配列 in l_リード群)
                    {
                        l_書き込み.WriteLine($"@r{l_位置}");
                        l_書き込み.WriteLine(l_配列);
                        l_書き込み.WriteLine(C_FASTQ品質区切り);
                        l_書き込み.WriteLine(new string('I', l_配列.Length));
                        l_位置++;
                    }
                }

                var l_検証器 = RepeatRMerVerifier.V_構築([l_パス, string.Empty], C_r長);
                var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL, p_r_mer検証器: l_検証器);
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
        /// A ・ B ・ R ・ C ・ D の unitig 配列と、それらから作った kmer 辞書を組み立てる
        /// </summary>
        /// <returns></returns>
        private static (List<string> A_unitig一覧, Dictionary<KmerKey, (int A_unitigID, int A_位置)> A_kmer辞書) V_構築()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
                A_スレッド数 = 1
            };
            List<string> l_unitig配列 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_unitigID, int A_位置)> l_kmer辞書 = [];
            var l_ID = 1;
            foreach (var l_配列 in new[]
            {
                C_入口unitig,
                C_代替入口unitig,
                C_反復unitig,
                C_出口unitig,
                C_代替出口unitig
            }

            )
            {
                l_unitig配列.Add(l_配列);
                l_unitig配列.Add(Util.V_逆相補(l_配列));
                for (var i = C_k長; i <= l_配列.Length; i++)
                {
                    var l_開始位置 = i - C_k長;
                    var l_キー = new KmerKey(l_配列.AsSpan(l_開始位置, C_k長));
                    V_登録(l_kmer辞書, l_キー, l_ID, l_開始位置);
                    V_登録(l_kmer辞書, l_キー.Get_逆相補(), -l_ID, l_配列.Length - i);
                }

                l_ID++;
            }

            return (l_unitig配列, l_kmer辞書);
        }

        /// <summary>
        /// k-mer を、それが載る unitig と開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_キー">登録する k-mer</param>
        /// <param name="p_ID">unitig ID</param>
        /// <param name="p_位置">unitig 内の開始位置</param>
        private static void V_登録(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存))
            {
                if (l_既存.Item1 is C_曖昧kmer番号 || l_既存.Item1 == p_ID)
                {
                    return;
                }

                p_辞書[p_キー] = (C_曖昧kmer番号, 0);
                return;
            }

            p_辞書[p_キー] = (p_ID, p_位置);
        }

        #endregion
    }
}
