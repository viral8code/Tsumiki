using Tsumiki.Commons;
using Tsumiki.Cores.Output;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// unitig グラフの GFA1 出力の検証
    /// </summary>
    public class GfaWriterTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// 項目 tsumiki gfa tests
        /// </summary>
        private const string C_項目_tsumiki_gfa_tests = "tsumiki_gfa_tests_";

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// ファイル名 graph gfa
        /// </summary>
        private const string C_ファイル名_graph_gfa = "graph.gfa";

        /// <summary>
        /// 項目 H VN Z 1 0
        /// </summary>
        private const string C_項目_H_VN_Z_1_0 = "H\tVN:Z:1.0";

        /// <summary>
        /// 項目 S
        /// </summary>
        private const string C_項目_S = "S\t";

        /// <summary>
        /// ファイル名 graph links gfa
        /// </summary>
        private const string C_ファイル名_graph_links_gfa = "graph_links.gfa";

        /// <summary>
        /// 項目 L
        /// </summary>
        private const string C_項目_L = "L\t";

        /// <summary>
        /// ファイル名 graph cn gfa
        /// </summary>
        private const string C_ファイル名_graph_cn_gfa = "graph_cn.gfa";

        /// <summary>
        /// 項目 S 2
        /// </summary>
        private const string C_項目_S_2 = "S\t2\t";

        /// <summary>
        /// 項目 CN i 2
        /// </summary>
        private const string C_項目_CN_i_2 = "CN:i:2";

        /// <summary>
        /// 項目 S 1
        /// </summary>
        private const string C_項目_S_1 = "S\t1\t";

        /// <summary>
        /// 項目 CN i 1
        /// </summary>
        private const string C_項目_CN_i_1 = "CN:i:1";

        /// <summary>
        /// ファイル名 graph nocn gfa
        /// </summary>
        private const string C_ファイル名_graph_nocn_gfa = "graph_nocn.gfa";

        /// <summary>
        /// 項目 CN i
        /// </summary>
        private const string C_項目_CN_i = "CN:i:";

        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int C_曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k長 = 8;

        /// <summary>
        /// 分岐元
        /// </summary>
        private const string C_unitigA = "TGGCAAGTCACTCTCGACCGA";

        /// <summary>
        /// 分岐先の片方
        /// </summary>
        private const string C_unitigB = "CGACCGAACGGCGCCGGATC";

        /// <summary>
        /// 分岐先のもう片方
        /// </summary>
        private const string C_unitigC = "CGACCGACTGTAATTCTACC";

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public GfaWriterTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), C_項目_tsumiki_gfa_tests + Guid.NewGuid().ToString(C_GUID書式));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// unitig ごとに順方向配列と長さを持つ S 行を 1 つ書くことを確かめる
        /// </summary>
        [Fact]
        public void V_出力_unitigごとに順方向配列と長さを持つS行を1つ書く()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_graph_gfa);
            GfaWriter.V_出力(l_パス, l_unitig一覧, l_グラフ, C_k長);
            var l_行 = File.ReadAllLines(l_パス);
            Assert.Equal(C_項目_H_VN_Z_1_0, l_行[0]);
            var l_S行 = l_行.Where(l => l.StartsWith(C_項目_S)).ToList();
            Assert.Equal(3, l_S行.Count);
            Assert.Contains($"S\t1\t{C_unitigA}\tLN:i:{C_unitigA.Length}", l_S行);
            Assert.Contains($"S\t2\t{C_unitigB}\tLN:i:{C_unitigB.Length}", l_S行);
            Assert.Contains($"S\t3\t{C_unitigC}\tLN:i:{C_unitigC.Length}", l_S行);
        }

        /// <summary>
        /// 物理的な隣接ごとに L 行を 1 本だけ書くことを確かめる
        /// </summary>
        [Fact]
        public void V_出力_物理的な隣接ごとにL行を1本だけ書く()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_graph_links_gfa);
            GfaWriter.V_出力(l_パス, l_unitig一覧, l_グラフ, C_k長);
            var l_L行 = File.ReadAllLines(l_パス).Where(l => l.StartsWith(C_項目_L)).ToList();
            Assert.Equal(2, l_L行.Count);
            Assert.Contains(l_L行, l => l == $"L\t1\t+\t2\t+\t{C_k長 - 1}M");
            Assert.Contains(l_L行, l => l == $"L\t1\t+\t3\t+\t{C_k長 - 1}M");
        }

        /// <summary>
        /// コピー数を渡せば CN タグを含めることを確かめる
        /// </summary>
        [Fact]
        public void V_出力_コピー数を渡せばCNタグを含める()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_graph_cn_gfa);
            Dictionary<int, int> l_コピー数 = new()
            {
                [1] = 1,
                [2] = 2,
                [3] = 1
            };
            GfaWriter.V_出力(l_パス, l_unitig一覧, l_グラフ, C_k長, l_コピー数);
            var l_S行 = File.ReadAllLines(l_パス).Where(l => l.StartsWith(C_項目_S)).ToList();
            Assert.Contains(l_S行, l => l.StartsWith(C_項目_S_2) && l.EndsWith(C_項目_CN_i_2));
            Assert.Contains(l_S行, l => l.StartsWith(C_項目_S_1) && l.EndsWith(C_項目_CN_i_1));
        }

        /// <summary>
        /// コピー数を渡さなければ CN タグを省くことを確かめる
        /// </summary>
        [Fact]
        public void V_出力_コピー数を渡さなければCNタグを省く()
        {
            var (l_unitig一覧, l_グラフ) = V_構築();
            var l_パス = Path.Combine(this._作業ディレクトリ, C_ファイル名_graph_nocn_gfa);
            GfaWriter.V_出力(l_パス, l_unitig一覧, l_グラフ, C_k長);
            var l_S行 = File.ReadAllLines(l_パス).Where(l => l.StartsWith(C_項目_S)).ToList();
            Assert.DoesNotContain(l_S行, l => l.Contains(C_項目_CN_i));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 分岐を持つ検証用の unitig グラフを組み立てる
        /// </summary>
        /// <returns>unitig 一覧とグラフ</returns>
        private static (List<string> A_unitig一覧, UnitigGraph A_グラフ) V_構築()
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
                A_スレッド数 = 1
            };
            List<string> l_unitig一覧 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_unitigID, int A_位置)> l_kmer辞書 = [];
            var l_ID = 1;
            foreach (var l_配列 in new[]
            {
                C_unitigA,
                C_unitigB,
                C_unitigC
            }

            )
            {
                l_unitig一覧.Add(l_配列);
                l_unitig一覧.Add(Util.V_逆相補(l_配列));
                for (var i = C_k長; i <= l_配列.Length; i++)
                {
                    var l_開始位置 = i - C_k長;
                    var l_キー = new KmerKey(l_配列.AsSpan(l_開始位置, C_k長));
                    V_登録(l_kmer辞書, l_キー, l_ID, l_開始位置);
                    V_登録(l_kmer辞書, l_キー.Get_逆相補(), -l_ID, l_配列.Length - i);
                }

                l_ID++;
            }

            return (l_unitig一覧, UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号));
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
