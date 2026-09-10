using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 環状の複製単位に同じ反復が 2 回現れる形の解きほぐしを固定する
    /// </summary>
    /// <remarks>
    /// 環状のゲノムに反復 R が 2 回あると、その間に挟まれる領域は必ず 2 つ
    /// (A と B) になり、ゲノムは A R B R を 1 周する形になる<br/>
    /// このとき R へ
    /// 入ってくるのも A と B、R から出ていくのも A と B で、同じ unitig が
    /// 入口と出口の両方に立つ<br/>
    /// これは退化した構造ではなく、細菌のゲノムで最もありふれた反復の形<br/>
    /// 「A の次は B」なのか「A の次は A」なのかはペアエンドで区別でき、
    /// 前者なら 1 本の環、後者なら 2 本の環になる<br/>
    /// 同じ unitig が 2 つの役回りに
    /// 立つことだけを理由に触らずにいると、証拠が揃っていても解けない
    /// </remarks>
    public class CircularRepeatResolutionTests
    {
        /// <summary>
        /// 曖昧 kmer
        /// </summary>
        private const int 曖昧kmer = int.MinValue;

        /// <summary>
        /// k 長
        /// </summary>
        private const int k長 = 8;

        // 反復
        // A と B はどちらも R の先頭 k-1 塩基で終わり、
        // R の末尾 k-1 塩基で始まる (= R が入次数 2・出次数 2 になる)

        /// <summary>
        /// ユニティグ R
        /// </summary>
        private const string ユニティグR = "CTCCGTCAGCTTGTTTGGAGCAGA";

        /// <summary>
        /// ユニティグ A
        /// </summary>
        private const string ユニティグA = "GAGCAGAGTCGTTCTGCGAGGACAGTTCGCGAGCCCTCCGTC";

        /// <summary>
        /// ユニティグ B
        /// </summary>
        private const string ユニティグB = "GAGCAGACCGTCTGTAACAGCTGTATTGAGGTCGTCTCCGTC";

        /// <summary>
        /// A・B・R のユニティグ配列と、それらから作った kmer 辞書を組み立てる
        /// </summary>
        private static (List<string> A_ユニティグ配列, Dictionary<KmerKey, (int, int)> A_kmer辞書) Get_構成()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            List<string> l_ユニティグ配列 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int, int)> l_kmer辞書 = [];

            var l_ID = 1;
            foreach (var l_配列 in new[] { ユニティグA, ユニティグB, ユニティグR })
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
        private static void V_登録(
            Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存))
            {
                if (l_既存.Item1 is not 曖昧kmer && l_既存.Item1 != p_ID)
                {
                    p_辞書[p_キー] = (曖昧kmer, 0);
                }
                return;
            }
            p_辞書[p_キー] = (p_ID, p_位置);
        }

        /// <summary>
        /// A R B R の環を組み、R が入次数 2・出次数 2 になっていることまで確かめる
        /// </summary>
        private static (UnitigGraph A_グラフ, List<string> A_ユニティグ配列, int A_a, int A_b, int A_r) Get_環()
        {
            var (l_ユニティグ配列, l_kmer辞書) = Get_構成();
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ配列, l_kmer辞書, k長, 曖昧kmer);

            var l_a = ContigMaker.Get_頂点番号(1);
            var l_b = ContigMaker.Get_頂点番号(2);
            var l_r = ContigMaker.Get_頂点番号(3);

            // この前提が崩れたら以降の検証は意味を持たない
            Assert.Equal(2, l_グラフ.A_出辺[l_r].Count);
            Assert.Equal(2, l_グラフ.Get_入次数(l_r));

            // 入口も出口も A と B の 2 つ、というのがこの形の要点
            Assert.Equal(
                new HashSet<int> { l_a, l_b },
                [.. l_グラフ.A_出辺[l_r]]);

            return (l_グラフ, l_ユニティグ配列, l_a, l_b, l_r);
        }

        /// <summary>
        /// A の次は B、B の次は A というペア証拠があれば、入口と出口が同じユニティグの環でも 1 本道に解ける
        /// </summary>
        [Fact]
        public void V_解決_短い反復_入口と出口に同じユニティグが立つ環でも解ける()
        {
            var (l_グラフ, l_ユニティグ配列, l_a, l_b, _) = Get_環();
            var l_頂点数 = l_グラフ.A_出辺.Count;

            // A の次は B、B の次は A(= 1 本の環) という証拠だけを与える
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_a, l_b)] = 30,
                [(l_b, l_a)] = 28,
            };

            var l_解決数 = l_グラフ.V_解決_短い反復(
                l_ユニティグ配列, [], l_ペア連結,
                p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.Equal(1, l_解決数);
            Assert.Equal(l_頂点数 + 2, l_グラフ.A_出辺.Count);
            Assert.Equal(ユニティグR, l_ユニティグ配列[l_頂点数]);

            // A →(コピー 1)→ B →(コピー 2)→ A の一本道になっていること
            var l_Aの次 = Assert.Single(l_グラフ.A_出辺[l_a]);
            var l_Bの次 = Assert.Single(l_グラフ.A_出辺[l_b]);
            Assert.NotEqual(l_Aの次, l_Bの次);
            Assert.Equal(ユニティグR, l_ユニティグ配列[l_Aの次]);
            Assert.Equal(ユニティグR, l_ユニティグ配列[l_Bの次]);
            Assert.Equal(l_b, Assert.Single(l_グラフ.A_出辺[l_Aの次]));
            Assert.Equal(l_a, Assert.Single(l_グラフ.A_出辺[l_Bの次]));
        }

        /// <summary>
        /// A の次は A、B の次は B というペア証拠があれば、2 本の独立した環に解ける
        /// </summary>
        [Fact]
        public void V_解決_短い反復_2本の環に分かれる対応付けも同じように解ける()
        {
            var (l_グラフ, l_ユニティグ配列, l_a, l_b, _) = Get_環();

            // A の次は A、B の次は B(= 2 本の独立した環) という証拠
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_a, l_a)] = 30,
                [(l_b, l_b)] = 28,
            };

            var l_解決数 = l_グラフ.V_解決_短い反復(
                l_ユニティグ配列, [], l_ペア連結,
                p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.Equal(1, l_解決数);

            var l_Aの次 = Assert.Single(l_グラフ.A_出辺[l_a]);
            var l_Bの次 = Assert.Single(l_グラフ.A_出辺[l_b]);
            Assert.NotEqual(l_Aの次, l_Bの次);
            Assert.Equal(l_a, Assert.Single(l_グラフ.A_出辺[l_Aの次]));
            Assert.Equal(l_b, Assert.Single(l_グラフ.A_出辺[l_Bの次]));
        }

        /// <summary>
        /// 1 本の環と 2 本の環の証拠が拮抗していれば解決せずグラフを変えない
        /// </summary>
        [Fact]
        public void V_解決_短い反復_どちらの対応付けとも決まらなければ触らない()
        {
            var (l_グラフ, l_ユニティグ配列, l_a, l_b, l_r) = Get_環();

            // 1 本の環と 2 本の環がほぼ拮抗している
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_a, l_b)] = 15,
                [(l_b, l_a)] = 14,
                [(l_a, l_a)] = 13,
                [(l_b, l_b)] = 12,
            };

            var l_頂点数 = l_グラフ.A_出辺.Count;
            var l_解決数 = l_グラフ.V_解決_短い反復(
                l_ユニティグ配列, [], l_ペア連結,
                p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.Equal(0, l_解決数);
            Assert.Equal(l_頂点数, l_グラフ.A_出辺.Count);
            Assert.Equal(2, l_グラフ.A_出辺[l_r].Count);
        }

        /// <summary>
        /// 反復を跨ぐペアの証拠数が最小証拠数に届かなければ解決せずグラフを変えない
        /// </summary>
        [Fact]
        public void V_解決_短い反復_跨いだペアが足りなければ触らない()
        {
            var (l_グラフ, l_ユニティグ配列, l_a, l_b, l_r) = Get_環();

            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_a, l_b)] = 2,
                [(l_b, l_a)] = 1,
            };

            var l_解決数 = l_グラフ.V_解決_短い反復(
                l_ユニティグ配列, [], l_ペア連結,
                p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 10);

            Assert.Equal(0, l_解決数);
            Assert.Equal(2, l_グラフ.A_出辺[l_r].Count);
        }

        /// <summary>
        /// 反復の長さが上限を超えていれば解決せずグラフを変えない
        /// </summary>
        [Fact]
        public void V_解決_短い反復_反復が長すぎれば触らない()
        {
            var (l_グラフ, l_ユニティグ配列, l_a, l_b, l_r) = Get_環();

            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_a, l_b)] = 30,
                [(l_b, l_a)] = 28,
            };

            var l_解決数 = l_グラフ.V_解決_短い反復(
                l_ユニティグ配列, [], l_ペア連結,
                p_反復長の上限: ユニティグR.Length - 1, p_優勢閾値: 0.8M, p_最小証拠数: 5);

            Assert.Equal(0, l_解決数);
            Assert.Equal(2, l_グラフ.A_出辺[l_r].Count);
        }
    }
}
