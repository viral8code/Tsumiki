using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 短い反復配列の解きほぐし (repeat resolution) の検証
    /// </summary>
    public class RepeatResolutionTests
    {
        #region 定数

        /// <summary>
        /// 曖昧塩基を含む窓を表す番号
        /// </summary>
        private const int C_曖昧kmer番号 = int.MinValue;

        /// <summary>
        /// この検証で使う k 長
        /// </summary>
        private const int C_k長 = 8;

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
        /// 跨いだペアが片方の対応付けだけを支持するとき、反復を 2 本の経路に解決する
        /// </summary>
        [Fact]
        public void V_跨いだペアが単一の対応付けを支持するとき反復を2本の経路に解決する()
        {
            var (l_unitig一覧, l_kmer辞書) = V_構築(C_入口unitig, C_代替入口unitig, C_反復unitig, C_出口unitig, C_代替出口unitig);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号);
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_反復配列 = ContigMaker.Get_頂点番号(3);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            var l_終端配列 = ContigMaker.Get_頂点番号(5);
            Assert.Equal(2, l_グラフ.A_出辺[l_反復配列].Count);
            Assert.Equal(2, l_グラフ.Get_入次数(l_反復配列));
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_末尾配列)] = 30UL,
                [(l_中間配列, l_終端配列)] = 28UL,
            };
            Dictionary<(int, int), ulong> l_支持 = [];
            var l_変更前頂点数 = l_グラフ.A_出辺.Count;
            var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(1, l_解決数);
            Assert.Equal(l_変更前頂点数 + 2, l_グラフ.A_出辺.Count);
            Assert.Equal(C_反復unitig, l_unitig一覧[l_変更前頂点数]);
            V_検証_もつれ解消(l_グラフ, l_unitig一覧, p_開始: l_先頭配列, p_終了: l_末尾配列, p_別開始: l_中間配列, p_別終了: l_終端配列);
        }

        /// <summary>
        /// 跨いだペアが交差した対応付けを示すときは、その通りに反復を解決する
        /// </summary>
        [Fact]
        public void V_跨いだペアが交差した対応付けを示すときはその通りに反復を解決する()
        {
            var (l_unitig一覧, l_kmer辞書) = V_構築(C_入口unitig, C_代替入口unitig, C_反復unitig, C_出口unitig, C_代替出口unitig);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号);
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            var l_終端配列 = ContigMaker.Get_頂点番号(5);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 25UL,
                [(l_中間配列, l_末尾配列)] = 31UL,
            };
            Dictionary<(int, int), ulong> l_支持 = [];
            var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(1, l_解決数);
            V_検証_もつれ解消(l_グラフ, l_unitig一覧, p_開始: l_先頭配列, p_終了: l_終端配列, p_別開始: l_中間配列, p_別終了: l_末尾配列);
        }

        /// <summary>
        /// 両方の対応付けが同程度に支持されている場合、どちらが正しいか判断できない
        /// </summary>
        [Fact]
        public void V_ペアがどちらの対応付けも優勢に支持しないときは反復をそのまま残す()
        {
            var (l_unitig一覧, l_kmer辞書) = V_構築(C_入口unitig, C_代替入口unitig, C_反復unitig, C_出口unitig, C_代替出口unitig);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号);
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_反復配列 = ContigMaker.Get_頂点番号(3);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            var l_終端配列 = ContigMaker.Get_頂点番号(5);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_末尾配列)] = 15UL,
                [(l_中間配列, l_終端配列)] = 14UL,
                [(l_先頭配列, l_終端配列)] = 13UL,
                [(l_中間配列, l_末尾配列)] = 16UL,
            };
            Dictionary<(int, int), ulong> l_支持 = [];
            var l_変更前頂点数 = l_グラフ.A_出辺.Count;
            var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(0, l_解決数);
            Assert.Equal(l_変更前頂点数, l_グラフ.A_出辺.Count);
            Assert.Equal(2, l_グラフ.A_出辺[l_反復配列].Count);
            Assert.Equal(2, l_グラフ.Get_入次数(l_反復配列));
        }

        /// <summary>
        /// フラグメントで跨げない長さの反復は、そもそも証拠が得られないので対象外
        /// </summary>
        [Fact]
        public void V_フラグメントが跨げない長さの反復は解決を見送る()
        {
            var (l_unitig一覧, l_kmer辞書) = V_構築(C_入口unitig, C_代替入口unitig, C_反復unitig, C_出口unitig, C_代替出口unitig);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号);
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            var l_終端配列 = ContigMaker.Get_頂点番号(5);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_末尾配列)] = 30UL,
                [(l_中間配列, l_終端配列)] = 28UL,
            };
            Dictionary<(int, int), ulong> l_支持 = [];
            var l_変更前頂点数 = l_グラフ.A_出辺.Count;
            var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, l_支持, l_ペア連結, p_反復長の上限: 10, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);
            Assert.Equal(0, l_解決数);
            Assert.Equal(l_変更前頂点数, l_グラフ.A_出辺.Count);
        }

        /// <summary>
        /// 跨いだペアが少なすぎる場合も、偶然の一致で繋いでしまわないよう見送る
        /// </summary>
        [Fact]
        public void V_跨いだペア数が最小証拠数に満たない場合は解決を見送る()
        {
            var (l_unitig一覧, l_kmer辞書) = V_構築(C_入口unitig, C_代替入口unitig, C_反復unitig, C_出口unitig, C_代替出口unitig);
            var l_グラフ = UnitigGraph.Get_グラフ(l_unitig一覧, l_kmer辞書, C_k長, C_曖昧kmer番号);
            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_末尾配列)] = 2UL
            };
            Dictionary<(int, int), ulong> l_支持 = [];
            var l_変更前頂点数 = l_グラフ.A_出辺.Count;
            var l_解決数 = l_グラフ.V_解決_短い反復(l_unitig一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 10UL);
            Assert.Equal(0, l_解決数);
            Assert.Equal(l_変更前頂点数, l_グラフ.A_出辺.Count);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 与えた unitig 配列群から、グラフ構築に使う unitig 一覧と kmer 辞書を組み立てる
        /// </summary>
        /// <param name="p_unitig配列">構築元にする unitig の配列</param>
        /// <returns>unitig 一覧と kmer 辞書の組</returns>
        private static (List<string> A_unitig一覧, Dictionary<KmerKey, (int A_unitigID, int A_位置)> A_kmer辞書) V_構築(params string[] p_unitig配列)
        {
            ConfigurationManager.A_実行時引数 = new Parameters
            {
                A_k長 = C_k長,
                A_スレッド数 = 1
            };
            List<string> l_unitig一覧 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_unitigID, int A_位置)> l_kmer辞書 = [];
            var l_ID = 1;
            foreach (var l_配列 in p_unitig配列)
            {
                l_unitig一覧.Add(l_配列);
                l_unitig一覧.Add(Util.V_逆相補(l_配列));
                for (var i = C_k長; i <= l_配列.Length; i++)
                {
                    var l_開始位置 = i - C_k長;
                    var l_kmer = new KmerKey(l_配列.AsSpan(l_開始位置, C_k長));
                    V_登録(l_kmer辞書, l_kmer, l_ID, l_開始位置);
                    V_登録(l_kmer辞書, l_kmer.Get_逆相補(), -l_ID, l_配列.Length - i);
                }

                l_ID++;
            }

            return (l_unitig一覧, l_kmer辞書);
        }

        /// <summary>
        /// k-mer を、それが載る unitig と開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_kmer">登録する k-mer</param>
        /// <param name="p_ID">unitig ID</param>
        /// <param name="p_位置">unitig 内の開始位置</param>
        private static void V_登録(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_kmer, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_kmer, out var l_既存))
            {
                if (l_既存.Item1 is C_曖昧kmer番号 || l_既存.Item1 == p_ID)
                {
                    return;
                }

                p_辞書[p_kmer] = (C_曖昧kmer番号, 0);
                return;
            }

            p_辞書[p_kmer] = (p_ID, p_位置);
        }

        /// <summary>
        /// from → (反復のコピー) → to と otherFrom → (別のコピー) → otherTo という 2 本の独立した一本道になっていることを検証する
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_unitig一覧"></param>
        /// <param name="p_開始"></param>
        /// <param name="p_終了"></param>
        /// <param name="p_別開始"></param>
        /// <param name="p_別終了"></param>
        private static void V_検証_もつれ解消(UnitigGraph p_グラフ, List<string> p_unitig一覧, int p_開始, int p_終了, int p_別開始, int p_別終了)
        {
            var l_経路1 = Assert.Single(p_グラフ.A_出辺[p_開始]);
            var l_経路2 = Assert.Single(p_グラフ.A_出辺[p_別開始]);
            Assert.NotEqual(l_経路1, l_経路2);
            Assert.Equal(C_反復unitig, p_unitig一覧[l_経路1]);
            Assert.Equal(C_反復unitig, p_unitig一覧[l_経路2]);
            Assert.Equal(1, p_グラフ.Get_入次数(l_経路1));
            Assert.Equal(1, p_グラフ.Get_入次数(l_経路2));
            Assert.Equal([p_終了], p_グラフ.A_出辺[l_経路1]);
            Assert.Equal([p_別終了], p_グラフ.A_出辺[l_経路2]);
            Assert.Contains(l_経路1 ^ 1, p_グラフ.A_出辺[p_終了 ^ 1]);
            Assert.Contains(l_経路2 ^ 1, p_グラフ.A_出辺[p_別終了 ^ 1]);
            Assert.Contains(p_開始 ^ 1, p_グラフ.A_出辺[l_経路1 ^ 1]);
            Assert.Contains(p_別開始 ^ 1, p_グラフ.A_出辺[l_経路2 ^ 1]);
        }

        #endregion
    }
}
