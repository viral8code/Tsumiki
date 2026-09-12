using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 短い反復配列の解きほぐし (repeat resolution) の検証
    /// </summary>
    /// <remarks>
    /// 反復配列 R がゲノム中に 2 回現れ、それぞれ A→R→C と B→R→D という文脈を持つ場合、de Bruijn グラフ上では R は 1 個の頂点に潰れて入次数 2 ・出次数 2 になる<br/>
    /// R の内部から読まれたリードはどちらのコピー由来か区別できないため、分岐でのリード支持は原理的に 5 割前後にしかならず解けない<br/>
    /// R を丸ごと跨いだフラグメントだけが手がかりになる<br/>
    /// 実データ (k=63) ではこの形の unitig が 151 本あり、うち 143 本がフラグメント長の中央値 (245 bp) より短かった
    /// </remarks>
    public class RepeatResolutionTests
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

        // k=8 で 5 本すべてを通じて重複する正規化 k-mer が無いことを確認済みの構成
        // A と B はどちらも R の先頭 k-1 塩基で終わり、C と D はどちらも
        // R の末尾 k-1 塩基で始まる (= R が入次数 2 ・出次数 2 の反復になる)

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
        /// 跨いだペアが片方の対応付けだけを支持するとき、反復を 2 本の経路に解決する
        /// </summary>
        [Fact]
        public void V_跨いだペアが単一の対応付けを支持するとき反復を2本の経路に解決する()
        {
            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(入口ユニティグ, 代替入口ユニティグ, 反復ユニティグ, 出口ユニティグ, 代替出口ユニティグ);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);

            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_反復配列 = ContigMaker.Get_頂点番号(3);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            var l_終端配列 = ContigMaker.Get_頂点番号(5);

            // 前提: R が入次数 2 ・出次数 2 の反復として構築されている
            Assert.Equal(2, l_グラフ.A_出辺[l_反復配列].Count);
            Assert.Equal(2, l_グラフ.Get_入次数(l_反復配列));

            // A-C と B-D を跨いだペアだけが観測された、という証拠を与える
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_末尾配列)] = 30UL,
                [(l_中間配列, l_終端配列)] = 28UL,
            };
            Dictionary<(int, int), ulong> l_支持 = [];
            var l_変更前頂点数 = l_グラフ.A_出辺.Count;

            var l_解決数 = l_グラフ.V_解決_短い反復(l_ユニティグ一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);

            Assert.Equal(1, l_解決数);

            // 反復が複製され、頂点が 2 つ (順鎖・逆鎖) 増えているはず
            Assert.Equal(l_変更前頂点数 + 2, l_グラフ.A_出辺.Count);
            Assert.Equal(反復ユニティグ, l_ユニティグ一覧[l_変更前頂点数]);

            // どちらの経路も「反復のコピーを 1 つだけ通る一本道」になっていること
            // 元の頂点 r と複製のどちらが A 側に残るかは辺の格納順に依存する
            // (unitig の番号とは無関係) ので、頂点の同一性ではなく
            // 経路の構造と対応付けを検証する
            V_検証_もつれ解消(l_グラフ, l_ユニティグ一覧, p_開始: l_先頭配列, p_終了: l_末尾配列, p_別開始: l_中間配列, p_別終了: l_終端配列);
        }

        /// <summary>
        /// 跨いだペアが交差した対応付けを示すときは、その通りに反復を解決する
        /// </summary>
        [Fact]
        public void V_跨いだペアが交差した対応付けを示すときはその通りに反復を解決する()
        {
            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(入口ユニティグ, 代替入口ユニティグ, 反復ユニティグ, 出口ユニティグ, 代替出口ユニティグ);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);

            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_中間配列 = ContigMaker.Get_頂点番号(2);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);
            var l_終端配列 = ContigMaker.Get_頂点番号(5);

            // 今度は A-D と B-C の組み合わせが支持されている
            Dictionary<(int, int), ulong> l_ペア連結 = new()
            {
                [(l_先頭配列, l_終端配列)] = 25UL,
                [(l_中間配列, l_末尾配列)] = 31UL,
            };
            Dictionary<(int, int), ulong> l_支持 = [];

            var l_解決数 = l_グラフ.V_解決_短い反復(l_ユニティグ一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);

            Assert.Equal(1, l_解決数);
            // 交差した対応付け: A は D へ、B は C へ繋がる
            V_検証_もつれ解消(l_グラフ, l_ユニティグ一覧, p_開始: l_先頭配列, p_終了: l_終端配列, p_別開始: l_中間配列, p_別終了: l_末尾配列);
        }

        /// <summary>
        /// 両方の対応付けが同程度に支持されている場合、どちらが正しいか判断できない
        /// </summary>
        /// <remarks>
        /// 誤った繋ぎ方は誤アセンブリを生むため、繋がずに残すのが正しい
        /// </remarks>
        [Fact]
        public void V_ペアがどちらの対応付けも優勢に支持しないときは反復をそのまま残す()
        {
            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(入口ユニティグ, 代替入口ユニティグ, 反復ユニティグ, 出口ユニティグ, 代替出口ユニティグ);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);

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

            var l_解決数 = l_グラフ.V_解決_短い反復(l_ユニティグ一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);

            Assert.Equal(0, l_解決数);
            Assert.Equal(l_変更前頂点数, l_グラフ.A_出辺.Count);
            Assert.Equal(2, l_グラフ.A_出辺[l_反復配列].Count);
            Assert.Equal(2, l_グラフ.Get_入次数(l_反復配列));
        }

        /// <summary>
        /// フラグメントで跨げない長さの反復は、そもそも証拠が得られないので対象外
        /// </summary>
        /// <remarks>
        /// (跨げていないのに偶然の対応付けで繋ぐと誤アセンブリになる)
        /// </remarks>
        [Fact]
        public void V_フラグメントが跨げない長さの反復は解決を見送る()
        {
            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(入口ユニティグ, 代替入口ユニティグ, 反復ユニティグ, 出口ユニティグ, 代替出口ユニティグ);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);

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

            // R は 24 bp なので、上限を 10 bp にすれば対象外になる
            var l_解決数 = l_グラフ.V_解決_短い反復(l_ユニティグ一覧, l_支持, l_ペア連結, p_反復長の上限: 10, p_優勢閾値: 0.8M, p_最小証拠数: 5UL);

            Assert.Equal(0, l_解決数);
            Assert.Equal(l_変更前頂点数, l_グラフ.A_出辺.Count);
        }

        /// <summary>
        /// 跨いだペアが少なすぎる場合も、偶然の一致で繋いでしまわないよう見送る
        /// </summary>
        /// <summary>
        /// 跨いだペア数が最小証拠数に満たない場合は、解決を見送る
        /// </summary>
        [Fact]
        public void V_跨いだペア数が最小証拠数に満たない場合は解決を見送る()
        {
            var (l_ユニティグ一覧, l_kmer辞書) = V_構築(入口ユニティグ, 代替入口ユニティグ, 反復ユニティグ, 出口ユニティグ, 代替出口ユニティグ);
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ一覧, l_kmer辞書, k長, 曖昧kmer番号);

            var l_先頭配列 = ContigMaker.Get_頂点番号(1);
            var l_末尾配列 = ContigMaker.Get_頂点番号(4);

            Dictionary<(int, int), ulong> l_ペア連結 = new() { [(l_先頭配列, l_末尾配列)] = 2UL };
            Dictionary<(int, int), ulong> l_支持 = [];
            var l_変更前頂点数 = l_グラフ.A_出辺.Count;

            var l_解決数 = l_グラフ.V_解決_短い反復(l_ユニティグ一覧, l_支持, l_ペア連結, p_反復長の上限: 500, p_優勢閾値: 0.8M, p_最小証拠数: 10UL);

            Assert.Equal(0, l_解決数);
            Assert.Equal(l_変更前頂点数, l_グラフ.A_出辺.Count);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 与えたユニティグ配列群から、グラフ構築に使うユニティグ一覧と kmer 辞書を組み立てる
        /// </summary>
        /// <param name="p_ユニティグ配列">構築元にするユニティグの配列</param>
        /// <returns>ユニティグ一覧と kmer 辞書の組</returns>
        private static (List<string> A_ユニティグ一覧, Dictionary<KmerKey, (int A_ユニティグID, int A_位置)> A_kmer辞書) V_構築(params string[] p_ユニティグ配列)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k長, A_スレッド数 = 1 };
            List<string> l_ユニティグ一覧 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_ユニティグID, int A_位置)> l_kmer辞書 = [];

            var l_ID = 1;
            foreach (var l_配列 in p_ユニティグ配列)
            {
                l_ユニティグ一覧.Add(l_配列);
                l_ユニティグ一覧.Add(Util.V_逆相補(l_配列));
                for (var i = k長; i <= l_配列.Length; i++)
                {
                    var l_開始位置 = i - k長;
                    var l_kmer = new KmerKey(l_配列.AsSpan(l_開始位置, k長));
                    V_登録(l_kmer辞書, l_kmer, l_ID, l_開始位置);
                    V_登録(l_kmer辞書, l_kmer.Get_逆相補(), -l_ID, l_配列.Length - i);
                }
                l_ID++;
            }
            return (l_ユニティグ一覧, l_kmer辞書);
        }

        /// <summary>
        /// k-mer を、それが載るユニティグと開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書">登録先の辞書</param>
        /// <param name="p_kmer">登録する k-mer</param>
        /// <param name="p_ID">ユニティグ ID</param>
        /// <param name="p_位置">ユニティグ内の開始位置</param>
        private static void V_登録(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_kmer, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_kmer, out var l_既存))
            {
                if (l_既存.Item1 is 曖昧kmer番号 || l_既存.Item1 == p_ID)
                {
                    return;
                }
                p_辞書[p_kmer] = (曖昧kmer番号, 0);
                return;
            }
            p_辞書[p_kmer] = (p_ID, p_位置);
        }

        /// <summary>
        /// from → (反復のコピー) → to と otherFrom → (別のコピー) → otherTo という 2 本の独立した一本道になっていることを検証する
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_ユニティグ一覧"></param>
        /// <param name="p_開始"></param>
        /// <param name="p_終了"></param>
        /// <param name="p_別開始"></param>
        /// <param name="p_別終了"></param>
        private static void V_検証_もつれ解消(UnitigGraph p_グラフ, List<string> p_ユニティグ一覧, int p_開始, int p_終了, int p_別開始, int p_別終了)
        {
            var l_経路1 = Assert.Single(p_グラフ.A_出辺[p_開始]);
            var l_経路2 = Assert.Single(p_グラフ.A_出辺[p_別開始]);

            // それぞれ別のコピーを通ること (同じ頂点を共有していたら解けていない)
            Assert.NotEqual(l_経路1, l_経路2);

            // 通る頂点はどちらも反復配列そのもの
            Assert.Equal(反復ユニティグ, p_ユニティグ一覧[l_経路1]);
            Assert.Equal(反復ユニティグ, p_ユニティグ一覧[l_経路2]);

            // 各コピーは入次数 1 ・出次数 1 の一本道
            Assert.Equal(1, p_グラフ.Get_入次数(l_経路1));
            Assert.Equal(1, p_グラフ.Get_入次数(l_経路2));
            Assert.Equal([p_終了], p_グラフ.A_出辺[l_経路1]);
            Assert.Equal([p_別終了], p_グラフ.A_出辺[l_経路2]);

            // 逆鎖側も対称であること (片側だけ付け替えるとグラフが壊れ、
            // 順鎖と逆鎖で別々の経路が組まれてしまう)
            Assert.Contains(l_経路1 ^ 1, p_グラフ.A_出辺[p_終了 ^ 1]);
            Assert.Contains(l_経路2 ^ 1, p_グラフ.A_出辺[p_別終了 ^ 1]);
            Assert.Contains(p_開始 ^ 1, p_グラフ.A_出辺[l_経路1 ^ 1]);
            Assert.Contains(p_別開始 ^ 1, p_グラフ.A_出辺[l_経路2 ^ 1]);
        }

        #endregion

    }
}
