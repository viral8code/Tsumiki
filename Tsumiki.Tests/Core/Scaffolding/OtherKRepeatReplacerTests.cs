using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終配列のうち、ほかの k の contig と錨の間の長さが食い違う所を、リードに続けて出てくる形に差し替える処理 (OtherKRepeatReplacer) の検証<br/>
    /// G = A (1,000 bp) + 単位 100 bp × 3 + B (1,000 bp)、リードは G の 250 bp の窓 (5 bp おき) を 2 回ずつ、リード長 250 (窓 175 bp)。<br/>
    /// 窓で比べる場面は、単位 12 bp × 17 の G で、ReadWindowSelector の窓は 212 bp
    /// </summary>
    public class OtherKRepeatReplacerTests
    {
        #region 定数

        /// <summary>
        /// 縦の反復の両脇の長さ
        /// </summary>
        private const int C_脇の長さ = 1000;

        /// <summary>
        /// 縦の反復の単位の長さ
        /// </summary>
        private const int C_単位長 = 100;

        /// <summary>
        /// G の縦の反復の単位の数
        /// </summary>
        private const int C_単位の数 = 3;

        /// <summary>
        /// ゲノムの乱数の種
        /// </summary>
        private const int C_ゲノムの種 = 7373;

        /// <summary>
        /// 代表リード長
        /// </summary>
        private const int C_リード長 = 250;

        /// <summary>
        /// リードの窓をずらす幅
        /// </summary>
        private const int C_リードのずれ = 5;

        /// <summary>
        /// ほかの k の contig に入れる、両脇の長さ
        /// </summary>
        private const int C_ほかのkの脇 = 500;

        /// <summary>
        /// 配列 ID S1
        /// </summary>
        private const string C_ID_S1 = "S1";

        /// <summary>
        /// 配列 ID S2
        /// </summary>
        private const string C_ID_S2 = "S2";

        /// <summary>
        /// 場面 縮んだ反復
        /// </summary>
        private const string C_場面_縮んだ反復 = "縮んだ反復";

        /// <summary>
        /// 場面 逆相補
        /// </summary>
        private const string C_場面_逆相補 = "逆相補";

        /// <summary>
        /// 場面 骨格が正しい
        /// </summary>
        private const string C_場面_骨格が正しい = "骨格が正しい";

        /// <summary>
        /// 場面 リード無し
        /// </summary>
        private const string C_場面_リード無し = "リード無し";

        /// <summary>
        /// 場面 左の片が短い
        /// </summary>
        private const string C_場面_片が短い = "片が短い";

        /// <summary>
        /// 場面 ほかの k が 2 通り
        /// </summary>
        private const string C_場面_2通り = "2 通り";

        /// <summary>
        /// 場面 最終配列が正しいがリードの欠けで通らず、ほかの k は単位が 1 つ足りない
        /// </summary>
        private const string C_場面_リードの欠け = "リードの欠け";

        /// <summary>
        /// リードの欠けの場面で、含むリードを除く G の窓の始まり
        /// </summary>
        private const int C_欠けの始め = 1050;

        /// <summary>
        /// 窓で比べる場面の、縦の反復の単位の長さ
        /// </summary>
        private const int C_短い単位長 = 12;

        /// <summary>
        /// 窓で比べる場面の、G の縦の反復の単位の数 (204 bp)
        /// </summary>
        private const int C_短い単位の数 = 17;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 差し替えるか: 単位 1 つに縮んだ最終配列 (A + 単位 + B) は、単位 3 つのほかの k の contig (逆相補でも) の形に差し替えて G にする<br/>
        /// 最終配列が G でほかの k が縮んだ形のとき、リードが無いとき、左の片が窓より短いとき、ほかの k が単位 2 つと 3 つの 2 通りでどちらもリードに続けて出るときは残し、理由ごとの数を数える<br/>
        /// 最終配列が G で、G[1050..1225) を含むリードを除いて最終配列の形が通らず、ほかの k が単位 2 つの形のときは、ほかの k の形だけが窓では続けて出るが、<br/>
        /// 左の錨から短いほうの区間と錨を合わせた長さ (両端の錨を含む) の配列がリードに無いので残す
        /// </summary>
        /// <param name="p_場面">場面</param>
        /// <param name="p_Is差し替える">ほかの k の形に差し替えるはずか</param>
        /// <param name="p_片が短い数">期待する片が短い数</param>
        /// <param name="p_骨格がリードにある数">期待する骨格がリードにある数</param>
        /// <param name="p_決まらない数">期待する決まらない数</param>
        /// <param name="p_錨を確かめられない数">期待する錨を確かめられない数</param>
        [Theory]
        [InlineData(C_場面_縮んだ反復, true, 0, 0, 0, 0)]
        [InlineData(C_場面_逆相補, true, 0, 0, 0, 0)]
        [InlineData(C_場面_骨格が正しい, false, 0, 1, 0, 0)]
        [InlineData(C_場面_リード無し, false, 0, 0, 1, 0)]
        [InlineData(C_場面_片が短い, false, 1, 0, 0, 0)]
        [InlineData(C_場面_2通り, false, 0, 0, 1, 0)]
        [InlineData(C_場面_リードの欠け, false, 0, 0, 0, 1)]
        public void 差し替えるか_最終配列の形がリードに無くほかのkの形がただ1つ続けて出るときだけ差し替える(string p_場面, bool p_Is差し替える, int p_片が短い数, int p_骨格がリードにある数, int p_決まらない数, int p_錨を確かめられない数)
        {
            var (l_A, l_単位, l_B) = V_生成_縦の反復();
            var l_G = l_A + V_繰り返す(l_単位, C_単位の数) + l_B;
            var l_ほかのkの正しい形 = l_A[C_ほかのkの脇..] + V_繰り返す(l_単位, C_単位の数) + l_B[..C_ほかのkの脇];
            var l_元の配列 = p_場面 switch
            {
                C_場面_骨格が正しい or C_場面_リードの欠け => l_G,
                C_場面_片が短い => l_A[(C_脇の長さ - 150)..] + l_単位 + l_B,
                _ => l_A + l_単位 + l_B,
            };
            List<string> l_ほかのk群 = p_場面 switch
            {
                C_場面_逆相補 => [Util.V_逆相補(l_ほかのkの正しい形)],
                C_場面_骨格が正しい => [l_A[C_ほかのkの脇..] + l_単位 + l_B[..C_ほかのkの脇]],
                C_場面_2通り => [l_ほかのkの正しい形, l_A[C_ほかのkの脇..] + V_繰り返す(l_単位, 2) + l_B[..C_ほかのkの脇]],
                C_場面_リードの欠け => [l_A[C_ほかのkの脇..] + V_繰り返す(l_単位, 2) + l_B[..C_ほかのkの脇]],
                _ => [l_ほかのkの正しい形],
            };
            var l_欠け = l_G.Substring(C_欠けの始め, C_リード長 * 7 / 10);
            List<string> l_リード = p_場面 switch
            {
                C_場面_リード無し => [],
                C_場面_リードの欠け => [.. V_リードの窓(l_G).Where(x => !x.Contains(l_欠け, StringComparison.Ordinal))],
                _ => V_リードの窓(l_G),
            };
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKRepeatReplacer.Get_差し替えた配列群([(C_ID_S1, l_元の配列)], l_ほかのk群, l_索引, C_リード長, out var l_統計);

            Assert.Equal(C_ID_S1, Assert.Single(l_結果).A_ID);
            Assert.Equal(p_Is差し替える ? l_G : l_元の配列, l_結果[0].A_配列);
            Assert.Equal((1, p_Is差し替える ? 1 : 0, 0, p_片が短い数, p_骨格がリードにある数, p_決まらない数, p_錨を確かめられない数), (l_統計.A_対象数, l_統計.A_差し替えた数, l_統計.A_窓で選んだ数, l_統計.A_片が短い数, l_統計.A_骨格がリードにある数, l_統計.A_決まらない数, l_統計.A_錨を確かめられない数));
        }

        /// <summary>
        /// 窓で比べる: 単位 12 bp の縦の反復で G は単位 17 個。最終配列が 16 個 (192 bp) なら 175 bp の窓ではリードに続けて出てくるが、<br/>
        /// 212 bp の窓で比べると最終配列の形にだけある窓がリードに無いので、17 個のほかの k の形に差し替えて G にする。<br/>
        /// 最終配列が 17 個でほかの k が 18 個なら、ほかの k の形にだけある窓がリードに無いので残す
        /// </summary>
        /// <param name="p_最終配列の単位の数">最終配列の縦の反復の単位の数</param>
        /// <param name="p_ほかのkの単位の数">ほかの k の contig の縦の反復の単位の数</param>
        /// <param name="p_Is差し替える">ほかの k の形に差し替えるはずか</param>
        [Theory]
        [InlineData(16, 17, true)]
        [InlineData(17, 18, false)]
        public void 窓で比べる_最終配列の形もリードに続けて出るときは窓で選んだほかのkの形だけに差し替える(int p_最終配列の単位の数, int p_ほかのkの単位の数, bool p_Is差し替える)
        {
            var l_A = V_生成_ゲノム(C_脇の長さ, C_ゲノムの種);
            var l_単位 = V_生成_ゲノム(C_短い単位長, C_ゲノムの種 + 3);
            var l_B = V_生成_ゲノム(C_脇の長さ, C_ゲノムの種 + 2);
            var l_G = l_A + V_繰り返す(l_単位, C_短い単位の数) + l_B;
            var l_元の配列 = l_A + V_繰り返す(l_単位, p_最終配列の単位の数) + l_B;
            var l_リード = V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKRepeatReplacer.Get_差し替えた配列群([(C_ID_S1, l_元の配列)], [l_A[C_ほかのkの脇..] + V_繰り返す(l_単位, p_ほかのkの単位の数) + l_B[..C_ほかのkの脇]], l_索引, C_リード長, out var l_統計);

            Assert.Equal(p_Is差し替える ? l_G : l_元の配列, l_結果[0].A_配列);
            Assert.Equal((1, p_Is差し替える ? 1 : 0, p_Is差し替える ? 1 : 0, p_Is差し替える ? 0 : 1), (l_統計.A_対象数, l_統計.A_差し替えた数, l_統計.A_窓で選んだ数, l_統計.A_骨格がリードにある数));
        }

        /// <summary>
        /// 確かめの限界: 最終配列が A + 単位 + B で、ほかの k が単位 4 つ (G は単位 3 つ) なら、単位 4 つの形に差し替える<br/>
        /// 単位 3 つの反復 (300 bp) はリードより長く、単位 4 つの形の窓もどれも G にあるので、リードでは単位の数を確かめられず、ほかの k の形に任せる
        /// </summary>
        [Fact]
        public void 確かめの限界_リードより長い縦の反復の単位の数はほかのkの形に任せる()
        {
            var (l_A, l_単位, l_B) = V_生成_縦の反復();
            var l_リード = V_リードの窓(l_A + V_繰り返す(l_単位, C_単位の数) + l_B);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKRepeatReplacer.Get_差し替えた配列群([(C_ID_S1, l_A + l_単位 + l_B)], [l_A[C_ほかのkの脇..] + V_繰り返す(l_単位, 4) + l_B[..C_ほかのkの脇]], l_索引, C_リード長, out var l_統計);

            Assert.Equal(l_A + V_繰り返す(l_単位, 4) + l_B, l_結果[0].A_配列);
            Assert.Equal((1, 1), (l_統計.A_対象数, l_統計.A_差し替えた数));
        }

        /// <summary>
        /// 錨の間の確かめ: 確かめる長さが窓長以下なら真、リード長より長ければ偽、その間なら左の錨からその長さの配列がリードに 2 か所以上あるときだけ真<br/>
        /// G の窓 (5 bp おき) を 2 回ずつのリードで、G[1000..1200) は真、G[1000..1101) + G[1102..1201) (1 塩基抜き) は偽
        /// </summary>
        /// <param name="p_確かめる長さ">確かめる配列の長さ</param>
        /// <param name="p_Is1塩基抜き">確かめる配列から 1 塩基を抜くか</param>
        /// <param name="p_期待">期待する結果</param>
        [Theory]
        [InlineData(175, true, true)]
        [InlineData(200, false, true)]
        [InlineData(200, true, false)]
        [InlineData(251, false, false)]
        public void 錨の間の確かめ_窓より長くリード長以下ならリードに2か所以上あるときだけ真(int p_確かめる長さ, bool p_Is1塩基抜き, bool p_期待)
        {
            var (l_A, l_単位, l_B) = V_生成_縦の反復();
            var l_G = l_A + V_繰り返す(l_単位, C_単位の数) + l_B;
            var l_リード = V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());
            var l_配列 = p_Is1塩基抜き ? l_G[..1101] + l_G[1102..] : l_G;

            Assert.Equal(p_期待, OtherKRepeatReplacer.Is錨の間がリードで確かめられる(l_配列, C_脇の長さ, p_確かめる長さ, l_索引, C_リード長, C_リード長 * 7 / 10));
        }

        /// <summary>
        /// 近い 2 つの食い違い: G = A + P (6 bp) + M (60 bp) + Q (6 bp) + B で、最終配列が A + M + B なら、ほかの k の G の形に差し替える<br/>
        /// 1 つずつでは長さの差 (6 bp) が下限 (11 bp) に届かないが、間の錨の並び (M) が短いので 1 つの区間 (差 12 bp) にまとめる
        /// </summary>
        [Fact]
        public void 近い2つの食い違い_間が短ければ1つの区間にまとめて差し替える()
        {
            var l_A = V_生成_ゲノム(C_脇の長さ, C_ゲノムの種);
            var l_P = V_生成_ゲノム(6, C_ゲノムの種 + 4);
            var l_M = V_生成_ゲノム(60, C_ゲノムの種 + 5);
            var l_Q = V_生成_ゲノム(6, C_ゲノムの種 + 6);
            var l_B = V_生成_ゲノム(C_脇の長さ, C_ゲノムの種 + 7);
            var l_G = l_A + l_P + l_M + l_Q + l_B;
            var l_リード = V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKRepeatReplacer.Get_差し替えた配列群([(C_ID_S1, l_A + l_M + l_B)], [l_G[C_ほかのkの脇..^C_ほかのkの脇]], l_索引, C_リード長, out var l_統計);

            Assert.Equal(l_G, l_結果[0].A_配列);
            Assert.Equal((1, 1), (l_統計.A_対象数, l_統計.A_差し替えた数));
        }

        /// <summary>
        /// 別の配列に当たる錨: G = A + 単位 + 3 塩基違う単位 + 単位 + B で、最終配列が A + 単位 + B と、3 塩基違う単位を含む別の短い配列 (S2) なら、<br/>
        /// ほかの k の contig (G の形) の上で、3 塩基違う単位の錨は S2 に当たるが、それを飛ばして S1 の錨の間で食い違いを見つけ、G に差し替える。S2 はそのまま
        /// </summary>
        [Fact]
        public void 別の配列に当たる錨_飛ばして同じ配列の錨の間で差し替える()
        {
            var (l_A, l_単位, l_B) = V_生成_縦の反復();
            var l_違う単位 = string.Concat(l_単位.Select((x, i) => i is 20 or 50 or 80 ? (x == Consts.塩基コード.A ? Consts.塩基コード.C : Consts.塩基コード.A) : x));
            var l_G = l_A + l_単位 + l_違う単位 + l_単位 + l_B;
            var l_S2 = V_生成_ゲノム(100, C_ゲノムの種 + 11) + l_違う単位 + V_生成_ゲノム(100, C_ゲノムの種 + 12);
            var l_リード = V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKRepeatReplacer.Get_差し替えた配列群([(C_ID_S1, l_A + l_単位 + l_B), (C_ID_S2, l_S2)], [l_G[C_ほかのkの脇..^C_ほかのkの脇]], l_索引, C_リード長, out var l_統計);

            Assert.Equal([l_G, l_S2], l_結果.Select(static x => x.A_配列));
            Assert.Equal((1, 1), (l_統計.A_対象数, l_統計.A_差し替えた数));
        }

        /// <summary>
        /// ほかの k の側の N: ほかの k の配列 (単位 3 つの形) の 2 つ目の単位に N があれば、食い違いの区間にせず、最終配列 (A + 単位 + B) はそのまま残す
        /// </summary>
        [Fact]
        public void ほかのkの側のN_区間にしない()
        {
            var (l_A, l_単位, l_B) = V_生成_縦の反復();
            var l_リード = V_リードの窓(l_A + V_繰り返す(l_単位, C_単位の数) + l_B);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());
            var l_元の配列 = l_A + l_単位 + l_B;

            var l_結果 = OtherKRepeatReplacer.Get_差し替えた配列群([(C_ID_S1, l_元の配列)], [l_A[C_ほかのkの脇..] + l_単位 + l_単位[..50] + Consts.塩基コード.N + l_単位[51..] + l_単位 + l_B[..C_ほかのkの脇]], l_索引, C_リード長, out var l_統計);

            Assert.Equal(l_元の配列, l_結果[0].A_配列);
            Assert.Equal((0, 0), (l_統計.A_対象数, l_統計.A_差し替えた数));
        }

        /// <summary>
        /// 小さい食い違い: 最終配列が G から B の中の 5 bp を除いた形なら、長さの差が下限 (11 bp) に届かないので区間にせず、そのまま残す
        /// </summary>
        [Fact]
        public void 小さい食い違い_長さの差が下限に届かなければ対象にしない()
        {
            var (l_A, l_単位, l_B) = V_生成_縦の反復();
            var l_G = l_A + V_繰り返す(l_単位, C_単位の数) + l_B;
            var l_元の配列 = l_G[..1600] + l_G[1605..];
            var l_リード = V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKRepeatReplacer.Get_差し替えた配列群([(C_ID_S1, l_元の配列)], [l_G[C_ほかのkの脇..^C_ほかのkの脇]], l_索引, C_リード長, out var l_統計);

            Assert.Equal(l_元の配列, l_結果[0].A_配列);
            Assert.Equal((0, 0), (l_統計.A_対象数, l_統計.A_差し替えた数));
        }

        /// <summary>
        /// 錨の表: 1 回だけの 31-mer は場所を持ち、逆相補で 2 回目が出る 31-mer は一意でない印 (-1) になる
        /// </summary>
        [Fact]
        public void 錨の表_逆相補で2回目が出るkmerは一意でない()
        {
            var l_X = V_生成_ゲノム(40, C_ゲノムの種 + 8);
            var l_Y = V_生成_ゲノム(40, C_ゲノムの種 + 9);

            var l_表 = OtherKRepeatReplacer.Get_一意な錨の表([(C_ID_S1, l_X), (C_ID_S2, l_Y + Util.V_逆相補(l_X[..31]))]);

            var (_, l_Xの頭, _) = OtherKRepeatReplacer.Get_錨群(l_X).First();
            var (_, l_Yの頭, l_Yの頭の順向き) = OtherKRepeatReplacer.Get_錨群(l_Y).First();
            Assert.Equal(-1L, l_表[l_Xの頭]);
            Assert.Equal((1L << 32) | (l_Yの頭の順向き ? 1L : 0L), l_表[l_Yの頭]);
        }

        /// <summary>
        /// 錨の列挙: 長さ 40 の配列から 10 個の 31-mer を列挙し、逆相補の配列からは同じ正準形が逆の順と逆の向きで出る。N を含む 31-mer は除く
        /// </summary>
        [Fact]
        public void 錨の列挙_逆相補では同じ正準形が逆の順と逆の向きで出る()
        {
            var l_X = V_生成_ゲノム(40, C_ゲノムの種 + 10);

            var l_順 = OtherKRepeatReplacer.Get_錨群(l_X).ToList();
            var l_逆 = OtherKRepeatReplacer.Get_錨群(Util.V_逆相補(l_X)).ToList();

            Assert.Equal(10, l_順.Count);
            Assert.Equal(l_順.Select(static x => (x.A_正準, !x.A_順向き)), l_逆.AsEnumerable().Reverse().Select(static x => (x.A_正準, x.A_順向き)));
            Assert.Equal(5, OtherKRepeatReplacer.Get_錨群(l_X[..20] + Consts.塩基コード.N + l_X[..35]).Count());
        }

        /// <summary>
        /// 区間群: 最終配列の上で重なる食い違いは 1 つの区間にまとめ、触れるだけ (終わり = 次の始め) の食い違いは分ける
        /// </summary>
        [Fact]
        public void 区間群_重なる食い違いだけをまとめる()
        {
            var l_区間群 = OtherKRepeatReplacer.Get_区間群([(0, 50, 80, "A"), (0, 10, 60, "C"), (0, 80, 90, "G")]);

            Assert.Equal([(10, 80, 2), (80, 90, 1)], l_区間群.Select(static x => (x.A_始め, x.A_終わり, x.A_形群.Count)));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 縦の反復を挟む両脇 A・B と、反復の単位を作る
        /// </summary>
        /// <returns>A、単位、B</returns>
        private static (string A_A, string A_単位, string A_B) V_生成_縦の反復()
        {
            return (V_生成_ゲノム(C_脇の長さ, C_ゲノムの種), V_生成_ゲノム(C_単位長, C_ゲノムの種 + 1), V_生成_ゲノム(C_脇の長さ, C_ゲノムの種 + 2));
        }

        /// <summary>
        /// 配列を指定の回数だけ繋ぐ
        /// </summary>
        /// <param name="p_配列">繋ぐ配列</param>
        /// <param name="p_回数">回数</param>
        /// <returns>繋いだ配列</returns>
        private static string V_繰り返す(string p_配列, int p_回数)
        {
            return string.Concat(Enumerable.Repeat(p_配列, p_回数));
        }

        /// <summary>
        /// 乱数の種を固定したゲノム (ACGT のみ) を作る
        /// </summary>
        /// <param name="p_長さ">長さ</param>
        /// <param name="p_種">乱数の種</param>
        /// <returns>ゲノム</returns>
        private static string V_生成_ゲノム(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[l_乱数.Next(4)]));
        }

        /// <summary>
        /// 配列の C_リード長 の窓を C_リードのずれ おきに切り出す<br/>
        /// 末尾まで届く窓も含める
        /// </summary>
        /// <param name="p_配列">切り出す元の配列</param>
        /// <returns>リードの一覧</returns>
        private static List<string> V_リードの窓(string p_配列)
        {
            return [.. Enumerable.Range(0, ((p_配列.Length - C_リード長) / C_リードのずれ) + 1).Select(i => p_配列.Substring(i * C_リードのずれ, C_リード長))];
        }

        #endregion
    }
}
