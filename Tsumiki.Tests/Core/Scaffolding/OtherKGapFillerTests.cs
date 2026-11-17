using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終配列に残った N の列を、ほかの k の contig の配列でリードに確かめて埋める処理 (OtherKGapFiller) の検証<br/>
    /// 乱数の 2,000 bp の G で、アセンブリは G[..1000) + N×50 + G[1100..)、リードは G の 250 bp の窓 (5 bp おき) を 2 回ずつ
    /// </summary>
    public class OtherKGapFillerTests
    {
        #region 定数

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 2000;

        /// <summary>
        /// ゲノムの乱数の種
        /// </summary>
        private const int C_ゲノムの種 = 8181;

        /// <summary>
        /// 代表リード長
        /// </summary>
        private const int C_リード長 = 250;

        /// <summary>
        /// リードの窓をずらす幅
        /// </summary>
        private const int C_リードのずれ = 5;

        /// <summary>
        /// N の列の始まり
        /// </summary>
        private const int C_N始め = 1000;

        /// <summary>
        /// N の列の終わり (含まない)
        /// </summary>
        private const int C_N終わり = 1100;

        /// <summary>
        /// アセンブリの N の列の長さ
        /// </summary>
        private const int C_N長 = 50;

        /// <summary>
        /// ほかの k の contig の始まり
        /// </summary>
        private const int C_ほかのk始め = 500;

        /// <summary>
        /// ほかの k の contig の終わり (含まない)
        /// </summary>
        private const int C_ほかのk終わり = 1700;

        /// <summary>
        /// 1 塩基違う 2 本目の contig で、塩基を変える G の位置
        /// </summary>
        private const int C_違う塩基の位置 = 1050;

        /// <summary>
        /// 錨の長さ
        /// </summary>
        private const int C_錨長 = 60;

        /// <summary>
        /// 配列 ID S1
        /// </summary>
        private const string C_ID_S1 = "S1";

        /// <summary>
        /// 配列 ID S2
        /// </summary>
        private const string C_ID_S2 = "S2";

        /// <summary>
        /// 場面 そのまま
        /// </summary>
        private const string C_場面_そのまま = "そのまま";

        /// <summary>
        /// 場面 逆相補
        /// </summary>
        private const string C_場面_逆相補 = "逆相補";

        /// <summary>
        /// 場面 リード無し
        /// </summary>
        private const string C_場面_リード無し = "リード無し";

        /// <summary>
        /// 場面 届かない
        /// </summary>
        private const string C_場面_届かない = "届かない";

        /// <summary>
        /// 場面 錨が一意でない
        /// </summary>
        private const string C_場面_錨が一意でない = "錨が一意でない";

        /// <summary>
        /// 場面 2 通り・リードは G だけ
        /// </summary>
        private const string C_場面_2通りGだけ = "2 通り・G だけ";

        /// <summary>
        /// 場面 2 通り・リードは両方
        /// </summary>
        private const string C_場面_2通り両方 = "2 通り・両方";

        /// <summary>
        /// 場面 両脇が反復の端
        /// </summary>
        private const string C_場面_反復の端 = "反復の端";

        /// <summary>
        /// 場面 左の片が反復を調べる長さより短い
        /// </summary>
        private const string C_場面_片が短い = "片が短い";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 埋めるか: ほかの k の contig が両脇の錨を順に含み、埋めた形がリードに続けて出てくる候補がただ 1 つなら、G の形で埋める<br/>
        /// 逆相補の contig でも埋める。リードが無い、contig が右の錨に届かない、錨が一意でない、リードで 2 通りに割れるときは埋めず、理由ごとの数を数える<br/>
        /// 右の片の錨のあとに左の片の 100 bp (G[800..900)) の写しを置いて両脇を反復の端にしたときと、左の片を 400 bp (G[600..1000)) にしたときも埋めない
        /// </summary>
        /// <param name="p_場面">場面</param>
        /// <param name="p_Is埋める">G の形で埋めるはずか</param>
        /// <param name="p_錨が取れない数">期待する錨が取れない数</param>
        /// <param name="p_反復の端の数">期待する反復の端の数</param>
        /// <param name="p_ほかのkに無い数">期待するほかの k に無い数</param>
        /// <param name="p_リードに無い数">期待するリードに無い数</param>
        [Theory]
        [InlineData(C_場面_そのまま, true, 0, 0, 0, 0)]
        [InlineData(C_場面_逆相補, true, 0, 0, 0, 0)]
        [InlineData(C_場面_リード無し, false, 0, 0, 0, 1)]
        [InlineData(C_場面_届かない, false, 0, 0, 1, 0)]
        [InlineData(C_場面_錨が一意でない, false, 1, 0, 0, 0)]
        [InlineData(C_場面_2通りGだけ, true, 0, 0, 0, 0)]
        [InlineData(C_場面_2通り両方, false, 0, 0, 0, 1)]
        [InlineData(C_場面_反復の端, false, 0, 1, 0, 0)]
        [InlineData(C_場面_片が短い, false, 0, 1, 0, 0)]
        public void 埋めるか_リードに続けて出てくる候補がただ1つのときだけほかのkの配列で埋める(string p_場面, bool p_Is埋める, int p_錨が取れない数, int p_反復の端の数, int p_ほかのkに無い数, int p_リードに無い数)
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_ほかのk = l_G[C_ほかのk始め..C_ほかのk終わり];
            var l_違う形のG = l_G[..C_違う塩基の位置] + (l_G[C_違う塩基の位置] == Consts.塩基コード.A ? Consts.塩基コード.C : Consts.塩基コード.A) + l_G[(C_違う塩基の位置 + 1)..];
            var l_N = new string(Consts.塩基コード.N, C_N長);
            var l_元の配列 = p_場面 switch
            {
                C_場面_反復の端 => l_G[..C_N始め] + l_N + l_G[C_N終わり..(C_N終わり + C_錨長)] + l_G[800..900] + l_G[(C_N終わり + C_錨長 + 100)..],
                C_場面_片が短い => l_G[600..C_N始め] + l_N + l_G[C_N終わり..],
                _ => l_G[..C_N始め] + l_N + l_G[C_N終わり..],
            };
            List<(string A_ID, string A_配列)> l_配列群 = [(C_ID_S1, l_元の配列)];
            if (p_場面 == C_場面_錨が一意でない)
            {
                l_配列群.Add((C_ID_S2, l_G[(C_N始め - C_錨長)..C_N始め] + V_生成_ゲノム(300, C_ゲノムの種 + 1)));
            }

            List<string> l_ほかのk群 = p_場面 switch
            {
                C_場面_逆相補 => [Util.V_逆相補(l_ほかのk)],
                C_場面_届かない => [l_G[C_ほかのk始め..(C_N始め + 40)]],
                C_場面_2通りGだけ or C_場面_2通り両方 => [l_ほかのk, l_違う形のG[C_ほかのk始め..C_ほかのk終わり]],
                _ => [l_ほかのk],
            };
            List<string> l_リード = p_場面 switch
            {
                C_場面_リード無し => [],
                C_場面_2通り両方 => [.. V_リードの窓(l_G), .. V_リードの窓(l_違う形のG)],
                _ => V_リードの窓(l_G),
            };
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKGapFiller.Get_埋めた配列群(l_配列群, l_ほかのk群, l_索引, C_リード長, out var l_統計, out var l_埋めた番号);

            Assert.Equal(l_配列群.Count, l_結果.Count);
            Assert.Equal(C_ID_S1, l_結果[0].A_ID);
            Assert.Equal(p_Is埋める ? l_G : l_元の配列, l_結果[0].A_配列);
            Assert.Equal((1, p_Is埋める ? 1 : 0, p_錨が取れない数, p_反復の端の数, p_ほかのkに無い数, p_リードに無い数), (l_統計.A_対象数, l_統計.A_埋めた数, l_統計.A_錨が取れない数, l_統計.A_反復の端の数, l_統計.A_ほかのkに無い数, l_統計.A_リードに無い数));
            string[] l_期待の鍵 = p_Is埋める ? [C_ID_S1] : [];
            Assert.Equal(l_期待の鍵, l_埋めた番号.Keys);
            if (p_Is埋める)
            {
                Assert.Equal([0], l_埋めた番号[C_ID_S1]);
            }
        }

        /// <summary>
        /// 埋めた番号: N の列が 2 つあり、2 つ目だけを埋めれば、埋めた番号は 1<br/>
        /// 1 つ目 (G[300..320) を N 20 個にした所) は左の片が 300 bp で反復を調べる長さより短いので、そのまま残し、反復の端の数は 1
        /// </summary>
        [Fact]
        public void 埋めた番号_2つ目のN列だけを埋めると番号は1()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_1つ目のN = new string(Consts.塩基コード.N, 20);
            List<(string A_ID, string A_配列)> l_配列群 = [(C_ID_S1, l_G[..300] + l_1つ目のN + l_G[320..C_N始め] + new string(Consts.塩基コード.N, C_N長) + l_G[C_N終わり..])];
            var l_リード = V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKGapFiller.Get_埋めた配列群(l_配列群, [l_G[C_ほかのk始め..C_ほかのk終わり]], l_索引, C_リード長, out var l_統計, out var l_埋めた番号);

            Assert.Equal(l_G[..300] + l_1つ目のN + l_G[320..], l_結果[0].A_配列);
            Assert.Equal((2, 1, 0, 1, 0, 0), (l_統計.A_対象数, l_統計.A_埋めた数, l_統計.A_錨が取れない数, l_統計.A_反復の端の数, l_統計.A_ほかのkに無い数, l_統計.A_リードに無い数));
            Assert.Equal([1], l_埋めた番号[C_ID_S1]);
        }

        /// <summary>
        /// 縦の反復の単位を抜いた形: G = A (1000 bp) + 単位 50 bp × 4 + B (1000 bp) で、ほかの k の contig が単位を 1 つ抜いた形 (× 3) なら埋めない<br/>
        /// 抜いた形の反復 (150 bp) はリード長の半分 (125 bp) の窓より長いので、その窓ではどれもリードに出てくるが、リード長の 7 割 (175 bp) の窓は反復の両脇を含みリードに無い
        /// </summary>
        [Fact]
        public void 縦の反復の単位を抜いた形_反復の両脇を含む窓がリードに無ければ埋めない()
        {
            var l_A = V_生成_ゲノム(1000, C_ゲノムの種);
            var l_単位 = V_生成_ゲノム(50, C_ゲノムの種 + 2);
            var l_B = V_生成_ゲノム(1000, C_ゲノムの種 + 3);
            var l_G = l_A + string.Concat(Enumerable.Repeat(l_単位, 4)) + l_B;
            var l_元の配列 = l_A + new string(Consts.塩基コード.N, C_N長) + l_B;
            List<(string A_ID, string A_配列)> l_配列群 = [(C_ID_S1, l_元の配列)];
            var l_リード = V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OtherKGapFiller.Get_埋めた配列群(l_配列群, [l_A[C_ほかのk始め..] + string.Concat(Enumerable.Repeat(l_単位, 3)) + l_B[..C_ほかのk始め]], l_索引, C_リード長, out var l_統計, out _);

            Assert.Equal(l_元の配列, l_結果[0].A_配列);
            Assert.Equal((1, 0, 0, 0, 0, 1), (l_統計.A_対象数, l_統計.A_埋めた数, l_統計.A_錨が取れない数, l_統計.A_反復の端の数, l_統計.A_ほかのkに無い数, l_統計.A_リードに無い数));
        }

        /// <summary>
        /// 埋めた後の番号: 埋めた番号を除き、それより後ろの番号を手前で埋めた数だけ詰める<br/>
        /// {0, 2, 3} から 0 と 2 を埋めれば {1}、{1, 4} から 0 と 2 を埋めれば {0, 2}
        /// </summary>
        /// <param name="p_長さ不明の番号">埋める前の、長さ不明のギャップの番号</param>
        /// <param name="p_埋めた番号">埋めた番号</param>
        /// <param name="p_期待">埋めた後の番号</param>
        [Theory]
        [InlineData(new[] { 0, 2, 3 }, new[] { 0, 2 }, new[] { 1 })]
        [InlineData(new[] { 1, 4 }, new[] { 0, 2 }, new[] { 0, 2 })]
        public void 埋めた後の番号_埋めた番号を除いて後ろを詰める(int[] p_長さ不明の番号, int[] p_埋めた番号, int[] p_期待)
        {
            var l_結果 = OtherKGapFiller.Get_埋めた後の番号(p_長さ不明の番号.ToHashSet(), p_埋めた番号.ToHashSet());

            Assert.Equal(p_期待.Order(), l_結果.Order());
        }

        #endregion

        #region 内部メソッド

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
