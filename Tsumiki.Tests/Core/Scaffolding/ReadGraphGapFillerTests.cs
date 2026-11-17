using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終配列の N の列を、リードに 2 本以上ある k-mer の 1 本の道で埋める処理 (ReadGraphGapFiller) の検証<br/>
    /// 乱数の G (4,000 bp) の N の列の両脇を欠けさせるか重ねて置く。リードは G の窓 (長さは場面ごと) を 5 bp おきに 2 本ずつ、リード長 150 (k = 105、75、52)。<br/>
    /// 縦の反復の場面は、乱数の単位 U を G の 2,000 bp 目に並べた T のリードで、道を比べる窓は 127 bp
    /// </summary>
    public class ReadGraphGapFillerTests
    {
        #region 定数

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 4000;

        /// <summary>
        /// ゲノムの乱数の種
        /// </summary>
        private const int C_ゲノムの種 = 8484;

        /// <summary>
        /// N の列の長さ
        /// </summary>
        private const int C_Nの列の長さ = 100;

        /// <summary>
        /// 代表リード長
        /// </summary>
        private const int C_リード長 = 150;

        /// <summary>
        /// リードの窓をずらす幅
        /// </summary>
        private const int C_リードのずれ = 5;

        /// <summary>
        /// 場面 埋める (G[..2000) + N + G[2060..)、リード 150 bp)
        /// </summary>
        private const string C_場面_埋める = "埋める";

        /// <summary>
        /// 場面 片の端の誤り (左の片の最後の塩基が G と違う)
        /// </summary>
        private const string C_場面_片の端の誤り = "片の端の誤り";

        /// <summary>
        /// 場面 リード長より長い (G[..2000) + N + G[2100..)。錨から錨まで 162 bp)
        /// </summary>
        private const string C_場面_リード長より長い = "リード長より長い";

        /// <summary>
        /// 場面 リードに無い (埋めると同じ配列で、リード 110 bp。錨から錨まで 122 bp を含むリードが無い)
        /// </summary>
        private const string C_場面_リードに無い = "リードに無い";

        /// <summary>
        /// 場面 道が 2 本 (埋めると同じ配列で、G[2030] を変えた G' のリードもある)
        /// </summary>
        private const string C_場面_道が2本 = "道が 2 本";

        /// <summary>
        /// 場面 道が無い (埋めると同じ配列で、G[2030] を含むリードが無い)
        /// </summary>
        private const string C_場面_道が無い = "道が無い";

        /// <summary>
        /// 場面 畳む形 (G[..2050) + N + G[2000..)。畳んだ形の深い錨から深い錨まで G[1969..2081) の 112 bp がリードにある)
        /// </summary>
        private const string C_場面_畳む形 = "畳む形";

        /// <summary>
        /// 場面 畳む形で錨から錨までが長い (G[..2100) + N + G[2000..)。深い錨から深い錨まで G[1969..2131) の 162 bp がリード長より長い)
        /// </summary>
        private const string C_場面_畳む形で錨から錨までが長い = "畳む形で錨から錨までが長い";

        /// <summary>
        /// 場面 k を下げて埋める (G[..2000) + N + G[2010..)、リード 100 bp。k 105 では端の k-mer がリードに無く、k 75 で埋める)
        /// </summary>
        private const string C_場面_kを下げて埋める = "k を下げて埋める";

        /// <summary>
        /// 場面 片が短い (G[1600..2000) + N + G[2060..))
        /// </summary>
        private const string C_場面_片が短い = "片が短い";

        /// <summary>
        /// 場面 窓で選んで埋める (U 10 bp × 12 の T のリード、配列 = G[..2000) + N + G[2000..))
        /// </summary>
        private const string C_場面_窓で選んで埋める = "窓で選んで埋める";

        /// <summary>
        /// 場面 反復が窓より長い (U 30 bp × 6 の T のリード。180 bp の反復が窓より長い)
        /// </summary>
        private const string C_場面_反復が窓より長い = "反復が窓より長い";

        /// <summary>
        /// 場面 両方の形のリードがある (U 10 bp × 12 の T と × 13 の T の両方のリード)
        /// </summary>
        private const string C_場面_両方の形のリードがある = "両方の形のリードがある";

        /// <summary>
        /// 縦の反復の単位の乱数の種
        /// </summary>
        private const int C_単位の種 = 4242;

        /// <summary>
        /// 縦の反復を挟む G の位置
        /// </summary>
        private const int C_反復の位置 = 2000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 埋めた配列群: 道が 1 本で、埋める形は錨から錨までが、畳む形は深い錨から深い錨までがリードで確かめられるときだけ、道で埋めて G の 1 本にする<br/>
        /// 片の端の k-mer がリードに無ければ端から内へずらし、k 105 で端の k-mer か道が無ければ k 75、52 と下げる。<br/>
        /// 錨から錨までがリード長より長いとき、それを含むリードが無いとき、道が 2 本か無いとき、片が短いときは埋めない
        /// </summary>
        /// <param name="p_場面">場面</param>
        /// <param name="p_Is埋める">埋めるはずか</param>
        /// <param name="p_理由">埋めない理由 (リードの道で埋めない理由 の値。埋めるなら -1)</param>
        /// <param name="p_k長">埋めた k か、最後に試した k</param>
        [Theory]
        [InlineData(C_場面_埋める, true, -1, 105)]
        [InlineData(C_場面_片の端の誤り, true, -1, 105)]
        [InlineData(C_場面_畳む形, true, -1, 105)]
        [InlineData(C_場面_kを下げて埋める, true, -1, 75)]
        [InlineData(C_場面_リード長より長い, false, (int)リードの道で埋めない理由.錨から錨までを確かめられない, 105)]
        [InlineData(C_場面_リードに無い, false, (int)リードの道で埋めない理由.錨から錨までを確かめられない, 105)]
        [InlineData(C_場面_道が2本, false, (int)リードの道で埋めない理由.道が決まらない, 105)]
        [InlineData(C_場面_道が無い, false, (int)リードの道で埋めない理由.道が無い, 52)]
        [InlineData(C_場面_畳む形で錨から錨までが長い, false, (int)リードの道で埋めない理由.畳む形を確かめられない, 105)]
        [InlineData(C_場面_片が短い, false, (int)リードの道で埋めない理由.片が短い, 105)]
        public void 埋めた配列群_道が1本でリードに確かめられるときだけ埋める(string p_場面, bool p_Is埋める, int p_理由, int p_k長)
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_N = new string(Consts.塩基コード.N, C_Nの列の長さ);
            var l_違う塩基 = l_G[1999] == Consts.塩基コード.A ? Consts.塩基コード.C : Consts.塩基コード.A;
            var l_配列 = p_場面 switch
            {
                C_場面_片の端の誤り => l_G[..1999] + l_違う塩基 + l_N + l_G[2060..],
                C_場面_リード長より長い => l_G[..2000] + l_N + l_G[2100..],
                C_場面_畳む形 => l_G[..2050] + l_N + l_G[2000..],
                C_場面_畳む形で錨から錨までが長い => l_G[..2100] + l_N + l_G[2000..],
                C_場面_kを下げて埋める => l_G[..2000] + l_N + l_G[2010..],
                C_場面_片が短い => l_G[1600..2000] + l_N + l_G[2060..],
                _ => l_G[..2000] + l_N + l_G[2060..],
            };
            var l_違う形のG = l_G[..2030] + (l_G[2030] == Consts.塩基コード.A ? Consts.塩基コード.C : Consts.塩基コード.A) + l_G[2031..];
            List<string> l_リード = p_場面 switch
            {
                C_場面_リードに無い => V_リードの窓(l_G, 110),
                C_場面_kを下げて埋める => V_リードの窓(l_G, 100),
                C_場面_道が2本 => [.. V_リードの窓(l_G, 150), .. V_リードの窓(l_違う形のG, 150)],
                C_場面_道が無い => [.. V_リードの窓(l_G, 150).Where((_, i) => i * C_リードのずれ > 2030 || (i * C_リードのずれ) + 150 <= 2030)],
                _ => V_リードの窓(l_G, 150),
            };
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = ReadGraphGapFiller.Get_埋めた配列群([("S", l_配列)], l_索引, C_リード長, out var l_統計, out var l_埋めた番号);

            Assert.Equal(p_Is埋める ? l_G[(p_場面 == C_場面_片が短い ? 1600 : 0)..] : l_配列, l_結果[0].A_配列);
            int[] l_理由ごと = [l_統計.A_片が短い数, l_統計.A_端がリードに無い数, l_統計.A_道が無い数, l_統計.A_道が決まらない数, l_統計.A_畳む形を確かめられない数, l_統計.A_錨から錨までを確かめられない数];
            Assert.Equal([.. Enumerable.Range(0, 6).Select(i => i == p_理由 ? 1 : 0)], l_理由ごと);
            Assert.Equal([105, 75, 52], l_統計.A_k長群);
            Assert.Equal((1, p_Is埋める ? 1 : 0, p_場面 == C_場面_畳む形 ? 1 : 0, p_Is埋める && p_k長 != 105 ? 1 : 0, 0), (l_統計.A_対象数, l_統計.A_埋めた数, l_統計.A_畳んだ数, l_統計.A_小さいkで埋めた数, l_統計.A_窓で選んだ数));
            Assert.Equal(p_Is埋める, l_埋めた番号.ContainsKey("S"));
        }

        /// <summary>
        /// 埋めた配列群: 道が 2 本以上の縦の反復では、短い順に並べた道をリード長の窓で比べ、1 本に決まるときだけその道で埋めて T の 1 本にする<br/>
        /// 反復が窓より長いと単位の数を分ける窓が無く、両方の形のリードがあると相手の形にだけある窓にリードがあるので、選ばずに道が決まらないとする
        /// </summary>
        /// <param name="p_場面">場面</param>
        /// <param name="p_Is埋める">埋めるはずか</param>
        [Theory]
        [InlineData(C_場面_窓で選んで埋める, true)]
        [InlineData(C_場面_反復が窓より長い, false)]
        [InlineData(C_場面_両方の形のリードがある, false)]
        public void 埋めた配列群_縦の反復は窓で1本に決まるときだけ埋める(string p_場面, bool p_Is埋める)
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_単位 = V_生成_ゲノム(p_場面 == C_場面_反復が窓より長い ? 30 : 10, C_単位の種);
            var l_単位の数 = p_場面 == C_場面_反復が窓より長い ? 6 : 12;
            var l_T = Get_反復を挟んだ配列(l_G, l_単位, l_単位の数);
            var l_配列 = l_G[..C_反復の位置] + new string(Consts.塩基コード.N, C_Nの列の長さ) + l_G[C_反復の位置..];
            List<string> l_リード = p_場面 == C_場面_両方の形のリードがある ? [.. V_リードの窓(l_T, 150), .. V_リードの窓(Get_反復を挟んだ配列(l_G, l_単位, l_単位の数 + 1), 150)] : V_リードの窓(l_T, 150);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = ReadGraphGapFiller.Get_埋めた配列群([("S", l_配列)], l_索引, C_リード長, out var l_統計, out _);

            Assert.Equal(p_Is埋める ? l_T : l_配列, l_結果[0].A_配列);
            Assert.Equal((1, p_Is埋める ? 1 : 0, p_Is埋める ? 1 : 0, p_Is埋める ? 0 : 1), (l_統計.A_対象数, l_統計.A_埋めた数, l_統計.A_窓で選んだ数, l_統計.A_道が決まらない数));
        }

        /// <summary>
        /// 短い順の道群: U 10 bp × 12 の縦の反復を挟む始まりと終わりの k-mer から、U × 11 から 10 bp ずつ長い道を 12 本、短い順に並べる
        /// </summary>
        [Fact]
        public void 短い順の道群_縦の反復の周回を短い順に12本並べる()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_単位 = V_生成_ゲノム(10, C_単位の種);
            var l_リード = V_リードの窓(Get_反復を挟んだ配列(l_G, l_単位, 12), 150);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());
            var l_始まり = l_G[(C_反復の位置 - 105)..C_反復の位置];
            var l_終わり = l_G[C_反復の位置..(C_反復の位置 + 105)];

            var l_道群 = ReadGraphGapFiller.Get_短い順の道群(l_始まり, l_終わり, (2 * 105) + 2000, l_索引, out var l_Is打ち切り);

            Assert.False(l_Is打ち切り);
            Assert.Equal([.. Enumerable.Range(11, 12).Select(n => l_始まり + string.Concat(Enumerable.Repeat(l_単位, n)) + l_終わり)], l_道群);
        }

        /// <summary>
        /// 埋める形: 左の片の最後の塩基が G と違うと、末尾の k-mer がリードに無いので 10 bp 内へずらし、道を G[1885..2165) にする
        /// </summary>
        [Fact]
        public void 埋める形_片の端のkmerがリードに無ければ10bpずらす()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_左 = l_G[..1999] + (l_G[1999] == Consts.塩基コード.A ? Consts.塩基コード.C : Consts.塩基コード.A);
            var l_リード = V_リードの窓(l_G, 150);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_形 = ReadGraphGapFiller.Get_埋める形(l_左, l_G[2060..], [], l_索引, C_リード長, 105, out _, out _);

            Assert.True(l_形.HasValue);
            Assert.Equal((10, l_G[1885..2165], 0), l_形.Value);
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
        /// 配列の窓 (長さ p_リードの長さ) を C_リードのずれ おきに切り出す
        /// </summary>
        /// <param name="p_配列">切り出す元の配列</param>
        /// <param name="p_リードの長さ">リードの長さ</param>
        /// <returns>リードの一覧</returns>
        private static List<string> V_リードの窓(string p_配列, int p_リードの長さ)
        {
            return [.. Enumerable.Range(0, ((p_配列.Length - p_リードの長さ) / C_リードのずれ) + 1).Select(i => p_配列.Substring(i * C_リードのずれ, p_リードの長さ))];
        }

        /// <summary>
        /// G の C_反復の位置 に、単位を並べた縦の反復を挟んだ配列を作る
        /// </summary>
        /// <param name="p_G">G</param>
        /// <param name="p_単位">反復の単位</param>
        /// <param name="p_単位の数">単位を並べる数</param>
        /// <returns>反復を挟んだ配列</returns>
        private static string Get_反復を挟んだ配列(string p_G, string p_単位, int p_単位の数)
        {
            return p_G[..C_反復の位置] + string.Concat(Enumerable.Repeat(p_単位, p_単位の数)) + p_G[C_反復の位置..];
        }

        #endregion
    }
}
