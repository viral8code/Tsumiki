using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアの橋渡しで反復の端として除いた辺を、ほかの k の contig の配列で繋ぐ形 (OtherKEdgeBridger) の検証<br/>
    /// 乱数の G (4,000 bp) で、u = G[..2800)、x = u の末尾 450 bp の写し + G[2850..)。リードは G の 250 bp の窓 (5 bp おき) を 2 回ずつ、リード長 250
    /// </summary>
    public class OtherKEdgeBridgerTests
    {
        #region 定数

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 4000;

        /// <summary>
        /// ゲノムの乱数の種
        /// </summary>
        private const int C_ゲノムの種 = 6262;

        /// <summary>
        /// u の長さ
        /// </summary>
        private const int C_u長 = 2800;

        /// <summary>
        /// x の写しの後ろが始まる G の位置
        /// </summary>
        private const int C_x開始 = 2850;

        /// <summary>
        /// x の先頭に置く、u の末尾の写しの長さ
        /// </summary>
        private const int C_写しの長さ = 450;

        /// <summary>
        /// 代表リード長
        /// </summary>
        private const int C_リード長 = 250;

        /// <summary>
        /// リードの窓をずらす幅
        /// </summary>
        private const int C_リードのずれ = 5;

        /// <summary>
        /// 場面 繋ぐ
        /// </summary>
        private const string C_場面_繋ぐ = "繋ぐ";

        /// <summary>
        /// 場面 ほかの k 無し
        /// </summary>
        private const string C_場面_ほかのk無し = "ほかの k 無し";

        /// <summary>
        /// 場面 リード無し
        /// </summary>
        private const string C_場面_リード無し = "リード無し";

        /// <summary>
        /// 場面 2 通り
        /// </summary>
        private const string C_場面_2通り = "2 通り";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 繋ぐ形: ほかの k の contig (G[1800..3850)) が u の末尾と x の先頭の一意な錨の間を繋ぎ、その配列がリードに続けて出てくれば、u と x を繋いで G にする形を返す<br/>
        /// ほかの k が無いとき、リードが無いとき、1 塩基違う 2 本目の contig があり両方の形のリードがあるときは返さない
        /// </summary>
        /// <param name="p_場面">場面</param>
        /// <param name="p_Is繋ぐ">繋ぐ形を返すはずか</param>
        [Theory]
        [InlineData(C_場面_繋ぐ, true)]
        [InlineData(C_場面_ほかのk無し, false)]
        [InlineData(C_場面_リード無し, false)]
        [InlineData(C_場面_2通り, false)]
        public void 繋ぐ形_ほかのkの配列がただ1つリードに続けて出るときだけ返す(string p_場面, bool p_Is繋ぐ)
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_u = l_G[..C_u長];
            var l_x = l_u[(C_u長 - C_写しの長さ)..] + l_G[C_x開始..];
            var l_違う形のG = l_G[..2825] + (l_G[2825] == Consts.塩基コード.A ? Consts.塩基コード.C : Consts.塩基コード.A) + l_G[2826..];
            List<string> l_ほかのk群 = p_場面 switch
            {
                C_場面_ほかのk無し => [],
                C_場面_2通り => [l_G[1800..3850], l_違う形のG[1800..3850]],
                _ => [l_G[1800..3850]],
            };
            List<string> l_リード = p_場面 switch
            {
                C_場面_リード無し => [],
                C_場面_2通り => [.. V_リードの窓(l_G), .. V_リードの窓(l_違う形のG)],
                _ => V_リードの窓(l_G),
            };
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());
            var l_錨表 = OtherKRepeatReplacer.Get_一意な錨の表([("u", l_u), ("x", l_x)]);

            var l_形 = OtherKEdgeBridger.Get_繋ぐ形(l_u, l_x, l_錨表, string.Join('#', l_ほかのk群), l_索引, C_リード長);

            Assert.Equal(p_Is繋ぐ, l_形.HasValue);
            if (l_形 is { } l_繋ぐ形)
            {
                Assert.Equal(l_G, l_u[..^l_繋ぐ形.A_u側] + l_繋ぐ形.A_埋める配列 + l_x[l_繋ぐ形.A_x側..]);
            }
        }

        /// <summary>
        /// 錨の位置群: どの 31-mer も一意な乱数の配列 (600 bp) では、前後 30 bp に始まる 31-mer も一意な所から、末尾からは 539, 514, 489, … と 25 bp おきに 20 個、先頭からは 30, 55, 80, … と 20 個を返す。<br/>
        /// 末尾 100 bp の写しが別の配列にあると、500 以降に始まる 31-mer は一意でないので、末尾からの最初は 469 になる
        /// </summary>
        [Fact]
        public void 錨の位置群_一意な錨を端に近い順に25bpおきに返す()
        {
            var l_X = V_生成_ゲノム(600, C_ゲノムの種 + 1);
            var l_表 = OtherKRepeatReplacer.Get_一意な錨の表([("X", l_X)]);
            var l_写しありの表 = OtherKRepeatReplacer.Get_一意な錨の表([("X", l_X), ("Y", l_X[500..])]);

            Assert.Equal([.. Enumerable.Range(0, 20).Select(i => 539 - (25 * i))], OtherKEdgeBridger.Get_錨の位置群(l_X, l_表, true));
            Assert.Equal([.. Enumerable.Range(0, 20).Select(i => 30 + (25 * i))], OtherKEdgeBridger.Get_錨の位置群(l_X, l_表, false));
            Assert.Equal(469, OtherKEdgeBridger.Get_錨の位置群(l_X, l_写しありの表, true)[0]);
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
