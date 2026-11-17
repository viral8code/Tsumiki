using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終配列の N の列を、深い錨とほかの k の contig の配列で繋ぐ処理 (OtherKDeepGapBridger) の検証<br/>
    /// 乱数の G (4,000 bp) の N の列の両脇を、重ねるか欠けさせて置く。リードは G の窓 (長さは場面ごと) を 5 bp おきに 2 本ずつ
    /// </summary>
    public class OtherKDeepGapBridgerTests
    {
        #region 定数

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 4000;

        /// <summary>
        /// ゲノムの乱数の種
        /// </summary>
        private const int C_ゲノムの種 = 7373;

        /// <summary>
        /// N の列の長さ
        /// </summary>
        private const int C_Nの列の長さ = 100;

        /// <summary>
        /// リードの窓をずらす幅
        /// </summary>
        private const int C_リードのずれ = 5;

        /// <summary>
        /// 場面 畳む (両脇の片が 50 bp 重なる。リード 150 bp)
        /// </summary>
        private const string C_場面_畳む = "畳む";

        /// <summary>
        /// 場面 リード長より長い (畳むと同じ配列。リード 100 bp で、錨から錨までの 112 bp がリード長より長い)
        /// </summary>
        private const string C_場面_リード長より長い = "リード長より長い";

        /// <summary>
        /// 場面 リードに無い (畳むと同じ配列。リード 110 bp、リード長 150 で、窓は通るが錨から錨までの 112 bp を含むリードが無い)
        /// </summary>
        private const string C_場面_リードに無い = "リードに無い";

        /// <summary>
        /// 場面 埋める (両脇の片の間の 100 bp が欠ける。リード 250 bp)
        /// </summary>
        private const string C_場面_埋める = "埋める";

        /// <summary>
        /// 場面 ほかの k 無し (畳むと同じ配列とリードで、ほかの k の contig が無い)
        /// </summary>
        private const string C_場面_ほかのk無し = "ほかの k 無し";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 繋いだ配列群: 錨から錨までがリード 1 本に収まり、リードに出てくるときだけ、ほかの k の錨の間の配列で繋いで G の 1 本にする<br/>
        /// 畳む形 (G[..2050) + N + G[2000..)、ほかの k は G[1500..2500)) も埋める形 (G[..2000) + N + G[2100..)、ほかの k は G[1500..2600)) も繋ぐ。<br/>
        /// 錨から錨までがリード長より長いとき、それを含むリードが無いとき、ほかの k が無いときは繋がない
        /// </summary>
        /// <param name="p_場面">場面</param>
        /// <param name="p_Is繋ぐ">繋ぐはずか</param>
        /// <param name="p_畳んだ数">畳んだ数</param>
        /// <param name="p_決まらない数">ほかの k で 1 通りに決まらない数</param>
        /// <param name="p_確かめられない数">錨から錨までをリードで確かめられない数</param>
        [Theory]
        [InlineData(C_場面_畳む, true, 1, 0, 0)]
        [InlineData(C_場面_リード長より長い, false, 0, 0, 1)]
        [InlineData(C_場面_リードに無い, false, 0, 0, 1)]
        [InlineData(C_場面_埋める, true, 0, 0, 0)]
        [InlineData(C_場面_ほかのk無し, false, 0, 1, 0)]
        public void 繋いだ配列群_錨から錨までがリードにあるときだけほかのkの配列で繋ぐ(string p_場面, bool p_Is繋ぐ, int p_畳んだ数, int p_決まらない数, int p_確かめられない数)
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_N = new string(Consts.塩基コード.N, C_Nの列の長さ);
            var l_配列 = p_場面 == C_場面_埋める ? l_G[..2000] + l_N + l_G[2100..] : l_G[..2050] + l_N + l_G[2000..];
            List<string> l_ほかのk群 = p_場面 switch
            {
                C_場面_ほかのk無し => [],
                C_場面_埋める => [l_G[1500..2600]],
                _ => [l_G[1500..2500]],
            };
            var (l_リードの長さ, l_リード長) = p_場面 switch
            {
                C_場面_リード長より長い => (100, 100),
                C_場面_リードに無い => (110, 150),
                C_場面_埋める => (250, 250),
                _ => (150, 150),
            };
            var l_索引 = V_構築_索引(l_G, l_リードの長さ);

            var l_結果 = OtherKDeepGapBridger.Get_繋いだ配列群([("S", l_配列)], l_ほかのk群, l_索引, l_リード長, out var l_統計, out var l_繋いだ番号);

            Assert.Equal(p_Is繋ぐ ? l_G : l_配列, l_結果[0].A_配列);
            Assert.Equal((1, p_Is繋ぐ ? 1 : 0, p_畳んだ数, p_決まらない数, p_確かめられない数), (l_統計.A_対象数, l_統計.A_繋いだ数, l_統計.A_畳んだ数, l_統計.A_決まらない数, l_統計.A_リードで確かめられない数));
            List<int> l_期待する番号 = p_Is繋ぐ ? [0] : [];
            List<int> l_番号 = l_繋いだ番号.TryGetValue("S", out var l_番号群) ? [.. l_番号群.Order()] : [];
            Assert.Equal(l_期待する番号, l_番号);
        }

        /// <summary>
        /// 繋いだ配列群: 2 つの N の列が 450 bp の片を挟んで続き (G[..2050) + N + G[2000..2450) + N + G[2400..))、どちらも両脇の片が 50 bp 重なるとき、<br/>
        /// 2 つ目の左の片は 1 つ目で使った右の錨から始め、2 つとも畳んで G の 1 本にする (ほかの k は G[1500..3000)、リード 150 bp)
        /// </summary>
        [Fact]
        public void 繋いだ配列群_続く2つのNの列を前から順に畳む()
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_N = new string(Consts.塩基コード.N, C_Nの列の長さ);
            var l_配列 = l_G[..2050] + l_N + l_G[2000..2450] + l_N + l_G[2400..];
            var l_索引 = V_構築_索引(l_G, 150);

            var l_結果 = OtherKDeepGapBridger.Get_繋いだ配列群([("S", l_配列)], [l_G[1500..3000]], l_索引, 150, out var l_統計, out var l_繋いだ番号);

            Assert.Equal(l_G, l_結果[0].A_配列);
            Assert.Equal((2, 2, 2, 0, 0), (l_統計.A_対象数, l_統計.A_繋いだ数, l_統計.A_畳んだ数, l_統計.A_決まらない数, l_統計.A_リードで確かめられない数));
            Assert.Equal([0, 1], l_繋いだ番号["S"].Order().ToList());
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
        /// 配列の窓 (長さ p_リードの長さ) を C_リードのずれ おきに 2 本ずつ切り出し、リードの索引を作る
        /// </summary>
        /// <param name="p_配列">切り出す元の配列</param>
        /// <param name="p_リードの長さ">リードの長さ</param>
        /// <returns>リードの索引</returns>
        private static ReadMinimizerIndex V_構築_索引(string p_配列, int p_リードの長さ)
        {
            List<string> l_リード = [.. Enumerable.Range(0, ((p_配列.Length - p_リードの長さ) / C_リードのずれ) + 1).Select(i => p_配列.Substring(i * C_リードのずれ, p_リードの長さ))];
            return ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());
        }

        #endregion
    }
}
