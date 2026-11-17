using Tsumiki.Commons;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終配列の N の列のうち、両脇の片が完全に重なる所をリードに確かめて畳む処理 (OverlapGapFolder) の検証<br/>
    /// リードは正しいゲノムの 250 bp の窓 (5 bp おき) を 2 回ずつ、リード長 250 (窓 175 bp)
    /// </summary>
    public class OverlapGapFolderTests
    {
        #region 定数

        /// <summary>
        /// ゲノムの長さ
        /// </summary>
        private const int C_ゲノム長 = 2000;

        /// <summary>
        /// ゲノムの乱数の種
        /// </summary>
        private const int C_ゲノムの種 = 5151;

        /// <summary>
        /// 代表リード長
        /// </summary>
        private const int C_リード長 = 250;

        /// <summary>
        /// リードの窓をずらす幅
        /// </summary>
        private const int C_リードのずれ = 5;

        /// <summary>
        /// 左の片の終わり (重なりの終わり) の G の位置
        /// </summary>
        private const int C_左の終わり = 1060;

        /// <summary>
        /// 右の片の始まりの G の位置
        /// </summary>
        private const int C_右の始め = 1000;

        /// <summary>
        /// アセンブリの N の列の長さ
        /// </summary>
        private const int C_N長 = 100;

        /// <summary>
        /// 縦の反復の両脇の長さ
        /// </summary>
        private const int C_脇の長さ = 1000;

        /// <summary>
        /// 配列 ID S1
        /// </summary>
        private const string C_ID_S1 = "S1";

        /// <summary>
        /// 場面 重なり
        /// </summary>
        private const string C_場面_重なり = "重なり";

        /// <summary>
        /// 場面 長い重なり
        /// </summary>
        private const string C_場面_長い重なり = "長い重なり";

        /// <summary>
        /// 場面 リード無し
        /// </summary>
        private const string C_場面_リード無し = "リード無し";

        /// <summary>
        /// 場面 片が短い
        /// </summary>
        private const string C_場面_片が短い = "片が短い";

        /// <summary>
        /// 場面 重なりが無い
        /// </summary>
        private const string C_場面_重なりが無い = "重なりが無い";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 畳むか: 乱数の G で、アセンブリが G[..1060) + N×100 + G[1000..) (重なり 60 bp) なら、畳んで G にする。重なりが 180 bp (G[..1180)) でも、重なりと両脇 1 塩基を含むリードがあるので G にする<br/>
        /// リードが無いとき、左の片が 400 bp で反復を調べる長さより短いとき、重なりが無い (G[..1000) + N×100 + G[1100..)) ときは畳まず、理由ごとの数を数える
        /// </summary>
        /// <param name="p_場面">場面</param>
        /// <param name="p_Is畳む">畳んで G にするはずか</param>
        /// <param name="p_片が短い数">期待する片が短い数</param>
        /// <param name="p_重なりが無い数">期待する重なりが無い数</param>
        /// <param name="p_リードに無い数">期待するリードに無い数</param>
        [Theory]
        [InlineData(C_場面_重なり, true, 0, 0, 0)]
        [InlineData(C_場面_長い重なり, true, 0, 0, 0)]
        [InlineData(C_場面_リード無し, false, 0, 0, 1)]
        [InlineData(C_場面_片が短い, false, 1, 0, 0)]
        [InlineData(C_場面_重なりが無い, false, 0, 1, 0)]
        public void 畳むか_両脇の片が重なり畳んだ形がリードに続けて出るときだけ畳む(string p_場面, bool p_Is畳む, int p_片が短い数, int p_重なりが無い数, int p_リードに無い数)
        {
            var l_G = V_生成_ゲノム(C_ゲノム長, C_ゲノムの種);
            var l_N = new string(Consts.塩基コード.N, C_N長);
            var l_元の配列 = p_場面 switch
            {
                C_場面_長い重なり => l_G[..1180] + l_N + l_G[C_右の始め..],
                C_場面_片が短い => l_G[660..C_左の終わり] + l_N + l_G[C_右の始め..],
                C_場面_重なりが無い => l_G[..C_右の始め] + l_N + l_G[1100..],
                _ => l_G[..C_左の終わり] + l_N + l_G[C_右の始め..],
            };
            List<string> l_リード = p_場面 == C_場面_リード無し ? [] : V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OverlapGapFolder.Get_畳んだ配列群([(C_ID_S1, l_元の配列)], l_索引, C_リード長, out var l_統計, out var l_畳んだ番号);

            Assert.Equal(C_ID_S1, Assert.Single(l_結果).A_ID);
            Assert.Equal(p_Is畳む ? l_G : l_元の配列, l_結果[0].A_配列);
            Assert.Equal((1, p_Is畳む ? 1 : 0, p_片が短い数, p_重なりが無い数, 0, p_リードに無い数, 0), (l_統計.A_対象数, l_統計.A_畳んだ数, l_統計.A_片が短い数, l_統計.A_重なりが無い数, l_統計.A_反復の端の数, l_統計.A_リードに無い数, l_統計.A_決まらない数));
            string[] l_期待の鍵 = p_Is畳む ? [C_ID_S1] : [];
            Assert.Equal(l_期待の鍵, l_畳んだ番号.Keys);
        }

        /// <summary>
        /// X・間・X: 正しいゲノムが A + X (60 bp) + 間 (113 bp) + X + B で、アセンブリが A + X + N×100 + X + B なら、<br/>
        /// 両脇の片は X で完全に重なるが、畳んだ形 (A + X + B) の、A の末尾・X・B の先頭を含む窓がリードに無いので畳まない (リードに無い 1)
        /// </summary>
        [Fact]
        public void X間X_畳んだ形の窓がリードに無ければ畳まない()
        {
            var l_A = V_生成_ゲノム(C_脇の長さ, C_ゲノムの種 + 1);
            var l_X = V_生成_ゲノム(60, C_ゲノムの種 + 2);
            var l_間 = V_生成_ゲノム(113, C_ゲノムの種 + 3);
            var l_B = V_生成_ゲノム(C_脇の長さ, C_ゲノムの種 + 4);
            var l_リード = V_リードの窓(l_A + l_X + l_間 + l_X + l_B);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());
            var l_元の配列 = l_A + l_X + new string(Consts.塩基コード.N, C_N長) + l_X + l_B;

            var l_結果 = OverlapGapFolder.Get_畳んだ配列群([(C_ID_S1, l_元の配列)], l_索引, C_リード長, out var l_統計, out _);

            Assert.Equal(l_元の配列, l_結果[0].A_配列);
            Assert.Equal((0, 1), (l_統計.A_畳んだ数, l_統計.A_リードに無い数));
        }

        /// <summary>
        /// 縦の反復の中: 正しいゲノムが A + 単位 100 bp × 5 + B で、アセンブリが (A + 単位 × 2) + N×100 + (単位 × 2 + B) なら畳まない<br/>
        /// 重なり 100 bp は、重なりを除いても両脇が単位を共有するので反復の端。重なり 200 bp (単位 2 つ) は、窓はどれもリードにあるが、<br/>
        /// 重なりと両脇 1 塩基 (A の末尾 + 単位 × 2 + B の先頭) を含むリードが無い (リードに無い 1)
        /// </summary>
        [Fact]
        public void 縦の反復の中_単位を抜いて畳む重なりは使わない()
        {
            var l_A = V_生成_ゲノム(C_脇の長さ, C_ゲノムの種 + 5);
            var l_単位 = V_生成_ゲノム(100, C_ゲノムの種 + 6);
            var l_B = V_生成_ゲノム(C_脇の長さ, C_ゲノムの種 + 7);
            var l_2つ = l_単位 + l_単位;
            var l_リード = V_リードの窓(l_A + string.Concat(Enumerable.Repeat(l_単位, 5)) + l_B);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());
            var l_元の配列 = l_A + l_2つ + new string(Consts.塩基コード.N, C_N長) + l_2つ + l_B;

            var l_結果 = OverlapGapFolder.Get_畳んだ配列群([(C_ID_S1, l_元の配列)], l_索引, C_リード長, out var l_統計, out _);

            Assert.Equal(l_元の配列, l_結果[0].A_配列);
            Assert.Equal((0, 0, 1), (l_統計.A_畳んだ数, l_統計.A_反復の端の数, l_統計.A_リードに無い数));
        }

        /// <summary>
        /// 畳んだ番号: N の列が 2 つあり、1 つ目 (G[600..700) を N×100 にした所) は重なりが無く、2 つ目 (重なり 60 bp) だけを畳めば、畳んだ番号は 1
        /// </summary>
        [Fact]
        public void 畳んだ番号_2つ目のN列だけを畳むと番号は1()
        {
            var l_G = V_生成_ゲノム(2600, C_ゲノムの種 + 8);
            var l_N = new string(Consts.塩基コード.N, C_N長);
            var l_リード = V_リードの窓(l_G);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード.Concat(l_リード).ToList());

            var l_結果 = OverlapGapFolder.Get_畳んだ配列群([(C_ID_S1, l_G[..600] + l_N + l_G[700..1660] + l_N + l_G[1600..])], l_索引, C_リード長, out var l_統計, out var l_畳んだ番号);

            Assert.Equal(l_G[..600] + l_N + l_G[700..], l_結果[0].A_配列);
            Assert.Equal((2, 1, 1), (l_統計.A_対象数, l_統計.A_畳んだ数, l_統計.A_重なりが無い数));
            Assert.Equal([1], l_畳んだ番号[C_ID_S1]);
        }

        /// <summary>
        /// 一致する重なり群: 左の片の末尾と右の片の先頭が Y + Z + Y (Y 25 bp、Z 15 bp) で終わり・始まるとき、Y (25 bp) と Y + Z + Y (65 bp) の 2 つの一致を短い順に返す
        /// </summary>
        [Fact]
        public void 一致する重なり群_下限以上の一致を短い順に返す()
        {
            var l_Y = V_生成_ゲノム(25, C_ゲノムの種 + 9);
            var l_Z = V_生成_ゲノム(15, C_ゲノムの種 + 10);
            var l_左 = V_生成_ゲノム(100, C_ゲノムの種 + 11) + l_Y + l_Z + l_Y;
            var l_右 = l_Y + l_Z + l_Y + V_生成_ゲノム(100, C_ゲノムの種 + 12);

            Assert.Equal([25, 65], OverlapGapFolder.Get_一致する重なり群(l_左, l_右));
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
