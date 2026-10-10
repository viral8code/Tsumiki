using Tsumiki.Commons;
using System.Text;
using Tsumiki.Cores.Pipeline;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終成果物で、未確認の繋ぎ目を長さ不明のギャップにし、AGP に書き分ける処理の検証
    /// </summary>
    public class UnknownGapTests
    {
        #region 定数

        /// <summary>
        /// GUID 書式
        /// </summary>
        private const string C_GUID書式 = "N";

        /// <summary>
        /// 項目 n
        /// </summary>
        private const string C_項目_n = "n";

        /// <summary>
        /// 塩基配列 ACGTnACNNNGT
        /// </summary>
        private const string C_塩基配列_ACGTnACNNNGT = "ACGTnACNNNGT";

        /// <summary>
        /// 項目 S1
        /// </summary>
        private const string C_項目_S1 = "S1";

        /// <summary>
        /// 項目 W
        /// </summary>
        private const string C_項目_W = "W";

        /// <summary>
        /// 項目 U
        /// </summary>
        private const string C_項目_U = "U";

        /// <summary>
        /// 項目 100
        /// </summary>
        private const string C_項目_100 = "100";

        /// <summary>
        /// 項目 3
        /// </summary>
        private const string C_項目_3 = "3";

        /// <summary>
        /// 項目 S1 1
        /// </summary>
        private const string C_項目_S1_1 = "S1_1";

        /// <summary>
        /// 項目 S1 2
        /// </summary>
        private const string C_項目_S1_2 = "S1_2";

        /// <summary>
        /// 項目 S1 3
        /// </summary>
        private const string C_項目_S1_3 = "S1_3";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 未確認の繋ぎ目の印だけを N 100 個にし、推定長の N はそのまま残す
        /// </summary>
        /// <param name="p_配列">置き換える前の配列</param>
        /// <param name="p_長さ不明の番号">長さ不明になるはずの N の連続の番号</param>
        [Theory]
        [InlineData("ACnGT", new[] { 0 })]
        [InlineData("ACNGTnA", new[] { 1 })]
        [InlineData("ACnGTNNNAnC", new[] { 0, 2 })]
        [InlineData("ACNGT", new int[0])]
        public void Get_長さ不明のギャップへ置換_未確認の繋ぎ目だけを置き換える(string p_配列, int[] p_長さ不明の番号)
        {
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(p_配列, l_番号群);
            Assert.Equal(p_長さ不明の番号, l_番号群.Order());
            Assert.DoesNotContain(Consts.未確認の繋ぎ目, l_結果);
            Assert.Equal(p_配列.Length + (p_長さ不明の番号.Length * (FinalAssemblyPipeline.C_長さ不明のギャップ長 - 1)), l_結果.Length);
            Assert.Equal(p_配列.Replace(C_GUID書式, string.Empty).Replace(C_項目_n, string.Empty), l_結果.Replace(C_GUID書式, string.Empty));
        }

        /// <summary>
        /// AGP では、長さ不明のギャップを U、推定長のギャップを N、配列の片を W として書く
        /// </summary>
        [Fact]
        public void Get_AGP行_ギャップの種類を書き分ける()
        {
            var l_番号群 = new HashSet<int>();
            var l_配列 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(C_塩基配列_ACGTnACNNNGT, l_番号群);
            var l_行群 = FinalAssemblyPipeline.Get_AGP行(C_項目_S1, l_配列, l_番号群).Select(x => x.Split('\t')).ToList();
            Assert.Equal([C_項目_W, C_項目_U, C_項目_W, C_GUID書式, C_項目_W], l_行群.Select(x => x[4]));
            Assert.Equal(C_項目_100, l_行群[1][5]);
            Assert.Equal(C_項目_3, l_行群[3][5]);
            Assert.Equal([C_項目_S1_1, C_項目_S1_2, C_項目_S1_3], l_行群.Where(x => x[4] == C_項目_W).Select(x => x[5]));
            Assert.Equal(l_配列.Length.ToString(), l_行群[^1][2]);
        }

        /// <summary>
        /// 未確認の繋ぎ目は、重なりをリードで確かめられれば畳み、確かめられなければ印のまま残す
        /// </summary>
        /// <param name="p_リード数">繋いだ配列を含むリードの本数</param>
        /// <param name="p_Is畳む">畳むはずか</param>
        [Theory]
        [InlineData(3, true)]
        [InlineData(1, false)]
        public void Get_確かめた繋ぎ目を畳んだ配列_リードで確かめられれば畳む(int p_リード数, bool p_Is畳む)
        {
            var l_乱数 = new Random(11);
            var l_共通 = Get_乱配列(l_乱数, 30);
            var l_左 = Get_乱配列(l_乱数, 200) + l_共通;
            var l_右 = l_共通 + Get_乱配列(l_乱数, 200);
            var l_繋いだ配列 = l_左 + l_右[30..];
            var l_リード群 = Enumerable.Range(0, p_リード数).Select(i => l_繋いだ配列.Substring(170 + i, 100)).ToList();
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);
            var l_総数 = 0;
            var l_畳んだ数 = 0;
            var l_埋めた数 = 0;
            var l_結果 = FinalAssemblyPipeline.Get_確かめた繋ぎ目を畳んだ配列(l_左 + Consts.未確認の繋ぎ目 + l_右, l_索引, 89, 100, 0, ref l_総数, ref l_畳んだ数, ref l_埋めた数);
            Assert.Equal(1, l_総数);
            Assert.Equal(p_Is畳む ? 1 : 0, l_畳んだ数);
            Assert.Equal(p_Is畳む ? l_繋いだ配列 : l_左 + Consts.未確認の繋ぎ目 + l_右, l_結果);
        }

        /// <summary>
        /// 重なりの無い短い隙間は、両側を跨ぐリードの続きで埋める (リードの向きによらない)
        /// </summary>
        [Fact]
        public void Get_確かめた繋ぎ目を畳んだ配列_跨ぐリードで隙間を埋める()
        {
            var l_乱数 = new Random(21);
            var l_左 = Get_乱配列(l_乱数, 200);
            var l_隙間 = Get_乱配列(l_乱数, 30);
            var l_右 = Get_乱配列(l_乱数, 200);
            var l_本当 = l_左 + l_隙間 + l_右;
            var l_リード群 = Enumerable.Range(0, 8).Select(i => l_本当.Substring(130 + (i * 3), 150)).Select((x, i) => i % 2 == 0 ? x : Util.V_逆相補(x)).ToList();
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);
            var l_総数 = 0;
            var l_畳んだ数 = 0;
            var l_埋めた数 = 0;
            var l_結果 = FinalAssemblyPipeline.Get_確かめた繋ぎ目を畳んだ配列(l_左 + Consts.未確認の繋ぎ目 + l_右, l_索引, 63, 150, 0, ref l_総数, ref l_畳んだ数, ref l_埋めた数);
            Assert.Equal((1, 0, 1), (l_総数, l_畳んだ数, l_埋めた数));
            Assert.Equal(l_本当, l_結果);
        }

        /// <summary>
        /// 縦に並んだ反復の単位で右の片が始まる繋ぎ目は、リードの続きに右の先頭が先に現れても埋めない (単位の数を取り違えるため)
        /// </summary>
        [Fact]
        public void Get_リードで埋めた繋ぎ目_縦に並んだ反復の中は埋めない()
        {
            var l_乱数 = new Random(22);
            var l_左 = Get_乱配列(l_乱数, 200);
            var l_単位 = Get_乱配列(l_乱数, 60);
            var l_右 = l_単位 + Get_乱配列(l_乱数, 200);
            var l_本当 = l_左 + l_単位 + l_右;
            var l_リード群 = Enumerable.Range(0, 12).Select(i => l_本当.Substring(100 + (i * 5), 200)).ToList();
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);
            Assert.Null(Scaffolder.Get_リードで埋めた繋ぎ目(l_索引, new StringBuilder(l_左), l_右, 200, 0));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 指定した乱数系列の塩基配列
        /// </summary>
        /// <param name="p_乱数"></param>
        /// <param name="p_長さ"></param>
        /// <returns></returns>
        private static string Get_乱配列(Random p_乱数, int p_長さ)
        {
            return new string([..Enumerable.Range(0, p_長さ).Select(_ => Consts.塩基文字[p_乱数.Next(4)])]);
        }

        #endregion
    }
}
