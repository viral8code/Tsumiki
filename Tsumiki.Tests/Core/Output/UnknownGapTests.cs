using Tsumiki.Commons;
using Tsumiki.Cores.Pipeline;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終成果物で、未確認の繋ぎ目を長さ不明のギャップにし、AGP に書き分ける処理の検証
    /// </summary>
    public class UnknownGapTests
    {
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
            Assert.Equal(p_配列.Replace("N", string.Empty).Replace("n", string.Empty), l_結果.Replace("N", string.Empty));
        }

        /// <summary>
        /// AGP では、長さ不明のギャップを U、推定長のギャップを N、配列の片を W として書く
        /// </summary>
        [Fact]
        public void Get_AGP行_ギャップの種類を書き分ける()
        {
            var l_番号群 = new HashSet<int>();
            var l_配列 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換("ACGTnACNNNGT", l_番号群);

            var l_行群 = FinalAssemblyPipeline.Get_AGP行("S1", l_配列, l_番号群).Select(x => x.Split('\t')).ToList();

            Assert.Equal(["W", "U", "W", "N", "W"], l_行群.Select(x => x[4]));
            Assert.Equal("100", l_行群[1][5]);
            Assert.Equal("3", l_行群[3][5]);
            Assert.Equal(["S1_1", "S1_2", "S1_3"], l_行群.Where(x => x[4] == "W").Select(x => x[5]));
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
            string Get_乱配列(int p_長さ) => new([.. Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)])]);
            var l_共通 = Get_乱配列(30);
            var l_左 = Get_乱配列(200) + l_共通;
            var l_右 = l_共通 + Get_乱配列(200);
            var l_繋いだ配列 = l_左 + l_右[30..];
            var l_リード群 = Enumerable.Range(0, p_リード数).Select(i => l_繋いだ配列.Substring(170 + i, 100)).ToList();
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);
            var l_総数 = 0;
            var l_畳んだ数 = 0;

            var l_結果 = FinalAssemblyPipeline.Get_確かめた繋ぎ目を畳んだ配列(l_左 + Consts.未確認の繋ぎ目 + l_右, l_索引, 89, 100, ref l_総数, ref l_畳んだ数);

            Assert.Equal(1, l_総数);
            Assert.Equal(p_Is畳む ? 1 : 0, l_畳んだ数);
            Assert.Equal(p_Is畳む ? l_繋いだ配列 : l_左 + Consts.未確認の繋ぎ目 + l_右, l_結果);
        }

        #endregion
    }
}
