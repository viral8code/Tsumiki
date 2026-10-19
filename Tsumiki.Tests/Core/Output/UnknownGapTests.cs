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
        /// 印の左右の片が長いとき、左の片の末尾 30 bp と右の片の先頭 30 bp を削り、N 100 個に置き換える
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_両側の片の端を削って置き換える()
        {
            var l_乱数 = new Random(41);
            var l_X = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_X + Consts.未確認の繋ぎ目 + l_Y, l_番号群);
            Assert.Equal(l_X[..170] + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + l_Y[30..], l_結果);
            Assert.Equal([0], l_番号群.Order());
        }

        /// <summary>
        /// 左の片の手前に別の N の連続があるとき、左の片はその N より後ろだけで数える (80 bp は 61 以上なので削る)
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_左の片は手前のギャップより後ろだけを数える()
        {
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(new string('A', 100) + "NNNN" + new string('C', 80) + "n" + new string('G', 200), l_番号群);
            Assert.Equal(new string('A', 100) + "NNNN" + new string('C', 50) + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + new string('G', 170), l_結果);
            Assert.Equal([1], l_番号群.Order());
        }

        /// <summary>
        /// 右の片の先頭が左の片の末尾より手前の配列の写しなら、その写しを読み飛ばす (写しの後ろの右の片は残す)
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_右の片の先頭が左の片の写しなら読み飛ばす()
        {
            var l_乱数 = new Random(51);
            var l_X = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            if (l_Y[0] == l_X[140])
            {
                l_Y = (l_X[140] == 'A' ? "C" : "A") + l_Y[1..];
            }

            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_X + Consts.未確認の繋ぎ目 + l_X[60..140] + l_Y, l_番号群);
            Assert.Equal(l_X[..170] + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + l_Y, l_結果);
            Assert.Equal([0], l_番号群.Order());
        }

        /// <summary>
        /// 左の片の末尾が右の片の先頭より後ろの配列の写しなら、その写しを左の片の末尾から除く (左の片の残りは残す)
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_左の片の末尾が右の片の写しなら除く()
        {
            var l_乱数 = new Random(52);
            var l_Z = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            if (l_Z[199] == l_Y[49])
            {
                l_Z = l_Z[..199] + (l_Y[49] == 'A' ? "C" : "A");
            }

            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_Z + l_Y[50..130] + Consts.未確認の繋ぎ目 + l_Y, l_番号群);
            Assert.Equal(l_Z + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + l_Y[30..], l_結果);
            Assert.Equal([0], l_番号群.Order());
        }

        /// <summary>
        /// 右の片の先頭の種が左の片の中に 2 か所あるときは、どちらの写しとも決めずに読み飛ばさない
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_種が左の片に2か所あれば読み飛ばさない()
        {
            var l_乱数 = new Random(53);
            var l_X = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            if (l_Y[0] == l_X[140])
            {
                l_Y = (l_X[140] == 'A' ? "C" : "A") + l_Y[1..];
            }

            var l_X写し = l_X[..10] + l_X[90..115] + l_X[35..];
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_X写し + Consts.未確認の繋ぎ目 + l_X写し[60..140] + l_Y, l_番号群);
            Assert.Equal(l_X写し[..170] + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + (l_X写し[60..140] + l_Y)[30..], l_結果);
            Assert.Equal([0], l_番号群.Order());
        }

        /// <summary>
        /// 右の片の先頭が左の片の末尾の写しで、1 塩基だけ違うなら、違いを許して左の片の末尾まで届いたところで N を書かずに畳む
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_1塩基違う写しは畳む()
        {
            var l_乱数 = new Random(61);
            var l_X = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            var l_X写し = Get_一部を変えた配列(l_X, 160);
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_X + Consts.未確認の繋ぎ目 + l_X写し[100..] + l_Y, l_番号群);
            Assert.Equal(l_X[..150] + l_X写し[150..200] + l_Y, l_結果);
            Assert.Empty(l_番号群);
        }

        /// <summary>
        /// 畳んだ印の後ろの N の連続の番号は、出力の中で何番目の N の連続かに合わせる (畳んだ印は数えない)
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_畳んだ後ろのNの連続の番号は出力に合わせる()
        {
            var l_乱数 = new Random(62);
            var l_X = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            var l_V = Get_乱配列(l_乱数, 200);
            var l_U = Get_乱配列(l_乱数, 200);
            var l_X写し = Get_一部を変えた配列(l_X, 160);
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_X + Consts.未確認の繋ぎ目 + l_X写し[100..] + l_Y + "NNNNN" + l_V + Consts.未確認の繋ぎ目 + l_U, l_番号群);
            Assert.Equal(l_X[..150] + l_X写し[150..200] + l_Y + "NNNNN" + l_V[..170] + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + l_U[30..], l_結果);
            Assert.Equal([1], l_番号群.Order());
        }

        /// <summary>
        /// 写しの違いが 3 か所あれば、違いを許しきれず左の片の末尾まで届かないので畳まず、右の片の先頭 27 塩基を削って N 100 個にする<br/>
        /// 手計算: X と X'' は位置 157・162・166 の 3 か所だけ違う (種 X[130..155) の外)。130〜156 の 27 塩基が一致する。157 で違い、その後の 10 塩基 (158〜167) に 162 の違いが入るので許さない。ℓ = 27 で、左の片の末尾 (170) に届かないので畳まない
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_違いが3か所あれば畳まない()
        {
            var l_乱数 = new Random(63);
            var l_X = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            var l_X写し = Get_一部を変えた配列(l_X, 157, 162, 166);
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_X + Consts.未確認の繋ぎ目 + l_X写し[100..] + l_Y, l_番号群);
            Assert.Equal(l_X[..170] + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + l_X写し[157..] + l_Y, l_結果);
            Assert.Equal([0], l_番号群.Order());
        }

        /// <summary>
        /// 右の片の先頭の種が右の片の中に 2 か所あるとき (縦列反復) は、左の片の末尾から削らず、左の片の末尾の種は右の片に無いので (b) でも削らない
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_右の片の種が縦列反復なら読み飛ばさない()
        {
            var l_乱数 = new Random(71);
            var l_X = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_X + Consts.未確認の繋ぎ目 + l_X[60..140] + l_X[90..115] + l_Y, l_番号群);
            Assert.Equal(l_X[..170] + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + l_X[90..140] + l_X[90..115] + l_Y, l_結果);
            Assert.Equal([0], l_番号群.Order());
        }

        /// <summary>
        /// 左の片の末尾の種が左の片の中に 2 か所あるとき (縦列反復) は、左の片の末尾から削らず、右の片の先頭の種は左の片に無いので (a) でも削らない
        /// </summary>
        [Fact]
        public void Get_長さ不明のギャップへ置換_左の片の種が縦列反復なら除かない()
        {
            var l_乱数 = new Random(72);
            var l_Z = Get_乱配列(l_乱数, 200);
            var l_Y = Get_乱配列(l_乱数, 200);
            var l_Z写し = l_Z[..100] + l_Y[75..100] + l_Z[125..];
            var l_番号群 = new HashSet<int>();
            var l_結果 = FinalAssemblyPipeline.Get_長さ不明のギャップへ置換(l_Z写し + l_Y[50..130] + Consts.未確認の繋ぎ目 + l_Y, l_番号群);
            Assert.Equal(l_Z写し + l_Y[50..100] + new string('N', FinalAssemblyPipeline.C_長さ不明のギャップ長) + l_Y[30..], l_結果);
            Assert.Equal([0], l_番号群.Order());
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
            var l_一致で畳んだ数 = 0;
            var l_結果 = FinalAssemblyPipeline.Get_確かめた繋ぎ目を畳んだ配列(l_左 + Consts.未確認の繋ぎ目 + l_右, l_索引, 89, 100, 0, ref l_総数, ref l_畳んだ数, ref l_埋めた数, ref l_一致で畳んだ数);
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
            var l_一致で畳んだ数 = 0;
            var l_結果 = FinalAssemblyPipeline.Get_確かめた繋ぎ目を畳んだ配列(l_左 + Consts.未確認の繋ぎ目 + l_右, l_索引, 63, 150, 0, ref l_総数, ref l_畳んだ数, ref l_埋めた数, ref l_一致で畳んだ数);
            Assert.Equal((1, 0, 1), (l_総数, l_畳んだ数, l_埋めた数));
            Assert.Equal(l_本当, l_結果);
        }

        /// <summary>
        /// 両側の片が長く完全に一致し、その長さがただ 1 通りなら、リードが無くても畳む
        /// </summary>
        /// <param name="p_共通">両側の片の末尾と先頭が一致する長さ</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_Is畳む">畳むはずか</param>
        [Theory]
        [InlineData(120, 135, true)]
        [InlineData(40, 135, false)]
        [InlineData(120, 101, false)]
        public void Get_確かめた繋ぎ目を畳んだ配列_両側の長い一致で畳む(int p_共通, int p_k長, bool p_Is畳む)
        {
            var l_乱数 = new Random(31);
            var l_共通 = Get_乱配列(l_乱数, p_共通);
            var l_左 = Get_乱配列(l_乱数, 300) + l_共通;
            var l_右 = l_共通 + Get_乱配列(l_乱数, 300);
            var l_索引 = ReadMinimizerIndex.V_構築(() => new List<string>());
            var l_総数 = 0;
            var l_畳んだ数 = 0;
            var l_埋めた数 = 0;
            var l_一致で畳んだ数 = 0;
            var l_結果 = FinalAssemblyPipeline.Get_確かめた繋ぎ目を畳んだ配列(l_左 + Consts.未確認の繋ぎ目 + l_右, l_索引, p_k長, 100, 0, ref l_総数, ref l_畳んだ数, ref l_埋めた数, ref l_一致で畳んだ数);
            Assert.Equal(1, l_総数);
            Assert.Equal(0, l_畳んだ数);
            Assert.Equal(0, l_埋めた数);
            Assert.Equal(p_Is畳む ? 1 : 0, l_一致で畳んだ数);
            Assert.Equal(p_Is畳む ? l_左 + l_右[p_共通..] : l_左 + Consts.未確認の繋ぎ目 + l_右, l_結果);
        }

        /// <summary>
        /// 一致する長さが複数あれば (縦の反復の中)、畳まない
        /// </summary>
        [Fact]
        public void Get_確かめた繋ぎ目を畳んだ配列_一致が複数あれば畳まない()
        {
            var l_乱数 = new Random(32);
            var l_単位 = Get_乱配列(l_乱数, 12);
            var l_共通 = string.Concat(Enumerable.Repeat(l_単位, 10));
            var l_左 = Get_乱配列(l_乱数, 300) + l_共通;
            var l_右 = l_共通 + Get_乱配列(l_乱数, 300);
            var l_索引 = ReadMinimizerIndex.V_構築(() => new List<string>());
            var l_総数 = 0;
            var l_畳んだ数 = 0;
            var l_埋めた数 = 0;
            var l_一致で畳んだ数 = 0;
            var l_結果 = FinalAssemblyPipeline.Get_確かめた繋ぎ目を畳んだ配列(l_左 + Consts.未確認の繋ぎ目 + l_右, l_索引, 135, 100, 0, ref l_総数, ref l_畳んだ数, ref l_埋めた数, ref l_一致で畳んだ数);
            Assert.Equal(0, l_一致で畳んだ数);
            Assert.Equal(l_左 + Consts.未確認の繋ぎ目 + l_右, l_結果);
        }

        /// <summary>
        /// 左の片は、最後のギャップ文字より後ろだけを見る
        /// </summary>
        [Fact]
        public void Get_配列の一致で決めた重なり長_最後のギャップより後ろだけを左の片とする()
        {
            var l_乱数 = new Random(33);
            var l_共通 = Get_乱配列(l_乱数, 120);
            var l_左 = Get_乱配列(l_乱数, 300) + l_共通;
            var l_右 = l_共通 + Get_乱配列(l_乱数, 300);
            var l_出力 = new StringBuilder("ACGT" + 'N' + Get_乱配列(l_乱数, 50) + l_左);
            Assert.Equal(120, Scaffolder.Get_配列の一致で決めた重なり長(l_出力, l_右, 135));
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
        /// 指定した位置の塩基だけを、別の塩基 (A なら C、それ以外なら A) に変えた配列
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_位置群">変える位置</param>
        /// <returns>変えた配列</returns>
        private static string Get_一部を変えた配列(string p_配列, params int[] p_位置群)
        {
            var l_文字群 = p_配列.ToCharArray();
            foreach (var l_位置 in p_位置群)
            {
                l_文字群[l_位置] = l_文字群[l_位置] == 'A' ? 'C' : 'A';
            }

            return new string(l_文字群);
        }

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
