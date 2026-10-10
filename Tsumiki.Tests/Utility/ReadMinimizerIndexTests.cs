using Tsumiki.Commons;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// minimizer 索引の「出てくるか」が、リードとその逆相補を愚直に探した結果と一致することの検証
    /// </summary>
    public class ReadMinimizerIndexTests
    {
        #region 定数

        /// <summary>
        /// 表区切り
        /// </summary>
        private const string C_表区切り = "|";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 愚直な照合と一致する
        /// </summary>
        /// <param name="p_種"></param>
        /// <param name="p_誤り率"></param>
        [Theory]
        [InlineData(1, 0.02D)]
        [InlineData(2, 0.0D)]
        [InlineData(3, 0.05D)]
        public void V_愚直な照合と一致する(int p_種, double p_誤り率)
        {
            var l_乱数 = new Random(p_種);
            var l_反復 = new string([..Enumerable.Range(0, 60).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_ゲノム = string.Concat(Enumerable.Range(0, 20).Select(i => i % 7 == 3 ? l_反復 : new string([..Enumerable.Range(0, 100).Select(_ => Consts.塩基文字[l_乱数.Next(4)])])));
            var l_リード群 = new List<string>();
            for (var i = 0; i < 400; i++)
            {
                var l_長さ = l_乱数.Next(20, 151);
                var l_開始 = l_乱数.Next(0, l_ゲノム.Length - l_長さ);
                var l_文字 = l_ゲノム.Substring(l_開始, l_長さ).ToCharArray();
                for (var j = 0; j < l_文字.Length; j++)
                {
                    if (l_乱数.NextDouble() < p_誤り率)
                    {
                        l_文字[j] = l_乱数.Next(10) == 0 ? 'N' : Consts.塩基文字[l_乱数.Next(4)];
                    }
                }

                var l_リード = new string(l_文字);
                l_リード群.Add(l_乱数.Next(2) == 0 ? l_リード : Get_逆相補(l_リード));
            }

            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);
            var l_全文 = string.Join(C_表区切り, l_リード群.SelectMany(static x => new[] { x, Get_逆相補(x) }));
            var l_調べる = new List<string>();
            for (var i = 0; i < 3_000; i++)
            {
                var l_長さ = l_乱数.Next(ReadMinimizerIndex.C_最短の問い合わせ長, 90);
                var l_元 = i % 3 == 0 ? l_ゲノム : l_リード群[l_乱数.Next(l_リード群.Count)];
                if (l_元.Length < l_長さ)
                {
                    continue;
                }

                var l_配列 = l_元.Substring(l_乱数.Next(0, l_元.Length - l_長さ + 1), l_長さ);
                if (!l_配列.Contains('N'))
                {
                    l_調べる.Add(i % 2 == 0 ? l_配列 : Get_逆相補(l_配列));
                }
            }

            foreach (var l_配列 in l_調べる)
            {
                Assert.True(l_全文.Contains(l_配列, StringComparison.Ordinal) == l_索引.Has出現(l_配列), l_配列);
            }
        }

        /// <summary>
        /// 窓の種数 17 の索引 (直接作ったものと、既定の索引から窓違いで作ったもの) も、愚直な照合と一致する (31〜89 塩基の問い合わせ)
        /// </summary>
        /// <param name="p_種"></param>
        /// <param name="p_誤り率"></param>
        [Theory]
        [InlineData(1, 0.02D)]
        [InlineData(2, 0.0D)]
        [InlineData(3, 0.05D)]
        public void V_窓17の索引も愚直な照合と一致する(int p_種, double p_誤り率)
        {
            var (l_リード群, l_調べる, l_全文) = Get_照合データ(p_種, p_誤り率, ReadMinimizerIndex.C_短い問い合わせの最短長);
            var l_窓17 = ReadMinimizerIndex.V_構築(() => l_リード群, 17);
            var l_既定から = ReadMinimizerIndex.V_構築(() => l_リード群).Get_窓違い(17);
            foreach (var l_配列 in l_調べる)
            {
                var l_期待 = l_全文.Contains(l_配列, StringComparison.Ordinal);
                Assert.True(l_期待 == l_窓17.Has出現(l_配列), l_配列);
                Assert.True(l_期待 == l_既定から.Has出現(l_配列), l_配列);
            }
        }

        /// <summary>
        /// 窓違いで作った窓 17 の索引と、窓 17 で直接作った索引は、同じ問い合わせに同じ出現数を返す
        /// </summary>
        [Fact]
        public void Get_窓違い17は窓17で直接作った索引と出現数が一致する()
        {
            var (l_リード群, l_調べる, _) = Get_照合データ(1, 0.02D, ReadMinimizerIndex.C_短い問い合わせの最短長);
            var l_窓17 = ReadMinimizerIndex.V_構築(() => l_リード群, 17);
            var l_既定から = ReadMinimizerIndex.V_構築(() => l_リード群).Get_窓違い(17);
            foreach (var l_配列 in l_調べる)
            {
                Assert.Equal(l_窓17.Get_出現数(l_配列, 100), l_既定から.Get_出現数(l_配列, 100));
            }
        }

        /// <summary>
        /// 窓 17 の索引は、最短の問い合わせ長 (31) に満たない 30 塩基の配列を拒む
        /// </summary>
        [Fact]
        public void V_窓17の索引は31塩基未満の配列を拒む()
        {
            var l_乱数 = new Random(4);
            var l_リード = new string([..Enumerable.Range(0, 100).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_索引 = ReadMinimizerIndex.V_構築(() => [l_リード], 17);
            _ = Assert.Throws<ArgumentException>(() => l_索引.Get_出現数(l_リード.AsSpan(0, 30), 1));
        }

        /// <summary>
        /// 既定の窓の索引は、最短の問い合わせ長 (41) に満たない 31〜40 塩基の配列を、今までどおり拒む
        /// </summary>
        /// <param name="p_長さ">問い合わせの長さ</param>
        [Theory]
        [InlineData(31)]
        [InlineData(35)]
        [InlineData(40)]
        public void V_既定の窓の索引は41塩基未満の配列を拒む(int p_長さ)
        {
            var l_乱数 = new Random(6);
            var l_リード = new string([..Enumerable.Range(0, 100).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_索引 = ReadMinimizerIndex.V_構築(() => [l_リード]);
            _ = Assert.Throws<ArgumentException>(() => l_索引.Get_出現数(l_リード.AsSpan(0, p_長さ), 1));
        }

        /// <summary>
        /// 錨の直前と直後を、リードが逆相補で入っていても錨の向きで返す
        /// </summary>
        [Fact]
        public void Get_前後群_錨の向きで前後を返す()
        {
            var l_乱数 = new Random(5);
            var l_ゲノム = new string([..Enumerable.Range(0, 300).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_索引 = ReadMinimizerIndex.V_構築(() => [l_ゲノム.Substring(50, 150), Get_逆相補(l_ゲノム.Substring(80, 150))]);
            var l_錨 = l_ゲノム.Substring(120, ReadMinimizerIndex.C_最短の問い合わせ長);
            var l_錨の終わり = 120 + ReadMinimizerIndex.C_最短の問い合わせ長;
            var l_前後群 = l_索引.Get_前後群(l_錨, 200, 10).OrderBy(x => x.A_前.Length).ToList();
            Assert.Equal([(l_ゲノム[80..120], l_ゲノム[l_錨の終わり..230]), (l_ゲノム[50..120], l_ゲノム[l_錨の終わり..200])], l_前後群);
        }

        /// <summary>
        /// 多くのリードに出てくる配列も、同じ場所を二重に数えずに上限まで数える
        /// </summary>
        /// <param name="p_上限">数え上げの上限</param>
        /// <param name="p_期待">期待する出現数</param>
        [Theory]
        [InlineData(100_000, 20_000)]
        [InlineData(5_000, 5_000)]
        public void Get_出現数_多く出てくる配列も上限まで重ねずに数える(int p_上限, int p_期待)
        {
            var l_乱数 = new Random(9);
            var l_単位 = new string([..Enumerable.Range(0, 60).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_索引 = ReadMinimizerIndex.V_構築(() => Enumerable.Repeat(l_単位, 20_000));
            Assert.Equal(p_期待, l_索引.Get_出現数(l_単位.AsSpan(5, ReadMinimizerIndex.C_最短の問い合わせ長), p_上限));
        }

        /// <summary>
        /// 全窓の判定は、窓ごとに「窓に N が無ければ Has 出現、あれば false」と一致する (窓 27 と窓 17 の索引、窓長 31〜60 のうち索引の最短以上)
        /// </summary>
        /// <param name="p_窓の種数">索引の窓の種数</param>
        [Theory]
        [InlineData(27)]
        [InlineData(17)]
        public void V_判定_出現_全窓は窓ごとのHas出現と一致する(int p_窓の種数)
        {
            var (l_リード群, l_調べる, _) = Get_照合データ(7, 0.02D, ReadMinimizerIndex.C_短い問い合わせの最短長);
            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群, p_窓の種数);
            var l_乱数 = new Random(11);
            for (var k = 0; k < 60; k++)
            {
                var l_長さ = l_乱数.Next(50, 401);
                var l_連結 = string.Concat(Enumerable.Range(0, 20).Select(_ => l_調べる[l_乱数.Next(l_調べる.Count)]));
                var l_文字 = l_連結.Substring(0, l_長さ).ToCharArray();
                for (var j = 0; j < l_文字.Length; j++)
                {
                    if (l_乱数.Next(40) == 0)
                    {
                        l_文字[j] = 'N';
                    }
                }

                var l_配列 = new string(l_文字);
                for (var l_窓 = Math.Max(31, l_索引.A_最短の問い合わせ長); l_窓 <= 60 && l_窓 <= l_配列.Length; l_窓++)
                {
                    var l_結果 = new bool[l_配列.Length - l_窓 + 1];
                    l_索引.V_判定_出現_全窓(l_配列.AsSpan(), l_窓, l_結果.AsSpan());
                    for (var i = 0; i < l_結果.Length; i++)
                    {
                        var l_窓の配列 = l_配列.AsSpan(i, l_窓);
                        var l_期待 = !l_窓の配列.Contains('N') && l_索引.Has出現(l_窓の配列);
                        Assert.True(l_期待 == l_結果[i], $"窓 {l_窓} の開始 {i}");
                    }
                }
            }
        }

        /// <summary>
        /// 全窓の判定は、窓長が索引の最短の問い合わせ長に満たない場合と、結果の長さが合わない場合に ArgumentException を投げる
        /// </summary>
        [Fact]
        public void V_判定_出現_全窓は窓の短さと結果の長さの違いを拒む()
        {
            var l_乱数 = new Random(12);
            var l_配列 = new string([..Enumerable.Range(0, 100).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_索引 = ReadMinimizerIndex.V_構築(() => [l_配列], 17);
            Assert.Throws<ArgumentException>(() => l_索引.V_判定_出現_全窓(l_配列.AsSpan(), 30, new bool[71]));
            Assert.Throws<ArgumentException>(() => l_索引.V_判定_出現_全窓(l_配列.AsSpan(), 41, new bool[10]));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 逆相補
        /// </summary>
        /// <param name="p_配列"></param>
        /// <returns></returns>
        private static string Get_逆相補(string p_配列)
        {
            return new string([..p_配列.Reverse().Select(static x => x switch
            {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => 'N'
            })]);
        }

        /// <summary>
        /// 愚直な照合と比べるための、リード群・問い合わせ群・全文を作る
        /// </summary>
        /// <param name="p_種">乱数の種</param>
        /// <param name="p_誤り率">リードに入れる誤りの率</param>
        /// <param name="p_最短">問い合わせの最短長</param>
        /// <returns>(リード群, 問い合わせ群, リードとその逆相補をつないだ全文)</returns>
        private static (List<string> A_リード群, List<string> A_調べる, string A_全文) Get_照合データ(int p_種, double p_誤り率, int p_最短)
        {
            var l_乱数 = new Random(p_種);
            var l_反復 = new string([..Enumerable.Range(0, 60).Select(_ => Consts.塩基文字[l_乱数.Next(4)])]);
            var l_ゲノム = string.Concat(Enumerable.Range(0, 20).Select(i => i % 7 == 3 ? l_反復 : new string([..Enumerable.Range(0, 100).Select(_ => Consts.塩基文字[l_乱数.Next(4)])])));
            var l_リード群 = new List<string>();
            for (var i = 0; i < 400; i++)
            {
                var l_長さ = l_乱数.Next(20, 151);
                var l_開始 = l_乱数.Next(0, l_ゲノム.Length - l_長さ);
                var l_文字 = l_ゲノム.Substring(l_開始, l_長さ).ToCharArray();
                for (var j = 0; j < l_文字.Length; j++)
                {
                    if (l_乱数.NextDouble() < p_誤り率)
                    {
                        l_文字[j] = l_乱数.Next(10) == 0 ? 'N' : Consts.塩基文字[l_乱数.Next(4)];
                    }
                }

                var l_リード = new string(l_文字);
                l_リード群.Add(l_乱数.Next(2) == 0 ? l_リード : Get_逆相補(l_リード));
            }

            var l_全文 = string.Join(C_表区切り, l_リード群.SelectMany(static x => new[] { x, Get_逆相補(x) }));
            var l_調べる = new List<string>();
            for (var i = 0; i < 3_000; i++)
            {
                var l_長さ = l_乱数.Next(p_最短, 90);
                var l_元 = i % 3 == 0 ? l_ゲノム : l_リード群[l_乱数.Next(l_リード群.Count)];
                if (l_元.Length < l_長さ)
                {
                    continue;
                }

                var l_配列 = l_元.Substring(l_乱数.Next(0, l_元.Length - l_長さ + 1), l_長さ);
                if (!l_配列.Contains('N'))
                {
                    l_調べる.Add(i % 2 == 0 ? l_配列 : Get_逆相補(l_配列));
                }
            }

            return (l_リード群, l_調べる, l_全文);
        }

        #endregion
    }
}
