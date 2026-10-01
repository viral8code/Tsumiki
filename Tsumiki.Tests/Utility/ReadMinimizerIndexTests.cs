using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// minimizer 索引の「出てくるか」が、リードとその逆相補を愚直に探した結果と一致することの検証
    /// </summary>
    public class ReadMinimizerIndexTests
    {
        private static string Get_逆相補(string p_配列)
        {
            return new string([.. p_配列.Reverse().Select(static x => x switch { 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', _ => 'N' })]);
        }

        [Theory]
        [InlineData(1, 0.02)]
        [InlineData(2, 0.0)]
        [InlineData(3, 0.05)]
        public void 愚直な照合と一致する(int p_種, double p_誤り率)
        {
            var l_乱数 = new Random(p_種);
            // 反復を含むゲノム (同じ 60 塩基を 3 回入れる)
            var l_反復 = new string([.. Enumerable.Range(0, 60).Select(_ => "ACGT"[l_乱数.Next(4)])]);
            var l_ゲノム = string.Concat(Enumerable.Range(0, 20).Select(i => i % 7 == 3 ? l_反復 : new string([.. Enumerable.Range(0, 100).Select(_ => "ACGT"[l_乱数.Next(4)])])));
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
                        l_文字[j] = l_乱数.Next(10) == 0 ? 'N' : "ACGT"[l_乱数.Next(4)];
                    }
                }
                var l_リード = new string(l_文字);
                l_リード群.Add(l_乱数.Next(2) == 0 ? l_リード : Get_逆相補(l_リード));
            }

            var l_索引 = ReadMinimizerIndex.V_構築(() => l_リード群);
            var l_全文 = string.Join("|", l_リード群.SelectMany(static x => new[] { x, Get_逆相補(x) }));

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
        /// 錨の直前と直後を、リードが逆相補で入っていても錨の向きで返す
        /// </summary>
        [Fact]
        public void Get_前後群_錨の向きで前後を返す()
        {
            var l_乱数 = new Random(5);
            var l_ゲノム = new string([.. Enumerable.Range(0, 300).Select(_ => "ACGT"[l_乱数.Next(4)])]);
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
            var l_単位 = new string([.. Enumerable.Range(0, 60).Select(_ => "ACGT"[l_乱数.Next(4)])]);
            var l_索引 = ReadMinimizerIndex.V_構築(() => Enumerable.Repeat(l_単位, 20_000));

            Assert.Equal(p_期待, l_索引.Get_出現数(l_単位.AsSpan(5, ReadMinimizerIndex.C_最短の問い合わせ長), p_上限));
        }
    }
}
