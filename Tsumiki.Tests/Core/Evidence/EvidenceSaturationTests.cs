using Tsumiki.Cores.Evidence;
using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 証拠の飽和と独立性の数え方を固定する
    /// </summary>
    /// <remarks>
    /// 狙いは、同じ話を何度も聞いたことが別の裏付けを得たことに化けないようにすること
    /// </remarks>
    public class EvidenceSaturationTests
    {
        /// <summary>
        /// 支持が無ければ飽和支持は0になること
        /// </summary>
        [Fact]
        public void Get_飽和支持_支持が無ければ0になる()
        {
            Assert.Equal(0, 証拠較正器.Get_飽和支持(0));
            Assert.Equal(0, 証拠較正器.Get_飽和支持(-5));
        }

        /// <summary>
        /// 飽和支持は本数の増加とともに単調に増えるが1を超えないこと
        /// </summary>
        [Fact]
        public void Get_飽和支持_単調に増えるが1を超えない()
        {
            var l_直前 = 0.0;
            foreach (var l_支持 in new[] { 1, 2, 3, 5, 10, 100 })
            {
                var l_値 = 証拠較正器.Get_飽和支持(l_支持);
                Assert.True(l_値 > l_直前, $"支持 {l_支持} で単調性が崩れた");
                Assert.True(l_値 < 1.0, $"支持 {l_支持} で 1 を超えた");
                l_直前 = l_値;
            }

            // 極端な本数では倍精度の分解能で 1 に到達するが、超えることはない
            Assert.True(証拠較正器.Get_飽和支持(10000) <= 1.0);
        }

        /// <summary>
        /// 飽和支持は本数を増やすほど伸びが小さくなり頭打ちになること
        /// </summary>
        [Fact]
        public void Get_飽和支持_本数を増やしても頭打ちになる()
        {
            // 3 本から 10 本への伸びより、1 本から 3 本への伸びのほうが大きい
            var l_1から3 = 証拠較正器.Get_飽和支持(3) - 証拠較正器.Get_飽和支持(1);
            var l_3から10 = 証拠較正器.Get_飽和支持(10) - 証拠較正器.Get_飽和支持(3);
            Assert.True(l_1から3 > l_3から10);

            // 1000 本と 10000 本はほぼ区別が付かない
            Assert.True(証拠較正器.Get_飽和支持(10000) - 証拠較正器.Get_飽和支持(1000) < 0.001);
        }

        /// <summary>
        /// 同じ距離を示す観測はまとめて1つの独立支持として数えること
        /// </summary>
        [Fact]
        public void Get_独立支持数_同じ距離を示す観測は1つに畳む()
        {
            // 全く同じ距離ばかりの 100 本は、別々の分子から得た裏付けとは言えない
            Assert.Equal(1, 証拠較正器.Get_独立支持数([.. Enumerable.Repeat(300, 100)]));
        }

        /// <summary>
        /// 距離が散らばっていれば独立支持数はその種類数になること
        /// </summary>
        [Fact]
        public void Get_独立支持数_距離が散らばっていればその数だけ数える()
        {
            Assert.Equal(4, 証拠較正器.Get_独立支持数([300, 305, 310, 305, 300, 320]));
        }

        /// <summary>
        /// 観測が無ければ独立支持数は0になること
        /// </summary>
        [Fact]
        public void Get_独立支持数_観測が無ければ0になる()
        {
            Assert.Equal(0, 証拠較正器.Get_独立支持数([]));
        }
    }
}
