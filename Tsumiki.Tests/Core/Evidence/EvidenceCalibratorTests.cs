using Tsumiki.Core.Evidence;
using Tsumiki.Core;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペア証拠を生カウントではなく期待本数との比で測るための較正器の検証。
    ///
    /// Scaffolder・ContigMaker(分岐選択)・BeamSearchExtender(先読みスコア)が
    /// 同じ較正を共有するために括り出したもの。ここでは較正器そのものの
    /// フォールバック規則と、「短い辺には厳しく、長い辺には緩く」という
    /// 固定閾値の偏りが正規化で解消されることを検証する。
    /// </summary>
    public class EvidenceCalibratorTests
    {
        /// <summary>
        /// 中央 400 付近に集まるフラグメント長の標本。
        /// </summary>
        private static List<int> Get_同一ユニティグ標本(int p_件数 = 2000)
        {
            var l_乱数 = new Random(20260908);
            return [.. Enumerable.Range(0, p_件数).Select(_ => 400 + l_乱数.Next(-50, 51))];
        }

        [Fact]
        public void Get_較正器_IsUnusable_WhenReadLengthIsUnknown()
        {
            var l_較正器 = 証拠較正器.Get_較正器(Get_同一ユニティグ標本(), p_リード長: null, [1000L, 2000L]);

            Assert.False(l_較正器.A_使えるか);
        }

        [Fact]
        public void Get_較正器_IsUnusable_WhenThereAreNoSamples()
        {
            var l_較正器 = 証拠較正器.Get_較正器([], p_リード長: 100, [1000L, 2000L]);

            Assert.False(l_較正器.A_使えるか);
        }

        [Fact]
        public void Get_較正器_IsUnusable_WhenEveryUnitigIsShorterThanTheFragment()
        {
            // すべての unitig がフラグメントより短いと、期待位置数の合計が 0 になり
            // 密度を較正できない(0 除算を避けて安全にフォールバックする)。
            var l_較正器 = 証拠較正器.Get_較正器(Get_同一ユニティグ標本(), p_リード長: 100, [50L, 80L]);

            Assert.False(l_較正器.A_使えるか);
        }

        [Fact]
        public void Get_正規化済み支持_ReturnsZero_WhenTheCalibratorIsUnusable()
        {
            var l_較正器 = 証拠較正器.Get_較正器([], p_リード長: 100, [1000L, 2000L]);

            Assert.Equal(0, l_較正器.Get_正規化済み支持(p_観測本数: 10, p_長さ1: 1000, p_長さ2: 1000, p_ギャップ長: 0));
        }

        /// <summary>
        /// 較正器の核心: 同じ観測本数でも、短い unitig への辺は正規化後の値が
        /// 大きくなる(=期待が小さいので少ない観測でも強い支持とみなされる)。
        /// これが「固定閾値10は短い辺には厳しく、長い辺には緩すぎる」という
        /// 提案Dの問題意識そのものへの解答になっている。
        /// </summary>
        [Fact]
        public void Get_正規化済み支持_IsHigherForTheSameRawCount_WhenTheFlankingUnitigIsShorter()
        {
            // 分岐元(片側)の長さは固定し、行き先側の長さだけを短い/長いで変える
            // (現実の分岐選択でも、変わるのは行き先の unitig 長のほうである)。
            var l_標本 = Get_同一ユニティグ標本();
            var l_較正器 = 証拠較正器.Get_較正器(l_標本, p_リード長: 100, [100_000L, 150L, 50_000L]);

            var l_短い辺への支持 = l_較正器.Get_正規化済み支持(p_観測本数: 8, p_長さ1: 100_000, p_長さ2: 150, p_ギャップ長: 0);
            var l_長い辺への支持 = l_較正器.Get_正規化済み支持(p_観測本数: 8, p_長さ1: 100_000, p_長さ2: 50_000, p_ギャップ長: 0);

            Assert.True(l_較正器.A_使えるか);
            Assert.True(l_短い辺への支持 > l_長い辺への支持,
                $"expected short-flank ratio ({l_短い辺への支持}) to exceed long-flank ratio ({l_長い辺への支持}) for the same raw count");
        }

        /// <summary>
        /// 観測本数がちょうど期待本数どおりなら比はおよそ 1.0 になること。
        /// </summary>
        [Fact]
        public void Get_正規化済み支持_IsAboutOne_WhenObservedMatchesIdeal()
        {
            var l_標本 = Get_同一ユニティグ標本(20_000);
            IReadOnlyList<long> l_ユニティグ長一覧 = [50_000L, 50_000L, 50_000L];
            var l_較正器 = 証拠較正器.Get_較正器(l_標本, p_リード長: 100, l_ユニティグ長一覧);

            // 密度較正に使ったのと同じ長さの unitig 同士の辺なら、
            // 「観測本数 = 密度 x 期待位置数」を代入すれば比はちょうど1になる。
            var l_モデル = new PairedDistanceModel(l_標本, p_リード長: 100);
            var l_期待位置数 = l_モデル.Get_期待位置数(50_000, 50_000, 0);
            var l_密度相当の観測本数 = (ulong)Math.Round(l_標本.Count / (3 * l_モデル.Get_期待位置数_単一(50_000)) * l_期待位置数);

            var l_比 = l_較正器.Get_正規化済み支持(l_密度相当の観測本数, 50_000, 50_000, 0);

            Assert.InRange(l_比, 0.9, 1.1);
        }
    }
}
