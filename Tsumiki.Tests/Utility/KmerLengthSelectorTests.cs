using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// リード長からの k 長自動選択の検証
    /// </summary>
    /// <remarks>
    /// 既定の k=31 は 150 bp リードに対して明確に短すぎ、実データで
    /// unitig の N50 が k=63 の場合の 1/5 にしかならなかった<br/>
    /// リード長は 75 bp から 300 bp まで大きく変わるため、固定値ではなく
    /// 実際のリード長から決める
    /// </remarks>
    public class KmerLengthSelectorTests
    {
        /// <summary>
        /// 推奨 k 長はリード長に比例して大きくなり、上限で頭打ちになる
        /// </summary>
        [Theory]
        [InlineData(150, 63)] // 現在の標準、0.6 倍は上限を超えるので頭打ち
        [InlineData(250, 63)] // MiSeq、同じく頭打ち
        [InlineData(100, 59)] // 0.6 倍 = 60、偶数なので 1 つ落とす
        [InlineData(75, 45)]
        [InlineData(50, 29)]  // 0.6 倍 = 30、偶数なので 1 つ落とす
        [InlineData(32, 19)]
        public void 推奨k長はリード長に比例し上限で頭打ちになる(int p_readLength, int p_expected)
        {
            Assert.Equal(p_expected, KmerLengthSelector.Get_推奨k長(p_readLength));
        }

        /// <summary>
        /// k が偶数だと k-mer 自身がその逆相補と一致しうる (回文) ため、
        /// 正規形が縮退して隣接判定が壊れる
        /// </summary>
        /// <remarks>
        /// どのリード長でも奇数を返すこと
        /// </remarks>
        [Fact]
        public void 推奨k長は常に奇数でリード長より短い()
        {
            for (var l_readLength = Consts.自動k長に必要な最小リード長; l_readLength <= 400; l_readLength++)
            {
                var l_suggestion = KmerLengthSelector.Get_推奨k長(l_readLength);
                Assert.NotNull(l_suggestion);
                Assert.True(l_suggestion.Value % 2 == 1, $"k={l_suggestion} for read length {l_readLength} is even");
                Assert.True(l_suggestion.Value < l_readLength, $"k={l_suggestion} is not shorter than read length {l_readLength}");
                Assert.True(l_suggestion.Value <= Consts.自動k長の上限);
            }
        }

        /// <summary>
        /// リードが短すぎる場合は null を返す
        /// </summary>
        [Fact]
        public void リードが短すぎる場合はnullを返す()
        {
            Assert.Null(KmerLengthSelector.Get_推奨k長(Consts.自動k長に必要な最小リード長 - 1));
        }

        /// <summary>
        /// k 長が未指定なら、推奨値を適用する
        /// </summary>
        [Fact]
        public void k長未指定なら推奨値を適用する()
        {
            var l_param = new Parameters();
            Assert.False(l_param.A_k長が明示指定されたか);

            KmerLengthSelector.V_解決_k長(l_param, 150);

            Assert.Equal(63, l_param.A_k長);
            // 自動適用は「明示指定された」扱いにしない
            Assert.False(l_param.A_k長が明示指定されたか);
        }

        /// <summary>
        /// 明示指定はユーザーの判断なので、推定値で上書きしてはいけない
        /// </summary>
        [Fact]
        public void 明示指定されたk長はそのまま残す()
        {
            var l_param = new Parameters { A_k長 = 31 };
            Assert.True(l_param.A_k長が明示指定されたか);

            KmerLengthSelector.V_解決_k長(l_param, 150);

            Assert.Equal(31, l_param.A_k長);
        }

        /// <summary>
        /// リード長が不明なら、既定値を維持する
        /// </summary>
        [Fact]
        public void リード長が不明なら既定値を維持する()
        {
            var l_param = new Parameters();

            KmerLengthSelector.V_解決_k長(l_param, null);

            Assert.Equal(Consts.k長の既定値, l_param.A_k長);
        }

        /// <summary>
        /// リードが短すぎて k 長を選べない場合は、既定値を維持する
        /// </summary>
        [Fact]
        public void リードが短すぎてk長を選べない場合は既定値を維持する()
        {
            var l_param = new Parameters();

            KmerLengthSelector.V_解決_k長(l_param, Consts.自動k長に必要な最小リード長 - 1);

            Assert.Equal(Consts.k長の既定値, l_param.A_k長);
        }
    }
}
