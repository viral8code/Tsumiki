using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// リード長からの k 長自動選択の検証。
    ///
    /// 既定の k=31 は 150bp リードに対して明確に短すぎ、実データで
    /// unitig の N50 が k=63 の場合の 1/5 にしかならなかった。
    /// リード長は 75bp から 300bp まで大きく変わるため、固定値ではなく
    /// 実際のリード長から決める。
    /// </summary>
    public class KmerLengthSelectorTests
    {
        [Theory]
        [InlineData(150, 63)] // 現在の標準。0.6倍は上限を超えるので頭打ち。
        [InlineData(250, 63)] // MiSeq。同じく頭打ち。
        [InlineData(100, 59)] // 0.6倍=60。偶数なので1つ落とす。
        [InlineData(75, 45)]
        [InlineData(50, 29)]  // 0.6倍=30。偶数なので1つ落とす。
        [InlineData(32, 19)]
        public void SuggestKmerLength_ScalesWithReadLength_AndIsCappedAtTheFastPathLimit(int readLength, int expected)
        {
            Assert.Equal(expected, KmerLengthSelector.Get_推奨k長(readLength));
        }

        /// <summary>
        /// k が偶数だと k-mer 自身がその逆相補と一致しうる(回文)ため、
        /// 正規形が縮退して隣接判定が壊れる。どのリード長でも奇数を返すこと。
        /// </summary>
        [Fact]
        public void SuggestKmerLength_IsAlwaysOddAndShorterThanTheRead()
        {
            for (var readLength = Consts.自動k長に必要な最小リード長; readLength <= 400; readLength++)
            {
                var suggestion = KmerLengthSelector.Get_推奨k長(readLength);
                Assert.NotNull(suggestion);
                Assert.True(suggestion.Value % 2 == 1, $"k={suggestion} for read length {readLength} is even");
                Assert.True(suggestion.Value < readLength, $"k={suggestion} is not shorter than read length {readLength}");
                Assert.True(suggestion.Value <= Consts.自動k長の上限);
            }
        }

        [Fact]
        public void SuggestKmerLength_ReadTooShort_ReturnsNull()
        {
            Assert.Null(KmerLengthSelector.Get_推奨k長(Consts.自動k長に必要な最小リード長 - 1));
        }

        [Fact]
        public void Resolve_WhenKmerLengthWasNotGiven_AppliesTheSuggestion()
        {
            var param = new Parameters();
            Assert.False(param.A_k長が明示指定されたか);

            KmerLengthSelector.V_解決_k長(param, 150);

            Assert.Equal(63, param.A_k長);
            // 自動適用は「明示指定された」扱いにしない。
            Assert.False(param.A_k長が明示指定されたか);
        }

        /// <summary>
        /// 明示指定はユーザーの判断なので、推定値で上書きしてはいけない。
        /// </summary>
        [Fact]
        public void Resolve_WhenKmerLengthWasGivenExplicitly_LeavesItAlone()
        {
            var param = new Parameters { A_k長 = 31 };
            Assert.True(param.A_k長が明示指定されたか);

            KmerLengthSelector.V_解決_k長(param, 150);

            Assert.Equal(31, param.A_k長);
        }

        [Fact]
        public void Resolve_WhenReadLengthIsUnknown_KeepsTheDefault()
        {
            var param = new Parameters();

            KmerLengthSelector.V_解決_k長(param, null);

            Assert.Equal(Consts.k長の既定値, param.A_k長);
        }

        [Fact]
        public void Resolve_WhenReadsAreTooShortToPickK_KeepsTheDefault()
        {
            var param = new Parameters();

            KmerLengthSelector.V_解決_k長(param, Consts.自動k長に必要な最小リード長 - 1);

            Assert.Equal(Consts.k長の既定値, param.A_k長);
        }
    }
}
