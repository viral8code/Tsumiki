using Tsumiki.Common;

namespace Tsumiki.Tests.Common
{
    /// <summary>
    /// 曖昧塩基の判定。候補の中身ではなく個数だけが必要な最内ループで使うため、
    /// List を返す版と同じ答えを返すことを固定する。
    /// </summary>
    public class AmbiguousBaseTests
    {
        [Theory]
        [InlineData('A')]
        [InlineData('C')]
        [InlineData('G')]
        [InlineData('T')]
        [InlineData('M')]
        [InlineData('V')]
        [InlineData('N')]
        [InlineData('H')]
        [InlineData('R')]
        [InlineData('D')]
        [InlineData('W')]
        [InlineData('S')]
        [InlineData('B')]
        [InlineData('Y')]
        [InlineData('K')]
        public void Get_曖昧塩基か_AgreesWithTheCandidateCount(char p_塩基文字)
        {
            Assert.Equal(
                Util.Get_塩基ID候補(p_塩基文字).Count > 1,
                Util.Get_曖昧塩基か(p_塩基文字));
        }

        /// <summary>
        /// 不正な文字を黙って通さないこと。List を返す版と同じく例外にする。
        /// </summary>
        [Theory]
        [InlineData('a')]
        [InlineData('Z')]
        [InlineData('-')]
        public void Get_曖昧塩基か_RejectsANonBaseCharacter(char p_文字)
        {
            _ = Assert.Throws<ArgumentException>(() => Util.Get_曖昧塩基か(p_文字));
        }
    }
}
