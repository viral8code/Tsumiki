using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    public class PhredSnifferTests
    {
        [Fact]
        public void Sample_ComputesMinMaxAsciiAcrossAllLines()
        {
            var sample = PhredSniffer.Get_標本(["hhhh", "IIhh", "!!!!"]);

            Assert.Equal('!', (char)sample.A_最小ASCII);
            Assert.Equal('h', (char)sample.A_最大ASCII);
            Assert.Equal(3, sample.A_標本リード数);
            Assert.Equal(12, sample.A_標本文字数);
        }

        [Fact]
        public void Sample_StopsAtMaxReadsToSample()
        {
            var lines = Enumerable.Repeat("hhhh", 100);

            var sample = PhredSniffer.Get_標本(lines, p_標本上限: 5);

            Assert.Equal(5, sample.A_標本リード数);
        }

        [Fact]
        public void IsUniform_TrueWhenEveryCharIsIdentical()
        {
            var sample = PhredSniffer.Get_標本(["hhhh", "hhhhhh"]);

            Assert.True(sample.A_一様か);
        }

        [Fact]
        public void IsUniform_FalseWhenQualityVaries()
        {
            var sample = PhredSniffer.Get_標本(["hhIh", "hhhh"]);

            Assert.False(sample.A_一様か);
        }

        [Fact]
        public void BuildWarning_RealisticPhred33Data_NoWarning()
        {
            // 典型的なPhred33品質文字('#'=Q2 〜 'J'=Q41相当)を模した、
            // ばらつきのあるサンプル。
            var sample = PhredSniffer.Get_標本(["#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"]);

            var warning = PhredSniffer.Get_警告文(sample, p_有効オフセット: 33);

            Assert.Null(warning);
        }

        [Fact]
        public void BuildWarning_Phred33WithAsciiImplausiblyHigh_WarnsToTryPhred64()
        {
            // 'h' = ASCII 104。Phred33なら Q=71 となり非現実的。
            var sample = PhredSniffer.Get_標本(["hhhh"]);

            var warning = PhredSniffer.Get_警告文(sample, p_有効オフセット: 33);

            Assert.NotNull(warning);
            Assert.Contains("Phred64", warning);
        }

        [Fact]
        public void BuildWarning_Phred64WithSameData_DecodesToPlausibleQ_ButStillFlagsUniformity()
        {
            // 'h' = ASCII 104。Phred64なら Q=40 で妥当な範囲だが、
            // 全く同一の値しか出ていない点は別途警告する。
            var sample = PhredSniffer.Get_標本(["hhhh", "hhhh", "hhhh"]);

            var warning = PhredSniffer.Get_警告文(sample, p_有効オフセット: 64);

            Assert.NotNull(warning);
            Assert.Contains("uniform", warning);
            Assert.DoesNotContain("Phred33", warning);
        }

        [Fact]
        public void BuildWarning_NegativeQ_WarnsRegardlessOfUniformity()
        {
            // '!' = ASCII 33。Phred64なら Q=-31 となり明らかに不正。
            var sample = PhredSniffer.Get_標本(["!!!!"]);

            var warning = PhredSniffer.Get_警告文(sample, p_有効オフセット: 64);

            Assert.NotNull(warning);
            Assert.Contains("Phred33", warning);
        }

        [Fact]
        public void BuildWarning_EmptySample_ReturnsNull()
        {
            var sample = PhredSniffer.Get_標本([]);

            Assert.Null(PhredSniffer.Get_警告文(sample, 33));
        }

        /// <summary>
        /// 実データ(Achromobacter の IS350 ライブラリ)で観測された ASCII 範囲
        /// [64, 104]。Phred33 と解釈すると Q[31, 71] となり上限がありえないが、
        /// Phred64 なら Q[0, 40] で完全に妥当。この判別ができないと、
        /// 「quality - Phred - QualityCutoff が負なら捨てる」という品質フィルタが
        /// 事実上まったく効かなくなる(Q0 の塩基が Q31 に見えるため)。
        /// </summary>
        [Fact]
        public void InferOffset_RealWorldPhred64Range_InfersPhred64()
        {
            var sample = PhredSniffer.Get_標本([new string((char)64, 4) + new string((char)104, 4)]);

            Assert.Equal(64, PhredSniffer.Get_推定オフセット(sample));
        }

        [Fact]
        public void InferOffset_TypicalPhred33Range_InfersPhred33()
        {
            // '!'(ASCII 33, Q0)から 'I'(ASCII 73, Q40)までの一般的な Phred33 範囲。
            // Phred64 と解釈すると Q が負になるため、33 side のみが妥当。
            var sample = PhredSniffer.Get_標本(["!!!!IIII"]);

            Assert.Equal(33, PhredSniffer.Get_推定オフセット(sample));
        }

        [Fact]
        public void InferOffset_AmbiguousRange_ReturnsNull()
        {
            // ASCII 66-70 は Phred33 なら Q[33,37]、Phred64 なら Q[2,6]。
            // どちらの解釈でも現実的な範囲に収まるため判別できない。
            var sample = PhredSniffer.Get_標本(["BCDEF"]);

            Assert.Null(PhredSniffer.Get_推定オフセット(sample));
        }

        [Fact]
        public void InferOffset_EmptySample_ReturnsNull()
        {
            Assert.Null(PhredSniffer.Get_推定オフセット(PhredSniffer.Get_標本([])));
        }
    }
}
