using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// クオリティ文字から Phred オフセットを推定する処理の検証
    /// </summary>
    public class PhredSnifferTests
    {
        /// <summary>
        /// 標本の全行を通じた最小・最大 ASCII を求めることを検証する
        /// </summary>
        [Fact]
        public void Sample_全行を通じた最小最大ASCIIを求める()
        {
            var sample = PhredSniffer.Get_標本(["hhhh", "IIhh", "!!!!"]);

            Assert.Equal('!', (char)sample.A_最小ASCII);
            Assert.Equal('h', (char)sample.A_最大ASCII);
            Assert.Equal(3, sample.A_標本リード数);
            Assert.Equal(12, sample.A_標本文字数);
        }

        /// <summary>
        /// 標本抽出が指定した上限リード数で打ち切られることを検証する
        /// </summary>
        [Fact]
        public void Sample_指定した上限リード数で打ち切る()
        {
            var lines = Enumerable.Repeat("hhhh", 100);

            var sample = PhredSniffer.Get_標本(lines, p_標本上限: 5);

            Assert.Equal(5, sample.A_標本リード数);
        }

        /// <summary>
        /// すべての文字が同一なら一様かが真になることを検証する
        /// </summary>
        [Fact]
        public void IsUniform_すべての文字が同一なら真になる()
        {
            var sample = PhredSniffer.Get_標本(["hhhh", "hhhhhh"]);

            Assert.True(sample.A_一様か);
        }

        /// <summary>
        /// クオリティがばらつけば一様かが偽になることを検証する
        /// </summary>
        [Fact]
        public void IsUniform_クオリティがばらつけば偽になる()
        {
            var sample = PhredSniffer.Get_標本(["hhIh", "hhhh"]);

            Assert.False(sample.A_一様か);
        }

        /// <summary>
        /// 現実的な Phred33 データでは警告文が出ないことを検証する
        /// </summary>
        [Fact]
        public void BuildWarning_現実的なPhred33データでは警告が出ない()
        {
            // 典型的な Phred33 品質文字 ('#'=Q2 〜 'J'=Q41 相当) を模した、
            // ばらつきのあるサンプル
            var sample = PhredSniffer.Get_標本(["#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"]);

            var warning = PhredSniffer.Get_警告文(sample, p_有効オフセット: 33);

            Assert.Null(warning);
        }

        /// <summary>
        /// Phred33 で ASCII が非現実的に高いときは Phred64 を試すよう促すことを検証する
        /// </summary>
        [Fact]
        public void BuildWarning_Phred33でASCIIが非現実的に高いときはPhred64を促す()
        {
            // 'h' = ASCII 104
            // Phred33 なら Q=71 となり非現実的
            var sample = PhredSniffer.Get_標本(["hhhh"]);

            var warning = PhredSniffer.Get_警告文(sample, p_有効オフセット: 33);

            Assert.NotNull(warning);
            Assert.Contains("Phred64", warning);
        }

        /// <summary>
        /// Phred64 なら妥当な Q になるが、一様性はなお警告することを検証する
        /// </summary>
        [Fact]
        public void BuildWarning_Phred64なら妥当なQになるが一様性はなお警告する()
        {
            // 'h' = ASCII 104
            // Phred64 なら Q=40 で妥当な範囲だが、
            // 全く同一の値しか出ていない点は別途警告する
            var sample = PhredSniffer.Get_標本(["hhhh", "hhhh", "hhhh"]);

            var warning = PhredSniffer.Get_警告文(sample, p_有効オフセット: 64);

            Assert.NotNull(warning);
            Assert.Contains("uniform", warning);
            Assert.DoesNotContain("Phred33", warning);
        }

        /// <summary>
        /// 負の Q になる場合は、一様性に関わらず警告することを検証する
        /// </summary>
        [Fact]
        public void BuildWarning_負のQになる場合は一様性に関わらず警告する()
        {
            // '!' = ASCII 33
            // Phred64 なら Q=-31 となり明らかに不正
            var sample = PhredSniffer.Get_標本(["!!!!"]);

            var warning = PhredSniffer.Get_警告文(sample, p_有効オフセット: 64);

            Assert.NotNull(warning);
            Assert.Contains("Phred33", warning);
        }

        /// <summary>
        /// 空の標本では警告文が null になることを検証する
        /// </summary>
        [Fact]
        public void BuildWarning_空の標本ではnullを返す()
        {
            var sample = PhredSniffer.Get_標本([]);

            Assert.Null(PhredSniffer.Get_警告文(sample, 33));
        }

        /// <summary>
        /// 実データ (Achromobacter の IS350 ライブラリ) で観測された ASCII 範囲
        /// [64, 104]
        /// </summary>
        /// <remarks>
        /// Phred33 と解釈すると Q[31, 71] となり上限がありえないが、
        /// Phred64 なら Q[0, 40] で完全に妥当<br/>
        /// この判別ができないと、
        /// 「quality - Phred - QualityCutoff が負なら捨てる」という品質フィルタが
        /// 事実上まったく効かなくなる (Q0 の塩基が Q31 に見えるため)
        /// </remarks>
        [Fact]
        public void InferOffset_実データのPhred64範囲ではPhred64と推定する()
        {
            var sample = PhredSniffer.Get_標本([new string((char)64, 4) + new string((char)104, 4)]);

            Assert.Equal(64, PhredSniffer.Get_推定オフセット(sample));
        }

        /// <summary>
        /// 典型的な Phred33 範囲では Phred33 と推定することを検証する
        /// </summary>
        [Fact]
        public void InferOffset_典型的なPhred33範囲ではPhred33と推定する()
        {
            // '!'(ASCII 33, Q0) から 'I'(ASCII 73, Q40) までの一般的な Phred33 範囲
            // Phred64 と解釈すると Q が負になるため、33 side のみが妥当
            var sample = PhredSniffer.Get_標本(["!!!!IIII"]);

            Assert.Equal(33, PhredSniffer.Get_推定オフセット(sample));
        }

        /// <summary>
        /// どちらの解釈でも現実的な範囲に収まる曖昧な範囲では null を返すことを検証する
        /// </summary>
        [Fact]
        public void InferOffset_曖昧な範囲ではnullを返す()
        {
            // ASCII 66-70 は Phred33 なら Q[33,37]、Phred64 なら Q[2,6]
            // どちらの解釈でも現実的な範囲に収まるため判別できない
            var sample = PhredSniffer.Get_標本(["BCDEF"]);

            Assert.Null(PhredSniffer.Get_推定オフセット(sample));
        }

        /// <summary>
        /// 空の標本では null を返すことを検証する
        /// </summary>
        [Fact]
        public void InferOffset_空の標本ではnullを返す()
        {
            Assert.Null(PhredSniffer.Get_推定オフセット(PhredSniffer.Get_標本([])));
        }
    }
}
