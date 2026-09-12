using Tsumiki.IO;

namespace Tsumiki.Tests.IO
{
    /// <summary>
    /// クオリティ文字から Phred オフセットを推定する処理の検証
    /// </summary>
    public class PhredSnifferTests
    {
        #region 公開メソッド

        /// <summary>
        /// 標本の全行を通じた最小・最大 ASCII を求めることを検証する
        /// </summary>
        [Fact]
        public void V_全行を通じた最小最大ASCIIを求める()
        {
            var l_標本 = PhredSniffer.Get_標本(["hhhh", "IIhh", "!!!!"]);

            Assert.Equal('!', (char)l_標本.A_最小ASCII);
            Assert.Equal('h', (char)l_標本.A_最大ASCII);
            Assert.Equal(3, l_標本.A_標本リード数);
            Assert.Equal(12, l_標本.A_標本文字数);
        }

        /// <summary>
        /// 標本抽出が指定した上限リード数で打ち切られることを検証する
        /// </summary>
        [Fact]
        public void V_指定した上限リード数で打ち切る()
        {
            var l_行一覧 = Enumerable.Repeat("hhhh", 100);

            var l_標本 = PhredSniffer.Get_標本(l_行一覧, p_標本上限: 5);

            Assert.Equal(5, l_標本.A_標本リード数);
        }

        /// <summary>
        /// すべての文字が同一なら一様かが真になることを検証する
        /// </summary>
        [Fact]
        public void V_すべての文字が同一なら真になる()
        {
            var l_標本 = PhredSniffer.Get_標本(["hhhh", "hhhhhh"]);

            Assert.True(l_標本.A_Is一様);
        }

        /// <summary>
        /// クオリティがばらつけば一様かが偽になることを検証する
        /// </summary>
        [Fact]
        public void V_クオリティがばらつけば偽になる()
        {
            var l_標本 = PhredSniffer.Get_標本(["hhIh", "hhhh"]);

            Assert.False(l_標本.A_Is一様);
        }

        /// <summary>
        /// 現実的な Phred33 データでは警告文が出ないことを検証する
        /// </summary>
        [Fact]
        public void V_現実的なPhred33データでは警告が出ない()
        {
            // 典型的な Phred33 品質文字 ('#'=Q2 〜 'J'=Q41 相当) を模した、
            // ばらつきのあるサンプル
            var l_標本 = PhredSniffer.Get_標本(["#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"]);

            var l_警告 = PhredSniffer.Get_警告文(l_標本, p_有効オフセット: 33);

            Assert.Null(l_警告);
        }

        /// <summary>
        /// Phred33 で ASCII が非現実的に高いときは Phred64 を試すよう促すことを検証する
        /// </summary>
        [Fact]
        public void V_Phred33でASCIIが非現実的に高いときはPhred64を促す()
        {
            // 'h' = ASCII 104
            // Phred33 なら Q=71 となり非現実的
            var l_標本 = PhredSniffer.Get_標本(["hhhh"]);

            var l_警告 = PhredSniffer.Get_警告文(l_標本, p_有効オフセット: 33);

            Assert.NotNull(l_警告);
            Assert.Contains("Phred64", l_警告);
        }

        /// <summary>
        /// Phred64 なら妥当な Q になるが、一様性はなお警告することを検証する
        /// </summary>
        [Fact]
        public void V_Phred64なら妥当なQになるが一様性はなお警告する()
        {
            // 'h' = ASCII 104
            // Phred64 なら Q=40 で妥当な範囲だが、
            // 全く同一の値しか出ていない点は別途警告する
            var l_標本 = PhredSniffer.Get_標本(["hhhh", "hhhh", "hhhh"]);

            var l_警告 = PhredSniffer.Get_警告文(l_標本, p_有効オフセット: 64);

            Assert.NotNull(l_警告);
            Assert.Contains("uniform", l_警告);
            Assert.DoesNotContain("Phred33", l_警告);
        }

        /// <summary>
        /// 負の Q になる場合は、一様性に関わらず警告することを検証する
        /// </summary>
        [Fact]
        public void V_負のQになる場合は一様性に関わらず警告する()
        {
            // '!' = ASCII 33
            // Phred64 なら Q=-31 となり明らかに不正
            var l_標本 = PhredSniffer.Get_標本(["!!!!"]);

            var l_警告 = PhredSniffer.Get_警告文(l_標本, p_有効オフセット: 64);

            Assert.NotNull(l_警告);
            Assert.Contains("Phred33", l_警告);
        }

        /// <summary>
        /// 空の標本では警告文が null になることを検証する
        /// </summary>
        [Fact]
        public void V_空の標本では警告文がnullになる()
        {
            var l_標本 = PhredSniffer.Get_標本([]);

            Assert.Null(PhredSniffer.Get_警告文(l_標本, 33));
        }

        /// <summary>
        /// 実データ (Achromobacter の IS350 ライブラリ) で観測された ASCII 範囲[64, 104]
        /// </summary>
        /// <remarks>
        /// Phred33 と解釈すると Q[31, 71] となり上限がありえないが、Phred64 なら Q[0, 40] で完全に妥当<br/>
        /// この判別ができないと、「quality - Phred - QualityCutoff が負なら捨てる」という品質フィルタが事実上まったく効かなくなる (Q0 の塩基が Q31 に見えるため)
        /// </remarks>
        [Fact]
        public void V_実データのPhred64範囲ではPhred64と推定する()
        {
            var l_標本 = PhredSniffer.Get_標本([new string((char)64, 4) + new string((char)104, 4)]);

            Assert.Equal(64, PhredSniffer.Get_推定オフセット(l_標本));
        }

        /// <summary>
        /// 典型的な Phred33 範囲では Phred33 と推定することを検証する
        /// </summary>
        [Fact]
        public void V_典型的なPhred33範囲ではPhred33と推定する()
        {
            // '!' (ASCII 33, Q0) から 'I' (ASCII 73, Q40) までの一般的な Phred33 範囲
            // Phred64 と解釈すると Q が負になるため、33 side のみが妥当
            var l_標本 = PhredSniffer.Get_標本(["!!!!IIII"]);

            Assert.Equal(33, PhredSniffer.Get_推定オフセット(l_標本));
        }

        /// <summary>
        /// どちらの解釈でも現実的な範囲に収まる曖昧な範囲では null を返すことを検証する
        /// </summary>
        [Fact]
        public void V_曖昧な範囲ではnullを返す()
        {
            // ASCII 66-70 は Phred33 なら Q[33,37]、Phred64 なら Q[2,6]
            // どちらの解釈でも現実的な範囲に収まるため判別できない
            var l_標本 = PhredSniffer.Get_標本(["BCDEF"]);

            Assert.Null(PhredSniffer.Get_推定オフセット(l_標本));
        }

        /// <summary>
        /// 空の標本では null を返すことを検証する
        /// </summary>
        [Fact]
        public void V_空の標本では推定オフセットがnullになる()
        {
            Assert.Null(PhredSniffer.Get_推定オフセット(PhredSniffer.Get_標本([])));
        }

        #endregion

    }
}
