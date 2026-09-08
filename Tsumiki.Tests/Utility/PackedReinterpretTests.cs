using Tsumiki.Common;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// ディスク上のパック済みバイト列を、塩基列へ展開せずそのままパック値として
    /// 読み替えられることを固定する。並びの規約が同じなので余りビット分の
    /// シフトだけで一致するはずで、ここがずれると信頼 k-mer 集合が丸ごと変わる。
    /// </summary>
    public class PackedReinterpretTests
    {
        private static byte[] Get_パック済み(byte[] p_塩基ID列)
        {
            // CountingDB.V_登録 と同じ手順。端数は A で埋める。
            var l_パック = new byte[(p_塩基ID列.Length + 3) / 4];
            var l_位置 = 0;
            for (var i = 0; i < p_塩基ID列.Length; i += 4)
            {
                var l_バイト = 0;
                for (var j = 0; j < 4; j++)
                {
                    var l_塩基ID = i + j < p_塩基ID列.Length ? p_塩基ID列[i + j] : Consts.塩基ID.A;
                    l_バイト <<= 2;
                    l_バイト |= l_塩基ID - 1;
                }
                l_パック[l_位置++] = (byte)l_バイト;
            }
            return l_パック;
        }

        private static byte[] Get_塩基ID列(int p_長さ, int p_種)
        {
            var l_乱数 = new Random(p_種);
            return [.. Enumerable.Range(0, p_長さ).Select(_ => (byte)(l_乱数.Next(4) + 1))];
        }

        [Theory]
        [InlineData(21)]
        [InlineData(31)]
        [InlineData(32)]
        public void Get_読み替え_小_MatchesPackingTheExpandedBases(int p_k長)
        {
            var l_塩基 = Get_塩基ID列(p_k長, 7 + p_k長);
            var l_パック = Get_パック済み(l_塩基);
            var l_余り = (8 * l_パック.Length) - (2 * p_k長);

            Assert.Equal(
                TrustedKmerIndex.Get_パック_小(l_塩基),
                TrustedKmerIndex.Get_読み替え_小(l_パック, l_余り));
        }

        [Theory]
        [InlineData(33)]
        [InlineData(41)]
        [InlineData(63)]
        [InlineData(64)]
        public void Get_読み替え_中_MatchesPackingTheExpandedBases(int p_k長)
        {
            var l_塩基 = Get_塩基ID列(p_k長, 11 + p_k長);
            var l_パック = Get_パック済み(l_塩基);
            var l_余り = (8 * l_パック.Length) - (2 * p_k長);

            Assert.Equal(
                TrustedKmerIndex.Get_パック_中(l_塩基),
                TrustedKmerIndex.Get_読み替え_中(l_パック, l_余り));
        }
    }
}
