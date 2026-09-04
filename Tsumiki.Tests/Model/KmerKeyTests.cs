using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Tests.Model
{
    public class KmerKeyTests
    {
        private static void SetKmerLength(int k)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = k };
        }

        [Theory]
        [InlineData(4)]   // 1つの ulong に収まる短いk-mer
        [InlineData(31)]  // デフォルトのk-mer長
        [InlineData(33)]  // 32境界をまたぐ長さ(Dataが複数ulongになる)
        [InlineData(64)]  // ちょうど2 ulong 分
        public void ByteConstructor_And_CharConstructor_ProduceEqualKeys(int k)
        {
            SetKmerLength(k);
            var bases = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"[..k];
            var byteKmer = new byte[k];
            for (var i = 0; i < k; i++)
            {
                byteKmer[i] = bases[i] switch
                {
                    'A' => Consts.塩基ID.A,
                    'C' => Consts.塩基ID.C,
                    'G' => Consts.塩基ID.G,
                    'T' => Consts.塩基ID.T,
                    _ => throw new InvalidOperationException(),
                };
            }

            var fromChar = new KmerKey(bases.AsSpan());
            var fromByte = new KmerKey(byteKmer);

            Assert.True(fromChar.Equals(fromByte));
            Assert.Equal(fromChar.GetHashCode(), fromByte.GetHashCode());
        }

        [Theory]
        [InlineData(4)]
        [InlineData(31)]
        [InlineData(33)]
        [InlineData(64)]
        public void ReverseComprement_MatchesStringBasedReverseComplement(int k)
        {
            SetKmerLength(k);
            var forward = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"[..k];

            var expected = new KmerKey(Util.V_逆相補(forward).AsSpan());
            var actual = new KmerKey(forward.AsSpan()).Get_逆相補();

            Assert.True(expected.Equals(actual), $"expected Data=[{string.Join(",", expected.A_パック済みデータ)}] actual Data=[{string.Join(",", actual.A_パック済みデータ)}]");
        }

        [Theory]
        [InlineData(4)]
        [InlineData(31)]
        [InlineData(33)]
        public void Canonical_IsSameForKmerAndItsReverseComplement(int k)
        {
            SetKmerLength(k);
            var forward = "ACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAA"[..k];
            var reverse = Util.V_逆相補(forward);

            var forwardKey = new KmerKey(forward.AsSpan());
            var reverseKey = new KmerKey(reverse.AsSpan());

            Assert.True(forwardKey.Get_正規形().Equals(reverseKey.Get_正規形()));
        }

        [Fact]
        public void Canonical_IsIdempotent()
        {
            SetKmerLength(31);
            var key = new KmerKey("ACGTGGCCTTAAACGTGGCCTTAAACGTG".PadRight(31, 'A').AsSpan());

            var canonical = key.Get_正規形();

            Assert.True(canonical.Equals(canonical.Get_正規形()));
        }

        [Fact]
        public void Canonical_DifferentKmers_RemainDistinct()
        {
            SetKmerLength(4);
            var a = new KmerKey("ACGT".AsSpan()).Get_正規形();
            var b = new KmerKey("TTTT".AsSpan()).Get_正規形();

            Assert.False(a.Equals(b));
        }
    }
}
