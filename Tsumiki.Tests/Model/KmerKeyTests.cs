using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Tests.Model
{
    public class KmerKeyTests
    {
        private static void SetKmerLength(int k)
        {
            ConfigurationManager.Arguments = new Parameters { Kmer = k };
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
                    'A' => Consts.NucleotideID.A,
                    'C' => Consts.NucleotideID.C,
                    'G' => Consts.NucleotideID.G,
                    'T' => Consts.NucleotideID.T,
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

            var expected = new KmerKey(Util.ReverseComprement(forward).AsSpan());
            var actual = new KmerKey(forward.AsSpan()).ReverseComprement();

            Assert.True(expected.Equals(actual), $"expected Data=[{string.Join(",", expected.Data)}] actual Data=[{string.Join(",", actual.Data)}]");
        }

        [Theory]
        [InlineData(4)]
        [InlineData(31)]
        [InlineData(33)]
        public void Canonical_IsSameForKmerAndItsReverseComplement(int k)
        {
            SetKmerLength(k);
            var forward = "ACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAAACGTGGCCTTAA"[..k];
            var reverse = Util.ReverseComprement(forward);

            var forwardKey = new KmerKey(forward.AsSpan());
            var reverseKey = new KmerKey(reverse.AsSpan());

            Assert.True(forwardKey.Canonical().Equals(reverseKey.Canonical()));
        }

        [Fact]
        public void Canonical_IsIdempotent()
        {
            SetKmerLength(31);
            var key = new KmerKey("ACGTGGCCTTAAACGTGGCCTTAAACGTG".PadRight(31, 'A').AsSpan());

            var canonical = key.Canonical();

            Assert.True(canonical.Equals(canonical.Canonical()));
        }

        [Fact]
        public void Canonical_DifferentKmers_RemainDistinct()
        {
            SetKmerLength(4);
            var a = new KmerKey("ACGT".AsSpan()).Canonical();
            var b = new KmerKey("TTTT".AsSpan()).Canonical();

            Assert.False(a.Equals(b));
        }
    }
}
