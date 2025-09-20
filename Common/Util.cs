namespace Tsumiki.Common
{
    internal class Util
    {
        private static class NucleotideID
        {
            public const byte A = 1;
            public const byte C = 2;
            public const byte G = 3;
            public const byte T = 4;
        }

        public static Span<byte> ReverseComprement(Span<byte> genome)
        {
            var buffer = new byte[genome.Length];
            for (var i = 0; i < genome.Length; i++)
            {
                buffer[genome.Length - 1 - i] = genome[i] switch
                {
                    NucleotideID.A => NucleotideID.T,
                    NucleotideID.C => NucleotideID.G,
                    NucleotideID.G => NucleotideID.C,
                    NucleotideID.T => NucleotideID.A,
                    _ => genome[i]
                };
            }
            return buffer.AsSpan();
        }

        public static Span<byte[]> ReverseComprement(Span<byte[]> genome)
        {
            var buffer = new byte[genome.Length][];
            for (var i = 0; i < genome.Length; i++)
            {
                var arr = genome[genome.Length - 1 - i];
                var newArr = new byte[arr.Length];
                for (var j = 0; j < arr.Length; j++)
                {
                    newArr[j] = arr[j] switch
                    {
                        NucleotideID.A => NucleotideID.T,
                        NucleotideID.C => NucleotideID.G,
                        NucleotideID.G => NucleotideID.C,
                        NucleotideID.T => NucleotideID.A,
                        _ => arr[j]
                    };
                }
                buffer[i] = newArr;
            }
            return buffer.AsSpan();
        }

        public static List<int> GetNucleotideIDs(char baseChar)
        {
            return baseChar switch
            {
                'A' => [NucleotideID.A],
                'M' => [NucleotideID.A, NucleotideID.C],
                'V' => [NucleotideID.A, NucleotideID.C, NucleotideID.G],
                'N' => [NucleotideID.A, NucleotideID.C, NucleotideID.G, NucleotideID.T],
                'H' => [NucleotideID.A, NucleotideID.C, NucleotideID.T],
                'R' => [NucleotideID.A, NucleotideID.G],
                'D' => [NucleotideID.A, NucleotideID.G, NucleotideID.T],
                'W' => [NucleotideID.A, NucleotideID.T],
                'C' => [NucleotideID.C],
                'S' => [NucleotideID.C, NucleotideID.G],
                'B' => [NucleotideID.C, NucleotideID.G, NucleotideID.T],
                'Y' => [NucleotideID.C, NucleotideID.T],
                'G' => [NucleotideID.G],
                'K' => [NucleotideID.G, NucleotideID.T],
                'T' => [NucleotideID.T],
                _ => throw new ArgumentException($"{baseChar} is not nucleotide base code")
            };
        }

        public static byte[] ByteToNucleotideSequence(byte read)
        {
            return [.. new[] { (read >>> 6) & 3, (read >>> 4) & 3, (read >>> 2) & 3, read & 3 }
                .Select(x => (x + 1) switch
                {
                    NucleotideID.A => NucleotideID.A,
                    NucleotideID.C => NucleotideID.C,
                    NucleotideID.G => NucleotideID.G,
                    NucleotideID.T => NucleotideID.T,
                    _ => throw new ArgumentException($"{x + 1} is not the expected value for a base")
                })];
        }

        public static string ByteToBaseString(byte read)
        {
            return read switch
            {
                NucleotideID.A => "A",
                NucleotideID.C => "C",
                NucleotideID.G => "G",
                NucleotideID.T => "T",
                _ => "N",
            };
        }

        public static List<byte[]> ToByteList(string read)
        {
            return [.. read.Select<char, byte[]>(c => c switch
            {
                'A' => [NucleotideID.A],
                'M' => [NucleotideID.A, NucleotideID.C],
                'V' => [NucleotideID.A, NucleotideID.C, NucleotideID.G],
                'N' => [NucleotideID.A, NucleotideID.C, NucleotideID.G, NucleotideID.T],
                'H' => [NucleotideID.A, NucleotideID.C, NucleotideID.T],
                'R' => [NucleotideID.A, NucleotideID.G],
                'D' => [NucleotideID.A, NucleotideID.G, NucleotideID.T],
                'W' => [NucleotideID.A, NucleotideID.T],
                'C' => [NucleotideID.C],
                'S' => [NucleotideID.C, NucleotideID.G],
                'B' => [NucleotideID.C, NucleotideID.G, NucleotideID.T],
                'Y' => [NucleotideID.C, NucleotideID.T],
                'G' => [NucleotideID.G],
                'K' => [NucleotideID.G, NucleotideID.T],
                'T' => [NucleotideID.T],
                _ => throw new ArgumentException($"{c} is not nucleotide base code")
            })];
        }

        public static ulong Pow(ulong value, long exp)
        {
            var ans = 1UL;
            while (exp > 0)
            {
                if ((exp & 1) > 0)
                {
                    ans *= value;
                }
                value *= value;
                exp >>= 1;
            }
            return ans;
        }
    }
}
