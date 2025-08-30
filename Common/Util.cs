using System.Text;

namespace Tsumiki.Common
{
    internal class Util
    {
        private static class NucleotideID
        {
            public const int A = 1;
            public const int C = 2;
            public const int G = 3;
            public const int T = 4;
        }

        public static string ReverseComprement(string genome)
        {
            var sb = new StringBuilder();
            foreach (var baseChar in genome.Reverse())
            {
                _ = sb.Append(baseChar switch
                {
                    'A' => 'T',
                    'B' => 'A',
                    'C' => 'G',
                    'D' => 'C',
                    'G' => 'C',
                    'H' => 'G',
                    'K' => 'M',
                    'M' => 'K',
                    'N' => 'N',
                    'R' => 'Y',
                    'S' => 'W',
                    'T' => 'A',
                    'V' => 'T',
                    'W' => 'S',
                    'Y' => 'R',
                    _ => throw new ArgumentException($"{baseChar} is not nucleotide base code")
                });
            }
            return sb.ToString();
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

        public static string ByteToNucleotideSequence(byte read)
        {
            return string.Join(string.Empty, new[] { (read >>> 6) & 3, (read >>> 4) & 3, (read >>> 2) & 3, read & 3 }
                .Select(x => (x + 1) switch
                {
                    NucleotideID.A => "A",
                    NucleotideID.C => "C",
                    NucleotideID.G => "G",
                    NucleotideID.T => "T",
                    _ => "N",
                }));
        }

        public static string ByteToBase(byte read)
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
    }
}
