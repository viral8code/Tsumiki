using System.Text;

namespace Tsumiki.Common
{
    internal class Util
    {
        private static class NucleotideID
        {
            public const int A = 0;
            public const int C = 1;
            public const int G = 2;
            public const int T = 3;
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

        public static string CanonicalRead(string read)
        {
            var reversedRead = Util.ReverseComprement(read);
            if (read.CompareTo(reversedRead) <= 0)
            {
                return read;
            }
            else
            {
                return reversedRead;
            }
        }
    }
}
