using System.Text;
using static Tsumiki.Common.Consts;

namespace Tsumiki.Common
{
    internal class Util
    {
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

        public static string ReverseComprement(string genome)
        {
            StringBuilder sb = new();
            for (var i = 0; i < genome.Length; i++)
            {
                _ = sb.Append(genome[i] switch
                {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => throw new ArgumentException($"{genome[i]} is not the expected value for a base")
                });
            }
            return string.Join(string.Empty, sb.ToString().Reverse());
        }

        /// <summary>
        /// FASTQ の生リードなど、N を含む曖昧塩基が混入しうる文字列向けの逆相補。
        /// A/C/G/T 以外の文字はそのまま(位置だけ反転して)出力し、例外を投げない。
        /// unitig/contig 配列(Bloom filter を通過した A/C/G/T のみの配列)には
        /// このメソッドを使わないこと。そちらは ReverseComprement(string) の方を使い、
        /// 想定外の文字が混入していた場合は例外で早期検知する。
        /// </summary>
        public static string ReverseComprementAllowAmbiguous(string genome)
        {
            StringBuilder sb = new();
            for (var i = 0; i < genome.Length; i++)
            {
                _ = sb.Append(genome[i] switch
                {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    var c => c,
                });
            }
            return string.Join(string.Empty, sb.ToString().Reverse());
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

        /// <summary>
        /// 単一の塩基文字を ID に変換する軽量版。
        /// A/C/G/T はそのまま NucleotideID を返し、それ以外(曖昧塩基)は
        /// 一律 Consts.InvalidBase を返す。曖昧塩基を無視する経路
        /// (LoadReadFileToBloomFilterIgnoreAmbiguity 等)専用で、
        /// GetNucleotideIDs のような List 確保を伴わないため高速。
        /// </summary>
        public static byte GetSimpleNucleotideID(char baseChar)
        {
            return baseChar switch
            {
                'A' => NucleotideID.A,
                'C' => NucleotideID.C,
                'G' => NucleotideID.G,
                'T' => NucleotideID.T,
                _ => InvalidBase,
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

        /// <summary>
        /// 曖昧塩基を無視する経路向けの軽量版。read の各文字を1バイトIDに変換する。
        /// A/C/G/T 以外は Consts.InvalidBase になる。ToByteList と異なり
        /// LINQ・per-char の byte[] アロケーションを行わないため大幅に高速。
        /// </summary>
        public static byte[] ToSimpleByteArray(string read)
        {
            var result = new byte[read.Length];
            for (var i = 0; i < read.Length; i++)
            {
                result[i] = GetSimpleNucleotideID(read[i]);
            }
            return result;
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

        public static bool HasNext(BinaryReader stream)
        {
            return stream.BaseStream.Position < stream.BaseStream.Length;
        }

        /// <summary>
        /// FASTQ のリード ID から、ペア判定に使うための「ベース部分」を取り出す。
        /// 対応する例:
        ///   "@READ001/1"                       -> "@READ001"
        ///   "@READ001/2"                       -> "@READ001"
        ///   "@INST:RUN:FLOWCELL:1:1:1:1 1:N:0:1" -> "@INST:RUN:FLOWCELL:1:1:1:1"
        ///   "@INST:RUN:FLOWCELL:1:1:1:1 2:N:0:1" -> "@INST:RUN:FLOWCELL:1:1:1:1"
        /// 上記どちらの記法にも当てはまらない場合は ID をそのまま返す
        /// (この場合、呼び出し側で「ペアかどうか」の確証が得られないことに注意)。
        /// </summary>
        public static string GetPairedReadBaseId(string id)
        {
            // Casava 1.8+ 形式: 空白区切りの後半が "1:..." または "2:..." で始まる。
            var spaceIndex = id.IndexOf(' ');
            if (spaceIndex >= 0 && spaceIndex + 1 < id.Length)
            {
                var suffix = id[(spaceIndex + 1)..];
                if (suffix.Length > 1 && suffix[1] == ':' && (suffix[0] == '1' || suffix[0] == '2'))
                {
                    return id[..spaceIndex];
                }
            }

            // 旧来の "/1", "/2" 形式。
            if (id.Length > 2 && id[^2] == '/' && (id[^1] == '1' || id[^1] == '2'))
            {
                return id[..^2];
            }

            // "/A", "/B" のような表記に対応する亜種も一応見ておく。
            return id.Length > 2 && id[^2] == '/' && (id[^1] == 'A' || id[^1] == 'B') ? id[..^2] : id;
        }
    }
}