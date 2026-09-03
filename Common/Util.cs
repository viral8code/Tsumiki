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
        /// "2G" / "512M" / "1.5g" / "2048" のようなサイズ指定をバイト数に変換する。
        ///
        /// 接尾辞 K/M/G/T は2進接頭辞として扱う(1K = 1024)。メモリ量の指定なので
        /// 1000 ではなく 1024 刻みのほうが直感に合う。
        /// 接尾辞が無い場合は MB とみなす(-mem を数値だけで指定したときの単位)。
        /// </summary>
        public static long ParseMemorySize(string text)
        {
            if (string.IsNullOrWhiteSpace(text))
            {
                throw new ArgumentException("Memory size must not be empty (e.g. 2G, 512M, 1024)");
            }

            var trimmed = text.Trim();
            // "2GB" のように B が付いていても受け付ける。
            if (trimmed.Length >= 2 && (trimmed[^1] is 'B' or 'b') && !char.IsDigit(trimmed[^2]))
            {
                trimmed = trimmed[..^1];
            }

            var multiplier = 1024L * 1024L; // 接尾辞なしは MB
            var lastChar = trimmed[^1];
            switch (lastChar)
            {
                case 'K' or 'k':
                    multiplier = 1024L;
                    trimmed = trimmed[..^1];
                    break;
                case 'M' or 'm':
                    multiplier = 1024L * 1024L;
                    trimmed = trimmed[..^1];
                    break;
                case 'G' or 'g':
                    multiplier = 1024L * 1024L * 1024L;
                    trimmed = trimmed[..^1];
                    break;
                case 'T' or 't':
                    multiplier = 1024L * 1024L * 1024L * 1024L;
                    trimmed = trimmed[..^1];
                    break;
                default:
                    break;
            }

            if (!double.TryParse(trimmed, System.Globalization.NumberStyles.Float,
                    System.Globalization.CultureInfo.InvariantCulture, out var value) || value <= 0)
            {
                throw new ArgumentException($"Could not read '{text}' as a memory size (e.g. 2G, 512M, 1024)");
            }

            var bytes = (long)(value * multiplier);
            if (bytes <= 0)
            {
                throw new ArgumentException($"Memory size '{text}' is too small");
            }
            return bytes;
        }

        /// <summary>
        /// バイト数を "2 G" のような読みやすい形に戻す(パラメータ表示用)。
        /// </summary>
        public static string FormatMemorySize(long bytes)
        {
            string[] units = ["", "K", "M", "G", "T"];
            double size = bytes;
            var unitIndex = 0;
            while (size >= 1024 && unitIndex < units.Length - 1)
            {
                size /= 1024;
                unitIndex++;
            }
            return $"{size:0.#} {units[unitIndex]}B";
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