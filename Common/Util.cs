using System.Text;

namespace Tsumiki.Common
{
    internal class Util
    {
        public static Span<byte> V_逆相補(Span<byte> p_塩基列)
        {
            var l_結果 = new byte[p_塩基列.Length];
            for (var i = 0; i < p_塩基列.Length; i++)
            {
                l_結果[p_塩基列.Length - 1 - i] = p_塩基列[i] switch
                {
                    Consts.塩基ID.A => Consts.塩基ID.T,
                    Consts.塩基ID.C => Consts.塩基ID.G,
                    Consts.塩基ID.G => Consts.塩基ID.C,
                    Consts.塩基ID.T => Consts.塩基ID.A,
                    _ => p_塩基列[i]
                };
            }
            return l_結果.AsSpan();
        }

        public static string V_逆相補(string p_配列)
        {
            StringBuilder l_結果 = new();
            for (var i = 0; i < p_配列.Length; i++)
            {
                _ = l_結果.Append(p_配列[i] switch
                {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => throw new ArgumentException($"{p_配列[i]} is not the expected value for a base")
                });
            }
            return string.Join(string.Empty, l_結果.ToString().Reverse());
        }

        /// <summary>
        /// 環状配列の開始位置を、辞書式順序で最小になる回転へ正規化する
        /// (Booth のアルゴリズム、O(n))。
        ///
        /// 環状に閉じた contig は開始位置が任意(walk がどこから始まったかの
        /// 産物でしかない)。決定的な基準を置かないと、同じ環状配列でも
        /// 実行のたびに(あるいは同じ実行内でも walk の起点が変われば)
        /// 別の文字列として出力され、下流の比較や再現性を損なう。
        /// </summary>
        public static string Get_最小回転(string p_配列)
        {
            if (p_配列.Length <= 1)
            {
                return p_配列;
            }
            var l_開始位置 = Get_最小回転の開始位置(p_配列);
            return l_開始位置 == 0 ? p_配列 : p_配列[l_開始位置..] + p_配列[..l_開始位置];
        }

        /// <summary>
        /// Booth のアルゴリズム。p_配列 を2つ繋げた仮想文字列の上で
        /// KMP の失敗関数に似た配列を作りながら、最小回転の開始位置を求める。
        /// </summary>
        private static int Get_最小回転の開始位置(string p_配列)
        {
            var l_長さ = p_配列.Length;
            var l_二重化 = p_配列 + p_配列;
            var l_失敗関数 = new int[l_二重化.Length];
            Array.Fill(l_失敗関数, -1);
            var k = 0;
            for (var j = 1; j < l_二重化.Length; j++)
            {
                var l_文字 = l_二重化[j];
                var i = l_失敗関数[j - k - 1];
                while (i != -1 && l_文字 != l_二重化[k + i + 1])
                {
                    if (l_文字 < l_二重化[k + i + 1])
                    {
                        k = j - i - 1;
                    }
                    i = l_失敗関数[i];
                }
                if (l_文字 != l_二重化[k + i + 1])
                {
                    if (l_文字 < l_二重化[k])
                    {
                        k = j;
                    }
                    l_失敗関数[j - k] = -1;
                }
                else
                {
                    l_失敗関数[j - k] = i + 1;
                }
            }
            return k % l_長さ;
        }

        /// <summary>
        /// 曖昧塩基が混入しうる文字列向けの逆相補。A/C/G/T 以外は位置だけ反転して通す。
        /// unitig/contig には使わないこと。そちらは V_逆相補(string) を使い、
        /// 想定外の文字を例外で早期検知する。
        /// </summary>
        public static string V_逆相補_曖昧塩基あり(string p_配列)
        {
            StringBuilder l_結果 = new();
            for (var i = 0; i < p_配列.Length; i++)
            {
                _ = l_結果.Append(p_配列[i] switch
                {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    var l_文字 => l_文字,
                });
            }
            return string.Join(string.Empty, l_結果.ToString().Reverse());
        }

        public static Span<byte[]> V_逆相補(Span<byte[]> p_塩基候補列)
        {
            var l_結果 = new byte[p_塩基候補列.Length][];
            for (var i = 0; i < p_塩基候補列.Length; i++)
            {
                var l_候補 = p_塩基候補列[p_塩基候補列.Length - 1 - i];
                var l_変換後 = new byte[l_候補.Length];
                for (var j = 0; j < l_候補.Length; j++)
                {
                    l_変換後[j] = l_候補[j] switch
                    {
                        Consts.塩基ID.A => Consts.塩基ID.T,
                        Consts.塩基ID.C => Consts.塩基ID.G,
                        Consts.塩基ID.G => Consts.塩基ID.C,
                        Consts.塩基ID.T => Consts.塩基ID.A,
                        _ => l_候補[j]
                    };
                }
                l_結果[i] = l_変換後;
            }
            return l_結果.AsSpan();
        }

        /// <summary>
        /// 塩基文字が曖昧(A/C/G/T のいずれでもない IUPAC コード)かどうか。
        /// 候補の中身ではなく個数だけが必要な場面で、List の確保を避ける。
        /// </summary>
        public static bool Get_曖昧塩基か(char p_塩基文字)
        {
            return p_塩基文字 switch
            {
                'A' or 'C' or 'G' or 'T' => false,
                'M' or 'V' or 'N' or 'H' or 'R' or 'D' or 'W' or 'S' or 'B' or 'Y' or 'K' => true,
                _ => throw new ArgumentException($"{p_塩基文字} is not nucleotide base code"),
            };
        }

        public static List<int> Get_塩基ID候補(char p_塩基文字)
        {
            return p_塩基文字 switch
            {
                'A' => [Consts.塩基ID.A],
                'M' => [Consts.塩基ID.A, Consts.塩基ID.C],
                'V' => [Consts.塩基ID.A, Consts.塩基ID.C, Consts.塩基ID.G],
                'N' => [Consts.塩基ID.A, Consts.塩基ID.C, Consts.塩基ID.G, Consts.塩基ID.T],
                'H' => [Consts.塩基ID.A, Consts.塩基ID.C, Consts.塩基ID.T],
                'R' => [Consts.塩基ID.A, Consts.塩基ID.G],
                'D' => [Consts.塩基ID.A, Consts.塩基ID.G, Consts.塩基ID.T],
                'W' => [Consts.塩基ID.A, Consts.塩基ID.T],
                'C' => [Consts.塩基ID.C],
                'S' => [Consts.塩基ID.C, Consts.塩基ID.G],
                'B' => [Consts.塩基ID.C, Consts.塩基ID.G, Consts.塩基ID.T],
                'Y' => [Consts.塩基ID.C, Consts.塩基ID.T],
                'G' => [Consts.塩基ID.G],
                'K' => [Consts.塩基ID.G, Consts.塩基ID.T],
                'T' => [Consts.塩基ID.T],
                _ => throw new ArgumentException($"{p_塩基文字} is not nucleotide base code")
            };
        }

        /// <summary>
        /// 単一の塩基文字を ID に変換する軽量版。曖昧塩基は一律 Consts.無効な塩基。
        /// List 確保を伴わないため、曖昧塩基を無視する経路ではこちらを使う。
        /// </summary>
        public static byte Get_塩基ID(char p_塩基文字)
        {
            return p_塩基文字 switch
            {
                'A' => Consts.塩基ID.A,
                'C' => Consts.塩基ID.C,
                'G' => Consts.塩基ID.G,
                'T' => Consts.塩基ID.T,
                _ => Consts.無効な塩基,
            };
        }

        public static byte[] V_変換_塩基列(byte p_パック済みバイト)
        {
            return [.. new[] { (p_パック済みバイト >>> 6) & 3, (p_パック済みバイト >>> 4) & 3, (p_パック済みバイト >>> 2) & 3, p_パック済みバイト & 3 }
                .Select(x => (x + 1) switch
                {
                    Consts.塩基ID.A => Consts.塩基ID.A,
                    Consts.塩基ID.C => Consts.塩基ID.C,
                    Consts.塩基ID.G => Consts.塩基ID.G,
                    Consts.塩基ID.T => Consts.塩基ID.T,
                    _ => throw new ArgumentException($"{x + 1} is not the expected value for a base")
                })];
        }

        public static string V_変換_塩基文字(byte p_塩基ID)
        {
            return p_塩基ID switch
            {
                Consts.塩基ID.A => "A",
                Consts.塩基ID.C => "C",
                Consts.塩基ID.G => "G",
                Consts.塩基ID.T => "T",
                _ => "N",
            };
        }

        /// <summary>
        /// 塩基IDを1文字へ変換する。文字列を返す版は連結のたびに確保が起きるため、
        /// 塩基列をまとめて文字列にする場面ではこちらを使う。
        /// </summary>
        public static char Get_塩基文字(byte p_塩基ID)
        {
            return p_塩基ID switch
            {
                Consts.塩基ID.A => 'A',
                Consts.塩基ID.C => 'C',
                Consts.塩基ID.G => 'G',
                Consts.塩基ID.T => 'T',
                _ => 'N',
            };
        }

        public static List<byte[]> V_変換_塩基候補列(string p_リード)
        {
            return [.. p_リード.Select<char, byte[]>(x => x switch
            {
                'A' => [Consts.塩基ID.A],
                'M' => [Consts.塩基ID.A, Consts.塩基ID.C],
                'V' => [Consts.塩基ID.A, Consts.塩基ID.C, Consts.塩基ID.G],
                'N' => [Consts.塩基ID.A, Consts.塩基ID.C, Consts.塩基ID.G, Consts.塩基ID.T],
                'H' => [Consts.塩基ID.A, Consts.塩基ID.C, Consts.塩基ID.T],
                'R' => [Consts.塩基ID.A, Consts.塩基ID.G],
                'D' => [Consts.塩基ID.A, Consts.塩基ID.G, Consts.塩基ID.T],
                'W' => [Consts.塩基ID.A, Consts.塩基ID.T],
                'C' => [Consts.塩基ID.C],
                'S' => [Consts.塩基ID.C, Consts.塩基ID.G],
                'B' => [Consts.塩基ID.C, Consts.塩基ID.G, Consts.塩基ID.T],
                'Y' => [Consts.塩基ID.C, Consts.塩基ID.T],
                'G' => [Consts.塩基ID.G],
                'K' => [Consts.塩基ID.G, Consts.塩基ID.T],
                'T' => [Consts.塩基ID.T],
                _ => throw new ArgumentException($"{x} is not nucleotide base code")
            })];
        }

        /// <summary>
        /// 曖昧塩基を無視する経路向けの軽量版。リードの各文字を1バイトIDに変換する。
        /// A/C/G/T 以外は Consts.無効な塩基 になる。V_変換_塩基候補列 と異なり
        /// LINQ・per-char の byte[] アロケーションを行わないため大幅に高速。
        /// </summary>
        public static byte[] V_変換_塩基列(string p_リード)
        {
            var l_結果 = new byte[p_リード.Length];
            for (var i = 0; i < p_リード.Length; i++)
            {
                l_結果[i] = Get_塩基ID(p_リード[i]);
            }
            return l_結果;
        }

        public static ulong V_累乗(ulong p_底, long p_指数)
        {
            var l_結果 = 1UL;
            while (p_指数 > 0)
            {
                if ((p_指数 & 1) > 0)
                {
                    l_結果 *= p_底;
                }
                p_底 *= p_底;
                p_指数 >>= 1;
            }
            return l_結果;
        }

        public static bool Get_続きがあるか(BinaryReader p_ストリーム)
        {
            return p_ストリーム.BaseStream.Position < p_ストリーム.BaseStream.Length;
        }

        /// <summary>
        /// "2G" / "512M" / "2048" のようなサイズ指定をバイト数に変換する。
        /// 接尾辞は2進接頭辞(1K = 1024)、接尾辞が無い場合は MB とみなす。
        /// </summary>
        public static long V_変換_メモリサイズ(string p_表記)
        {
            if (string.IsNullOrWhiteSpace(p_表記))
            {
                throw new ArgumentException("Memory size must not be empty (e.g. 2G, 512M, 1024)");
            }

            var l_本体 = p_表記.Trim();
            // "2GB" のように B が付いていても受け付ける。
            if (l_本体.Length >= 2 && (l_本体[^1] is 'B' or 'b') && !char.IsDigit(l_本体[^2]))
            {
                l_本体 = l_本体[..^1];
            }

            var l_倍率 = 1024L * 1024L; // 接尾辞なしは MB
            switch (l_本体[^1])
            {
                case 'K' or 'k':
                    l_倍率 = 1024L;
                    l_本体 = l_本体[..^1];
                    break;
                case 'M' or 'm':
                    l_倍率 = 1024L * 1024L;
                    l_本体 = l_本体[..^1];
                    break;
                case 'G' or 'g':
                    l_倍率 = 1024L * 1024L * 1024L;
                    l_本体 = l_本体[..^1];
                    break;
                case 'T' or 't':
                    l_倍率 = 1024L * 1024L * 1024L * 1024L;
                    l_本体 = l_本体[..^1];
                    break;
                default:
                    break;
            }

            if (!double.TryParse(l_本体, System.Globalization.NumberStyles.Float,
                    System.Globalization.CultureInfo.InvariantCulture, out var l_数値) || l_数値 <= 0)
            {
                throw new ArgumentException($"Could not read '{p_表記}' as a memory size (e.g. 2G, 512M, 1024)");
            }

            var l_バイト数 = (long)(l_数値 * l_倍率);
            return l_バイト数 <= 0 ? throw new ArgumentException($"Memory size '{p_表記}' is too small") : l_バイト数;
        }

        /// <summary>
        /// バイト数を "2 GB" のような読みやすい形に戻す(パラメータ表示用)。
        /// </summary>
        public static string Get_表示用メモリサイズ(long p_バイト数)
        {
            string[] l_単位 = ["", "K", "M", "G", "T"];
            double l_サイズ = p_バイト数;
            var l_単位位置 = 0;
            while (l_サイズ >= 1024 && l_単位位置 < l_単位.Length - 1)
            {
                l_サイズ /= 1024;
                l_単位位置++;
            }
            return $"{l_サイズ:0.#} {l_単位[l_単位位置]}B";
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
        public static string Get_ペア共通ID(string p_ID)
        {
            // Casava 1.8+ 形式: 空白区切りの後半が "1:..." または "2:..." で始まる。
            var l_空白位置 = p_ID.IndexOf(' ');
            if (l_空白位置 >= 0 && l_空白位置 + 1 < p_ID.Length)
            {
                var l_後半 = p_ID[(l_空白位置 + 1)..];
                if (l_後半.Length > 1 && l_後半[1] == ':' && (l_後半[0] == '1' || l_後半[0] == '2'))
                {
                    return p_ID[..l_空白位置];
                }
            }

            // 旧来の "/1", "/2" 形式。
            if (p_ID.Length > 2 && p_ID[^2] == '/' && (p_ID[^1] == '1' || p_ID[^1] == '2'))
            {
                return p_ID[..^2];
            }

            // "/A", "/B" のような表記に対応する亜種も一応見ておく。
            return p_ID.Length > 2 && p_ID[^2] == '/' && (p_ID[^1] == 'A' || p_ID[^1] == 'B') ? p_ID[..^2] : p_ID;
        }
    }
}
