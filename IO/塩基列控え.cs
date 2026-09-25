using System.Collections.Concurrent;
using System.Text;
using Tsumiki.Utilities;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTQ の塩基列だけを詰めた控えを作り、2 回目からは FASTQ を読まずに塩基列を返す
    /// </summary>
    internal static class 塩基列控え
    {
        #region 定数

        /// <summary>
        /// 控えの拡張子
        /// </summary>
        private const string C_控えの拡張子 = ".bases";

        /// <summary>
        /// 1 バイトに詰める塩基数
        /// </summary>
        private const int C_バイトあたりの塩基数 = 4;

        /// <summary>
        /// 読み書きのバッファの大きさ
        /// </summary>
        private const int C_バッファの大きさ = 1 << 20;

        /// <summary>
        /// 詰めたバイト列をスタックに置く大きさの上限
        /// </summary>
        private const int C_スタックに置くバイト数の上限 = 1024;

        /// <summary>
        /// 詰めずにそのまま書いた記録の目印 (長さの最下位ビット)
        /// </summary>
        private const ulong C_そのままの目印 = 1UL;

        /// <summary>
        /// 詰めた塩基の並び (2 bit の値の順)
        /// </summary>
        private const string C_塩基の並び = "ACGT";

        /// <summary>
        /// パスの大文字小文字を区別するかは OS で違う
        /// </summary>
        private static readonly StringComparer C_パスの比較 = OperatingSystem.IsWindows() ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;

        /// <summary>
        /// 詰めた 1 バイトを 4 文字に戻す表
        /// </summary>
        private static readonly char[] C_展開表 = Get_展開表();

        #endregion

        #region 内部変数

        /// <summary>
        /// 元のパスごとの、作り終えた控えの場所と、作ったときの元ファイルの版
        /// </summary>
        private static readonly ConcurrentDictionary<string, (string A_控えパス, (long A_長さ, DateTime A_更新時刻) A_版)> _控え群 = new(C_パスの比較);

        /// <summary>
        /// 控えを作っている途中の元のパス
        /// </summary>
        private static readonly ConcurrentDictionary<string, byte> _作成中 = new(C_パスの比較);

        /// <summary>
        /// 控えの名前に振る通し番号
        /// </summary>
        private static int _通し番号;

        #endregion

        #region プロパティ

        /// <summary>
        /// 控えを置くディレクトリ、null なら控えを作らない
        /// </summary>
        public static string? A_置き場所 { get; set; }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// FASTQ の塩基列を先頭から順に返す
        /// </summary>
        /// <param name="p_パス">FASTQ のパス</param>
        /// <returns>塩基列の並び</returns>
        public static IEnumerable<string> Get_塩基列(string p_パス)
        {
            if (A_置き場所 is null)
            {
                return Get_FASTQの塩基列(p_パス);
            }

            var l_キー = Path.GetFullPath(p_パス);
            if (_控え群.TryGetValue(l_キー, out var l_控え))
            {
                if (l_控え.A_版 == Get_版(p_パス))
                {
                    return Get_控えの塩基列(l_控え.A_控えパス);
                }

                V_無効化(p_パス);
            }

            return _作成中.TryAdd(l_キー, 0) ? Get_控えながらの塩基列(p_パス, l_キー) : Get_FASTQの塩基列(p_パス);
        }

        /// <summary>
        /// 元のファイルが書き換わる・消えるときに、その控えを捨てる
        /// </summary>
        /// <param name="p_パス">元のパス</param>
        public static void V_無効化(string p_パス)
        {
            if (!_控え群.IsEmpty && _控え群.TryRemove(Path.GetFullPath(p_パス), out var l_控え))
            {
                中間データ置き場.V_削除(l_控え.A_控えパス);
            }
        }

        /// <summary>
        /// すべての控えを捨てる
        /// </summary>
        public static void V_全消去()
        {
            foreach (var l_キー in _控え群.Keys)
            {
                V_無効化(l_キー);
            }
        }

        #endregion

        #region テストメソッド

        /// <summary>
        /// 作り終えた控えがあるか
        /// </summary>
        /// <param name="p_パス">元のパス</param>
        /// <returns></returns>
        public static bool Has控え(string p_パス)
        {
            return _控え群.ContainsKey(Path.GetFullPath(p_パス));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// FASTQ を読んで塩基列を返す
        /// </summary>
        /// <param name="p_パス">FASTQ のパス</param>
        /// <returns></returns>
        private static IEnumerable<string> Get_FASTQの塩基列(string p_パス)
        {
            using var l_読み込み = new FastqReader(p_パス);
            while (l_読み込み.Has続き())
            {
                yield return l_読み込み.Get_次のレコード().A_配列;
            }
        }

        /// <summary>
        /// FASTQ を読んで塩基列を返しながら控えを書き、最後まで読めたら控えを登録する
        /// </summary>
        /// <param name="p_パス">FASTQ のパス</param>
        /// <param name="p_キー">控え群のキー</param>
        /// <returns></returns>
        private static IEnumerable<string> Get_控えながらの塩基列(string p_パス, string p_キー)
        {
            var l_版 = Get_版(p_パス);
            var l_控えパス = Path.Combine(A_置き場所!, $"{Interlocked.Increment(ref _通し番号)}{C_控えの拡張子}");
            var l_Is完了 = false;
            try
            {
                if (!中間データ置き場.A_Is有効)
                {
                    _ = Directory.CreateDirectory(A_置き場所!);
                }

                using (var l_書き込み = new BufferedStream(中間データ置き場.Get_書込ストリーム(l_控えパス), C_バッファの大きさ))
                {
                    foreach (var l_配列 in Get_FASTQの塩基列(p_パス))
                    {
                        V_書込_1本(l_書き込み, l_配列);
                        yield return l_配列;
                    }
                }

                l_Is完了 = true;
            }
            finally
            {
                if (l_Is完了 && Get_版(p_パス) == l_版)
                {
                    _控え群[p_キー] = (l_控えパス, l_版);
                }
                else
                {
                    中間データ置き場.V_削除(l_控えパス);
                }

                _ = _作成中.TryRemove(p_キー, out _);
            }
        }

        /// <summary>
        /// 控えから塩基列を返す
        /// </summary>
        /// <param name="p_控えパス">控えのパス</param>
        /// <returns></returns>
        private static IEnumerable<string> Get_控えの塩基列(string p_控えパス)
        {
            using var l_読み込み = new BufferedStream(中間データ置き場.Get_読込ストリーム(p_控えパス), C_バッファの大きさ);
            var l_バッファ = new byte[256];
            while (Get_可変長整数(l_読み込み) is { } l_頭)
            {
                var l_長さ = (int)(l_頭 >> 1);
                var l_Isそのまま = (l_頭 & C_そのままの目印) != 0;
                var l_バイト数 = l_Isそのまま ? l_長さ : (l_長さ + C_バイトあたりの塩基数 - 1) / C_バイトあたりの塩基数;
                if (l_バッファ.Length < l_バイト数)
                {
                    l_バッファ = new byte[l_バイト数];
                }

                l_読み込み.ReadExactly(l_バッファ, 0, l_バイト数);
                yield return l_Isそのまま ? Encoding.UTF8.GetString(l_バッファ, 0, l_バイト数) : Get_展開済み(l_バッファ, l_長さ);
            }
        }

        /// <summary>
        /// 1 本の塩基列を控えに書く
        /// </summary>
        /// <param name="p_書き込み">書き込み先</param>
        /// <param name="p_配列">塩基列</param>
        private static void V_書込_1本(Stream p_書き込み, string p_配列)
        {
            if (!Is詰められる(p_配列))
            {
                var l_バイト列 = Encoding.UTF8.GetBytes(p_配列);
                V_書込_可変長整数(p_書き込み, ((ulong)l_バイト列.Length << 1) | C_そのままの目印);
                p_書き込み.Write(l_バイト列);
                return;
            }

            V_書込_可変長整数(p_書き込み, (ulong)p_配列.Length << 1);
            var l_バイト数 = (p_配列.Length + C_バイトあたりの塩基数 - 1) / C_バイトあたりの塩基数;
            var l_詰め = l_バイト数 <= C_スタックに置くバイト数の上限 ? stackalloc byte[l_バイト数] : new byte[l_バイト数];
            l_詰め.Clear();
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_ずらし = 2 * (C_バイトあたりの塩基数 - 1 - (i % C_バイトあたりの塩基数));
                l_詰め[i / C_バイトあたりの塩基数] |= (byte)(Get_2bit値(p_配列[i]) << l_ずらし);
            }

            p_書き込み.Write(l_詰め);
        }

        /// <summary>
        /// 大文字の ACGT だけでできているか
        /// </summary>
        /// <param name="p_配列">塩基列</param>
        /// <returns></returns>
        private static bool Is詰められる(string p_配列)
        {
            return !p_配列.AsSpan().ContainsAnyExcept(C_塩基の並び);
        }

        /// <summary>
        /// 塩基 1 文字の 2 bit の値
        /// </summary>
        /// <param name="p_塩基">A・C・G・T のどれか</param>
        /// <returns></returns>
        private static int Get_2bit値(char p_塩基)
        {
            return p_塩基 switch
            {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                _ => 3,
            };
        }

        /// <summary>
        /// 詰めたバイト列を塩基列に戻す
        /// </summary>
        /// <param name="p_詰め">詰めたバイト列</param>
        /// <param name="p_長さ">塩基数</param>
        /// <returns></returns>
        private static string Get_展開済み(byte[] p_詰め, int p_長さ)
        {
            return string.Create(p_長さ, p_詰め, static (l_文字列, l_詰め) =>
            {
                for (var i = 0; i < l_文字列.Length; i++)
                {
                    l_文字列[i] = C_展開表[(l_詰め[i / C_バイトあたりの塩基数] * C_バイトあたりの塩基数) + (i % C_バイトあたりの塩基数)];
                }
            });
        }

        /// <summary>
        /// 1 バイトの 256 通りそれぞれを 4 文字に戻した表
        /// </summary>
        /// <returns></returns>
        private static char[] Get_展開表()
        {
            var l_表 = new char[256 * C_バイトあたりの塩基数];
            for (var l_値 = 0; l_値 < 256; l_値++)
            {
                for (var j = 0; j < C_バイトあたりの塩基数; j++)
                {
                    l_表[(l_値 * C_バイトあたりの塩基数) + j] = C_塩基の並び[(l_値 >> (2 * (C_バイトあたりの塩基数 - 1 - j))) & 3];
                }
            }

            return l_表;
        }

        /// <summary>
        /// 可変長の整数を書く
        /// </summary>
        /// <param name="p_書き込み">書き込み先</param>
        /// <param name="p_値">値</param>
        private static void V_書込_可変長整数(Stream p_書き込み, ulong p_値)
        {
            while (p_値 >= 0x80)
            {
                p_書き込み.WriteByte((byte)(p_値 | 0x80));
                p_値 >>= 7;
            }

            p_書き込み.WriteByte((byte)p_値);
        }

        /// <summary>
        /// 可変長の整数を読む
        /// </summary>
        /// <param name="p_読み込み">読み込み元</param>
        /// <returns>読んだ値、末尾に達していれば null</returns>
        private static ulong? Get_可変長整数(Stream p_読み込み)
        {
            var l_値 = 0UL;
            for (var l_ずらし = 0; ; l_ずらし += 7)
            {
                var l_バイト = p_読み込み.ReadByte();
                if (l_バイト < 0)
                {
                    return l_ずらし == 0 ? null : throw new InvalidDataException("塩基列の控えが途中で切れている");
                }

                l_値 |= (ulong)(l_バイト & 0x7F) << l_ずらし;
                if (l_バイト < 0x80)
                {
                    return l_値;
                }
            }
        }

        /// <summary>
        /// 元ファイルの版 (ディスクに無ければ固定値)
        /// </summary>
        /// <param name="p_パス">元のパス</param>
        /// <returns></returns>
        private static (long A_長さ, DateTime A_更新時刻) Get_版(string p_パス)
        {
            var l_情報 = new FileInfo(p_パス);
            return l_情報.Exists ? (l_情報.Length, l_情報.LastWriteTimeUtc) : (-1L, default);
        }

        #endregion
    }
}
