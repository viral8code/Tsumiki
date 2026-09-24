using System.Collections.Concurrent;
using System.IO.Compression;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 中間データ (前処理・訂正済みリード、断片、k-mer 計数の一時ファイル) を読み書きする口
    /// </summary>
    /// <remarks>
    /// オンメモリモードではパスをキーにしてメモリ上に置き、ディスクへは書かない<br/>
    /// パスのまま受け渡す作りを変えずに済むよう、置き場に無いパスはそのままファイルとして扱う
    /// </remarks>
    internal static class 中間データ置き場
    {
        #region 定数

        /// <summary>
        /// 1 つの塊の大きさ
        /// </summary>
        /// <remarks>
        /// 1 本の配列に詰めると 2 GB を超えるデータを持てないので、塊に分けて繋ぐ
        /// </remarks>
        private const int 塊の大きさ = 16 << 20;

        /// <summary>
        /// ファイルとして開くときのバッファの大きさ
        /// </summary>
        private const int ファイルのバッファ = 1 << 20;

        /// <summary>
        /// パスの大文字小文字を区別するかは OS で違う
        /// </summary>
        private static readonly StringComparer パスの比較 = OperatingSystem.IsWindows() ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;

        /// <summary>
        /// 名前の比較 (パスの比較と揃える)
        /// </summary>
        private static readonly StringComparison 名前の比較 = OperatingSystem.IsWindows() ? StringComparison.OrdinalIgnoreCase : StringComparison.Ordinal;

        #endregion

        #region 内部変数

        /// <summary>
        /// パスごとの圧縮済みの中身
        /// </summary>
        private static readonly ConcurrentDictionary<string, (List<byte[]> A_塊群, long A_元の長さ)> _置き場 = new(パスの比較);

        /// <summary>
        /// 取り込んだ入力のキー
        /// </summary>
        /// <remarks>
        /// 利用者のファイルなので、置き場から消した後もディスク上のものは消さない
        /// </remarks>
        private static readonly ConcurrentDictionary<string, byte> _取り込み済み = new(パスの比較);

        #endregion

        #region プロパティ

        /// <summary>
        /// 中間データをメモリに置くか
        /// </summary>
        public static bool A_Is有効 { get; set; }

        /// <summary>
        /// いま置き場が抱えている圧縮済みの大きさ (バイト)
        /// </summary>
        public static long A_使用量 => _置き場.Values.Sum(x => x.A_塊群.Sum(y => (long)y.Length));

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 中間データを書き込むストリームを開く
        /// </summary>
        /// <param name="p_パス">書き込み先</param>
        /// <returns>閉じた時点で書き込みが確定するストリーム</returns>
        public static Stream Get_書込ストリーム(string p_パス)
        {
            if (!A_Is有効)
            {
                return new FileStream(p_パス, FileMode.Create, FileAccess.Write, FileShare.Read, ファイルのバッファ, FileOptions.SequentialScan);
            }

            var l_キー = Get_キー(p_パス);
            _ = _置き場.TryRemove(l_キー, out _);
            var l_塊 = new 塊書込ストリーム();
            return new 計数ストリーム(new DeflateStream(l_塊, CompressionLevel.Fastest), l_元の長さ => _置き場[l_キー] = (l_塊.A_塊群, l_元の長さ));
        }

        /// <summary>
        /// ディスク上のファイルを置き場へ写す
        /// </summary>
        /// <param name="p_パス">写すファイル</param>
        /// <remarks>
        /// 入力リードは工程ごとに何度も読み直されるので、メモリに置くときは最初に 1 回だけ読む<br/>
        /// 写した後もディスク上のファイルには触らない (削除しても置き場から消えるだけ)
        /// </remarks>
        public static void V_取り込み(string p_パス)
        {
            if (!A_Is有効 || _置き場.ContainsKey(Get_キー(p_パス)))
            {
                return;
            }

            _取り込み済み[Get_キー(p_パス)] = 0;
            using var l_元 = new FileStream(p_パス, FileMode.Open, FileAccess.Read, FileShare.Read, ファイルのバッファ, FileOptions.SequentialScan);
            using var l_先 = Get_書込ストリーム(p_パス);
            l_元.CopyTo(l_先);
        }

        /// <summary>
        /// 読み込むストリームを開く
        /// </summary>
        /// <param name="p_パス">読み込み元</param>
        /// <returns>置き場にあればその中身、無ければファイル</returns>
        public static Stream Get_読込ストリーム(string p_パス)
        {
            return _置き場.TryGetValue(Get_キー(p_パス), out var l_項目)
                ? new 計数ストリーム(new DeflateStream(new 塊読込ストリーム(l_項目.A_塊群), CompressionMode.Decompress), l_項目.A_元の長さ)
                : new FileStream(p_パス, FileMode.Open, FileAccess.Read, FileShare.Read, ファイルのバッファ, FileOptions.SequentialScan);
        }

        /// <summary>
        /// 置き場かディスクのどちらかにあるか
        /// </summary>
        /// <param name="p_パス">調べるパス</param>
        /// <returns></returns>
        public static bool Is存在(string? p_パス)
        {
            return !string.IsNullOrWhiteSpace(p_パス) && (_置き場.ContainsKey(Get_キー(p_パス)) || File.Exists(p_パス));
        }

        /// <summary>
        /// 置き場かディスクから消す
        /// </summary>
        /// <param name="p_パス">消すパス</param>
        public static void V_削除(string p_パス)
        {
            var l_キー = Get_キー(p_パス);
            if (!_置き場.TryRemove(l_キー, out _) && !_取り込み済み.ContainsKey(l_キー) && File.Exists(p_パス))
            {
                File.Delete(p_パス);
            }
        }

        /// <summary>
        /// ディレクトリ直下で、名前が接頭辞で始まり接尾辞で終わるものを置き場とディスクの両方から列挙する
        /// </summary>
        /// <param name="p_ディレクトリ">探す場所</param>
        /// <param name="p_接頭辞">名前の先頭</param>
        /// <param name="p_接尾辞">名前の末尾</param>
        /// <returns></returns>
        public static IReadOnlyList<string> Get_一覧(string p_ディレクトリ, string p_接頭辞, string p_接尾辞)
        {
            var l_場所 = Get_キー(p_ディレクトリ);
            var l_置き場 = _置き場.Keys.Where(x => string.Equals(Path.GetDirectoryName(x), l_場所, 名前の比較) && Is名前が一致(Path.GetFileName(x), p_接頭辞, p_接尾辞));
            IEnumerable<string> l_ディスク = Directory.Exists(p_ディレクトリ)
                ? Directory.EnumerateFiles(p_ディレクトリ).Where(x => Is名前が一致(Path.GetFileName(x), p_接頭辞, p_接尾辞))
                : [];
            return [.. l_置き場, .. l_ディスク];
        }

        /// <summary>
        /// 置き場を空にする
        /// </summary>
        public static void V_消去()
        {
            _置き場.Clear();
            _取り込み済み.Clear();
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 同じ場所を指す別表記を 1 つにまとめたキー
        /// </summary>
        /// <param name="p_パス"></param>
        /// <returns></returns>
        private static string Get_キー(string p_パス)
        {
            return Path.GetFullPath(p_パス);
        }

        /// <summary>
        /// 名前が接頭辞で始まり接尾辞で終わるか
        /// </summary>
        /// <param name="p_名前"></param>
        /// <param name="p_接頭辞"></param>
        /// <param name="p_接尾辞"></param>
        /// <returns></returns>
        private static bool Is名前が一致(string p_名前, string p_接頭辞, string p_接尾辞)
        {
            return p_名前.StartsWith(p_接頭辞, 名前の比較) && p_名前.EndsWith(p_接尾辞, 名前の比較);
        }

        #endregion

        #region 内部クラス

        /// <summary>
        /// 展開後の読み書きの量を数え、長さと位置を答える
        /// </summary>
        /// <remarks>
        /// 圧縮ストリームは長さも位置も持たないが、呼び出し側には位置と長さで終端を判定するものがある
        /// </remarks>
        private sealed class 計数ストリーム : Stream
        {
            /// <summary>
            /// 包んでいるストリーム
            /// </summary>
            private readonly Stream _中身;

            /// <summary>
            /// 書き終えたときに展開後の長さを渡す先 (読み込み用なら null)
            /// </summary>
            private readonly Action<long>? _確定;

            /// <summary>
            /// 展開後の長さ (読み込み用のときだけ分かる)
            /// </summary>
            private readonly long _長さ;

            /// <summary>
            /// ここまでに読み書きした展開後の量
            /// </summary>
            private long _位置;

            /// <summary>
            /// 閉じたか
            /// </summary>
            private bool _Is確定済み;

            /// <summary>
            /// 書き込み用
            /// </summary>
            /// <param name="p_中身">包む圧縮ストリーム</param>
            /// <param name="p_確定">閉じたときに展開後の長さを渡す先</param>
            public 計数ストリーム(Stream p_中身, Action<long> p_確定)
            {
                this._中身 = p_中身;
                this._確定 = p_確定;
            }

            /// <summary>
            /// 読み込み用
            /// </summary>
            /// <param name="p_中身">包む展開ストリーム</param>
            /// <param name="p_長さ">展開後の長さ</param>
            public 計数ストリーム(Stream p_中身, long p_長さ)
            {
                this._中身 = p_中身;
                this._長さ = p_長さ;
            }

            public override bool CanRead => this._確定 is null;
            public override bool CanSeek => false;
            public override bool CanWrite => this._確定 is not null;
            public override long Length => this._確定 is null ? this._長さ : this._位置;
            public override long Position { get => this._位置; set => throw new NotSupportedException(); }

            public override int Read(byte[] p_バッファ, int p_開始, int p_長さ)
            {
                return this.Read(p_バッファ.AsSpan(p_開始, p_長さ));
            }

            public override int Read(Span<byte> p_バッファ)
            {
                var l_読んだ = this._中身.Read(p_バッファ);
                this._位置 += l_読んだ;
                return l_読んだ;
            }

            public override void Write(byte[] p_バッファ, int p_開始, int p_長さ)
            {
                this.Write(p_バッファ.AsSpan(p_開始, p_長さ));
            }

            public override void Write(ReadOnlySpan<byte> p_バッファ)
            {
                this._中身.Write(p_バッファ);
                this._位置 += p_バッファ.Length;
            }

            public override void Flush()
            {
                this._中身.Flush();
            }

            public override long Seek(long p_位置, SeekOrigin p_起点) => throw new NotSupportedException();
            public override void SetLength(long p_長さ) => throw new NotSupportedException();

            protected override void Dispose(bool p_Is明示)
            {
                if (!this._Is確定済み)
                {
                    this._Is確定済み = true;
                    this._中身.Dispose();
                    this._確定?.Invoke(this._位置);
                }
                base.Dispose(p_Is明示);
            }
        }

        /// <summary>
        /// 塊の列へ追記する
        /// </summary>
        private sealed class 塊書込ストリーム : Stream
        {
            /// <summary>
            /// 書き終えた塊 (閉じた後は書きかけの塊も含む)
            /// </summary>
            public List<byte[]> A_塊群 { get; } = [];

            /// <summary>
            /// 書きかけの塊
            /// </summary>
            private byte[] _今の塊 = new byte[塊の大きさ];

            /// <summary>
            /// 書きかけの塊の使用量
            /// </summary>
            private int _今の位置;

            /// <summary>
            /// 閉じたか
            /// </summary>
            private bool _Is確定済み;

            public override bool CanRead => false;
            public override bool CanSeek => false;
            public override bool CanWrite => true;
            public override long Length => throw new NotSupportedException();
            public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }

            public override void Write(byte[] p_バッファ, int p_開始, int p_長さ)
            {
                this.Write(p_バッファ.AsSpan(p_開始, p_長さ));
            }

            public override void Write(ReadOnlySpan<byte> p_バッファ)
            {
                while (!p_バッファ.IsEmpty)
                {
                    if (this._今の位置 == this._今の塊.Length)
                    {
                        this.A_塊群.Add(this._今の塊);
                        this._今の塊 = new byte[塊の大きさ];
                        this._今の位置 = 0;
                    }
                    var l_書く分 = Math.Min(p_バッファ.Length, this._今の塊.Length - this._今の位置);
                    p_バッファ[..l_書く分].CopyTo(this._今の塊.AsSpan(this._今の位置));
                    this._今の位置 += l_書く分;
                    p_バッファ = p_バッファ[l_書く分..];
                }
            }

            public override void Flush()
            {
            }

            public override int Read(byte[] p_バッファ, int p_開始, int p_長さ) => throw new NotSupportedException();
            public override long Seek(long p_位置, SeekOrigin p_起点) => throw new NotSupportedException();
            public override void SetLength(long p_長さ) => throw new NotSupportedException();

            protected override void Dispose(bool p_Is明示)
            {
                if (!this._Is確定済み)
                {
                    this._Is確定済み = true;
                    if (this._今の位置 > 0)
                    {
                        this.A_塊群.Add(this._今の塊.AsSpan(0, this._今の位置).ToArray());
                    }
                }
                base.Dispose(p_Is明示);
            }
        }

        /// <summary>
        /// 塊の列を先頭から順に読む
        /// </summary>
        /// <param name="p_塊群">読む塊の列</param>
        private sealed class 塊読込ストリーム(List<byte[]> p_塊群) : Stream
        {
            /// <summary>
            /// 読んでいる塊の番号
            /// </summary>
            private int _塊番号;

            /// <summary>
            /// 読んでいる塊の中の位置
            /// </summary>
            private int _位置;

            public override bool CanRead => true;
            public override bool CanSeek => false;
            public override bool CanWrite => false;
            public override long Length => throw new NotSupportedException();
            public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }

            public override int Read(byte[] p_バッファ, int p_開始, int p_長さ)
            {
                return this.Read(p_バッファ.AsSpan(p_開始, p_長さ));
            }

            public override int Read(Span<byte> p_バッファ)
            {
                var l_読んだ = 0;
                while (!p_バッファ.IsEmpty && this._塊番号 < p_塊群.Count)
                {
                    var l_塊 = p_塊群[this._塊番号];
                    if (this._位置 == l_塊.Length)
                    {
                        this._塊番号++;
                        this._位置 = 0;
                        continue;
                    }
                    var l_読む分 = Math.Min(p_バッファ.Length, l_塊.Length - this._位置);
                    l_塊.AsSpan(this._位置, l_読む分).CopyTo(p_バッファ);
                    this._位置 += l_読む分;
                    l_読んだ += l_読む分;
                    p_バッファ = p_バッファ[l_読む分..];
                }
                return l_読んだ;
            }

            public override void Flush()
            {
            }

            public override long Seek(long p_位置, SeekOrigin p_起点) => throw new NotSupportedException();
            public override void SetLength(long p_長さ) => throw new NotSupportedException();
            public override void Write(byte[] p_バッファ, int p_開始, int p_長さ) => throw new NotSupportedException();
        }

        #endregion
    }
}
