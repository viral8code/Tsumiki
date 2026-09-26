using System.Buffers.Binary;
using System.Runtime.InteropServices;
using Tsumiki.Commons;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 外部ソートで k-mer の出現回数を数えるデータベース
    /// </summary>
    internal class CountingDB : IDisposable
    {
        #region 定数

        /// <summary>
        /// 辞書が倍々に伸びるために、容量が件数に対して膨らみうる倍率
        /// </summary>
        private const int C_容量の膨らみ = 2;

        /// <summary>
        /// IO バッファサイズ
        /// </summary>
        private const int C_IOバッファサイズ = 1 << 20;

        /// <summary>実際の種類数が分かる前に確保する辞書容量の上限</summary>
        private const int C_初期容量の上限 = 16_384;

        /// <summary>
        /// 値を ulong で持つ k の上限
        /// </summary>
        private const int C_小さい値のk上限 = 32;

        /// <summary>
        /// 値を UInt128 で持つ k の上限
        /// </summary>
        private const int C_中くらいの値のk上限 = 64;

        #endregion

        #region 内部変数

        /// <summary>
        /// 比較器
        /// </summary>
        private readonly ByteArrayComparer _比較器;

        /// <summary>
        /// 等価比較器
        /// </summary>
        private readonly ByteArrayEqualityComparer _等価比較器;

        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        /// <summary>
        /// ファイル接頭辞
        /// </summary>
        private readonly string _ファイル接頭辞;

        /// <summary>
        /// k 長
        /// </summary>
        private readonly int _k長;

        /// <summary>
        /// パック長
        /// </summary>
        private readonly int _パック長;

        /// <summary>
        /// フラッシュ閾値
        /// </summary>
        private readonly int _フラッシュ閾値;

        /// <summary>
        /// 次に書き出す一時ファイルの連番
        /// </summary>
        private int _ファイル連番;

        /// <summary>
        /// まだディスクへ書き出していない k-mer と出現回数 (パック済みバイト列で受けたもの)
        /// </summary>
        private Dictionary<byte[], ulong> _バッファ;

        /// <summary>
        /// まだディスクへ書き出していない k-mer と出現回数 (128 塩基までのパック値で受けたもの)
        /// </summary>
        private Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong> _値バッファ;

        /// <summary>
        /// まだディスクへ書き出していない k-mer と出現回数 (k &lt;= 32 の値)
        /// </summary>
        private Dictionary<ulong, ulong> _値バッファ_小;

        /// <summary>
        /// まだディスクへ書き出していない k-mer と出現回数 (32 &lt; k &lt;= 64 の値)
        /// </summary>
        private Dictionary<UInt128, ulong> _値バッファ_中;

        /// <summary>
        /// フラッシュ済みファイル
        /// </summary>
        private readonly List<string> _フラッシュ済みファイル = [];

        #endregion

        #region コンストラクタ

        /// <summary>
        /// p_シャード数 には、同時に生きている CountingDB の総数を渡す
        /// </summary>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_シャード数"></param>
        public CountingDB(string p_一時ディレクトリ, int p_シャード数 = 1)
        {
            this._ファイル接頭辞 = Guid.NewGuid().ToString("N");
            this._比較器 = new();
            this._等価比較器 = new();
            this._一時ディレクトリ = p_一時ディレクトリ;
            this._k長 = ConfigurationManager.A_実行時引数.A_k長;
            this._パック長 = (this._k長 + 3) / 4;
            var l_総予算 = ConfigurationManager.A_実行時引数.A_メモリ予算バイト数;
            var l_シャードあたりの予算 = l_総予算 / Math.Max(1, p_シャード数);
            this._フラッシュ閾値 = (int)Math.Max(1_024L, Math.Min(int.MaxValue, l_シャードあたりの予算 / Get_エントリあたりのバイト数(this._k長)));
            this._バッファ = new Dictionary<byte[], ulong>(Math.Min(this._フラッシュ閾値, C_初期容量の上限), this._等価比較器);
            this._値バッファ = [];
            this._値バッファ_小 = [];
            this._値バッファ_中 = [];
            this.V_用意_値バッファ();
            this._ファイル連番 = 0;
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// k-mer を 1 件数える
        /// </summary>
        /// <param name="p_kmer">数える k-mer</param>
        public void V_登録(Span<byte> p_kmer)
        {
            var l_パック済み = new byte[(p_kmer.Length + 3) / 4];
            var l_書き込み位置 = 0;
            for (var i = 0; i < p_kmer.Length; i += 4)
            {
                var l_バイト = 0;
                for (var j = 0; j < 4; j++)
                {
                    var l_塩基ID = i + j < p_kmer.Length ? p_kmer[i + j] : Consts.塩基ID.A;
                    l_バイト <<= 2;
                    l_バイト |= l_塩基ID - 1;
                }

                l_パック済み[l_書き込み位置++] = (byte)l_バイト;
            }

            this.V_登録_パック済み(l_パック済み);
        }

        /// <summary>
        /// k-mer を 1 件登録する
        /// </summary>
        /// <param name="p_パック済みkmer"></param>
        public void V_登録_パック済み(byte[] p_パック済みkmer)
        {
            if (this._バッファ.TryGetValue(p_パック済みkmer, out var l_出現回数))
            {
                this._バッファ[p_パック済みkmer] = l_出現回数 + 1UL;
                return;
            }

            this._バッファ[p_パック済みkmer] = 1UL;
            if (this.Is閾値到達())
            {
                this.V_フラッシュ();
            }
        }

        /// <summary>
        /// 右詰めのパック値で表した k-mer を 1 件登録する (k &lt;= 128)
        /// </summary>
        /// <param name="p_値"></param>
        public void V_登録_値((UInt128 A_上位, UInt128 A_下位) p_値)
        {
            bool l_Is既存;
            if (this._k長 <= C_小さい値のk上限)
            {
                CollectionsMarshal.GetValueRefOrAddDefault(this._値バッファ_小, (ulong)p_値.A_下位, out l_Is既存)++;
            }
            else if (this._k長 <= C_中くらいの値のk上限)
            {
                CollectionsMarshal.GetValueRefOrAddDefault(this._値バッファ_中, p_値.A_下位, out l_Is既存)++;
            }
            else
            {
                CollectionsMarshal.GetValueRefOrAddDefault(this._値バッファ, p_値, out l_Is既存)++;
            }

            if (!l_Is既存 && this.Is閾値到達())
            {
                this.V_フラッシュ();
            }
        }

        /// <summary>
        /// このシャードのフラッシュ済みファイルをすべて 1 本にマージし、そのパスを返す
        /// </summary>
        /// <returns></returns>
        public string Get_統合ファイル()
        {
            this.V_フラッシュ();

            var l_対象ファイル = new List<string>(this._フラッシュ済みファイル);

            if (l_対象ファイル.Count == 0)
            {
                return Get_空ファイル(this._一時ディレクトリ, this._ファイル接頭辞);
            }

            var l_連番 = this._ファイル連番 + 1;
            while (l_対象ファイル.Count > 1)
            {
                var l_出力先 = Path.Combine(this._一時ディレクトリ, $"{this._ファイル接頭辞}_merged_{l_連番++}");
                V_マージ_2ファイル(l_対象ファイル[0], l_対象ファイル[1], l_出力先, this._パック長);
                l_対象ファイル.RemoveRange(0, 2);
                l_対象ファイル.Add(l_出力先);
            }

            var l_最終ファイル = l_対象ファイル[0];
            _ = this._フラッシュ済みファイル.Remove(l_最終ファイル);
            return l_最終ファイル;
        }

        /// <summary>
        /// 保持している資源を解放する
        /// </summary>
        public void Dispose()
        {
            foreach (var l_ファイル in this._フラッシュ済みファイル)
            {
                中間データ置き場.V_削除(l_ファイル);
            }
        }

        /// <summary>
        /// 右詰めのパック値を、ファイル上の並び (2 bit/塩基・先頭塩基が最上位・余りビットは下位) のバイト列にする
        /// </summary>
        /// <param name="p_値"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_出力">パック長ぶんの書き込み先</param>
        internal static void V_変換_パック済みバイト列((UInt128 A_上位, UInt128 A_下位) p_値, int p_k長, Span<byte> p_出力)
        {
            var l_余りビット = (8 * p_出力.Length) - (2 * p_k長);
            var (l_上位, l_下位) = p_値;
            if (l_余りビット > 0)
            {
                l_上位 = (l_上位 << l_余りビット) | (l_下位 >> (128 - l_余りビット));
                l_下位 <<= l_余りビット;
            }

            for (var i = p_出力.Length - 1; i >= 0; i--)
            {
                p_出力[i] = (byte)l_下位;
                l_下位 = (l_下位 >> 8) | (l_上位 << 120);
                l_上位 >>= 8;
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// メモリ上に 1 件置くのに要るバイト数
        /// </summary>
        /// <param name="p_k長">k 長</param>
        /// <returns></returns>
        private static int Get_エントリあたりのバイト数(int p_k長)
        {
            if (p_k長 <= TrustedKmerIndex.C_パック値のk上限)
            {
                const int l_容量に比例する分 = (32 + 8 + 4 + 4) + 4;

                const int l_一時配列 = 32 + 8;
                return (l_容量に比例する分 * C_容量の膨らみ) + l_一時配列;
            }

            const int l_辞書の分 = (8 + 8 + 4 + 4) + 4;

            var l_鍵の実体 = 24 + (((((p_k長 + 3) / 4) + 7) / 8) * 8);

            const int l_一時配列_大 = 16;
            return (l_辞書の分 * C_容量の膨らみ) + l_鍵の実体 + l_一時配列_大;
        }

        /// <summary>
        /// メモリ上の件数がフラッシュ閾値に達したか
        /// </summary>
        /// <returns></returns>
        private bool Is閾値到達()
        {
            return this._バッファ.Count + this._値バッファ.Count + this._値バッファ_小.Count + this._値バッファ_中.Count >= this._フラッシュ閾値;
        }

        /// <summary>
        /// k に合った幅の値バッファだけを、空で用意する
        /// </summary>
        private void V_用意_値バッファ()
        {
            var l_容量 = Math.Min(this._フラッシュ閾値, C_初期容量の上限);
            this._値バッファ_小 = this._k長 <= C_小さい値のk上限 ? new Dictionary<ulong, ulong>(l_容量) : [];
            this._値バッファ_中 = this._k長 is > C_小さい値のk上限 and <= C_中くらいの値のk上限 ? new Dictionary<UInt128, ulong>(l_容量) : [];
            this._値バッファ = this._k長 > C_中くらいの値のk上限 ? new Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>(l_容量) : [];
        }

        /// <summary>
        /// 値と出現回数を値の昇順に並べて、パック済みバイト列でファイルへ書き出す
        /// </summary>
        /// <param name="p_バッファ">書き出す値と出現回数</param>
        /// <param name="p_変換">値を右詰めのパック値に直す</param>
        private void V_書出_値バッファ<T>(Dictionary<T, ulong> p_バッファ, Func<T, (UInt128 A_上位, UInt128 A_下位)> p_変換) where T : notnull
        {
            var l_ファイル名 = this.Get_次のファイル名();
            var l_キー = new T[p_バッファ.Count];
            var l_回数 = new ulong[l_キー.Length];
            var l_位置 = 0;
            foreach (var (l_値, l_出現回数) in p_バッファ)
            {
                l_キー[l_位置] = l_値;
                l_回数[l_位置++] = l_出現回数;
            }

            Array.Sort(l_キー, l_回数);

            var l_バイト列 = new byte[this._パック長];
            using (var l_書き込み = new BinaryWriter(Get_書き込みストリーム(l_ファイル名)))
            {
                for (var i = 0; i < l_キー.Length; i++)
                {
                    V_変換_パック済みバイト列(p_変換(l_キー[i]), this._k長, l_バイト列);
                    l_書き込み.Write(l_バイト列);
                    l_書き込み.Write(l_回数[i]);
                }
            }

            this._フラッシュ済みファイル.Add(l_ファイル名);
        }

        /// <summary>
        /// 書き込み用のストリームを開いて返す
        /// </summary>
        /// <param name="p_ファイル名">開くファイル名</param>
        /// <returns>書き込み用のストリーム</returns>
        private static Stream Get_書き込みストリーム(string p_ファイル名)
        {
            return 中間データ置き場.A_Is有効
                ? 中間データ置き場.Get_書込ストリーム(p_ファイル名)
                : new FileStream(p_ファイル名, FileMode.Create, FileAccess.Write, FileShare.None, C_IOバッファサイズ, FileOptions.SequentialScan);
        }

        /// <summary>
        /// 読み込み用のストリームを開いて返す
        /// </summary>
        /// <param name="p_ファイル名">開くファイル名</param>
        /// <returns>読み込み用のストリーム</returns>
        private static Stream Get_読み込みストリーム(string p_ファイル名)
        {
            return 中間データ置き場.Get_読込ストリーム(p_ファイル名);
        }

        /// <summary>
        /// メモリ上の集約済みカウントをキー順にソートしてディスクへ書き出す
        /// </summary>
        private void V_フラッシュ()
        {
            if (this._バッファ.Count > 0)
            {
                var l_ファイル名 = this.Get_次のファイル名();
                var l_エントリ = this._バッファ.ToArray();
                Array.Sort(l_エントリ, (x, y) => this._比較器.Compare(x.Key, y.Key));

                using (var l_書き込み = new BinaryWriter(Get_書き込みストリーム(l_ファイル名)))
                {
                    foreach (var l_項目 in l_エントリ)
                    {
                        l_書き込み.Write(l_項目.Key);
                        l_書き込み.Write(l_項目.Value);
                    }
                }

                this._フラッシュ済みファイル.Add(l_ファイル名);
                this._バッファ = new Dictionary<byte[], ulong>(Math.Min(this._フラッシュ閾値, C_初期容量の上限), this._等価比較器);
            }

            if (this._値バッファ.Count + this._値バッファ_小.Count + this._値バッファ_中.Count > 0)
            {
                if (this._値バッファ_小.Count > 0)
                {
                    this.V_書出_値バッファ(this._値バッファ_小, static x => ((UInt128)0, (UInt128)x));
                }

                if (this._値バッファ_中.Count > 0)
                {
                    this.V_書出_値バッファ(this._値バッファ_中, static x => ((UInt128)0, x));
                }

                if (this._値バッファ.Count > 0)
                {
                    this.V_書出_値バッファ(this._値バッファ, static x => x);
                }

                this.V_用意_値バッファ();
            }
        }

        /// <summary>
        /// 次に書き出す一時ファイルのパス
        /// </summary>
        /// <returns></returns>
        private string Get_次のファイル名()
        {
            this._ファイル連番 += 1;
            return Path.Combine(this._一時ディレクトリ, $"{this._ファイル接頭辞}_{this._ファイル連番}");
        }

        /// <summary>
        /// ソート済み・集約済みの 2 ファイルを 1 本にマージする
        /// </summary>
        /// <param name="p_ファイル1"></param>
        /// <param name="p_ファイル2"></param>
        /// <param name="p_出力先"></param>
        /// <param name="p_パック長"></param>
        private static void V_マージ_2ファイル(string p_ファイル1, string p_ファイル2, string p_出力先, int p_パック長)
        {
            using (var l_読み込み1 = new エントリ読み込み(Get_読み込みストリーム(p_ファイル1), p_パック長))
            {
                using var l_読み込み2 = new エントリ読み込み(Get_読み込みストリーム(p_ファイル2), p_パック長);
                using var l_書き込み = new BufferedStream(Get_書き込みストリーム(p_出力先), C_IOバッファサイズ);
                while (l_読み込み1.Has項目 && l_読み込み2.Has項目)
                {
                    var l_比較結果 = l_読み込み1.A_キー.SequenceCompareTo(l_読み込み2.A_キー);
                    if (l_比較結果 == 0)
                    {
                        V_書込_エントリ(l_書き込み, l_読み込み1.A_キー, l_読み込み1.A_出現回数 + l_読み込み2.A_出現回数);
                        l_読み込み1.V_進む();
                        l_読み込み2.V_進む();
                    }
                    else if (l_比較結果 < 0)
                    {
                        V_書込_エントリ(l_書き込み, l_読み込み1.A_キー, l_読み込み1.A_出現回数);
                        l_読み込み1.V_進む();
                    }
                    else
                    {
                        V_書込_エントリ(l_書き込み, l_読み込み2.A_キー, l_読み込み2.A_出現回数);
                        l_読み込み2.V_進む();
                    }
                }

                foreach (var l_残り in new[] { l_読み込み1, l_読み込み2 })
                {
                    while (l_残り.Has項目)
                    {
                        V_書込_エントリ(l_書き込み, l_残り.A_キー, l_残り.A_出現回数);
                        l_残り.V_進む();
                    }
                }
            }

            中間データ置き場.V_削除(p_ファイル1);
            中間データ置き場.V_削除(p_ファイル2);
        }

        /// <summary>
        /// k-mer と出現回数の 1 件を書く
        /// </summary>
        /// <param name="p_書き込み">書き込み先</param>
        /// <param name="p_キー">パック済みの k-mer</param>
        /// <param name="p_出現回数">出現回数</param>
        private static void V_書込_エントリ(Stream p_書き込み, ReadOnlySpan<byte> p_キー, ulong p_出現回数)
        {
            Span<byte> l_回数 = stackalloc byte[sizeof(ulong)];
            BinaryPrimitives.WriteUInt64LittleEndian(l_回数, p_出現回数);
            p_書き込み.Write(p_キー);
            p_書き込み.Write(l_回数);
        }

        /// <summary>
        /// 空のソート済みファイルを作って、そのパスを返す
        /// </summary>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_接頭辞"></param>
        /// <returns></returns>
        private static string Get_空ファイル(string p_一時ディレクトリ, string p_接頭辞)
        {
            var l_ファイル名 = Path.Combine(p_一時ディレクトリ, $"{p_接頭辞}_empty");
            using (Get_書き込みストリーム(l_ファイル名))
            {
            }

            return l_ファイル名;
        }

        #endregion

    }
}
