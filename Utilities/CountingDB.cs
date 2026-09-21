using Tsumiki.Commons;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 外部ソートで k-mer の出現回数を数えるデータベース
    /// </summary>
    /// <remarks>
    /// メモリ上の Dictionary で集約しつつ、閾値を超えたらソート済みファイルへフラッシュし、最後にペアワイズマージして 1 本の整列済みファイルへ統合する
    /// </remarks>
    internal class CountingDB : IDisposable
    {
        #region 定数

        /// <summary>
        /// 辞書が倍々に伸びるために、容量が件数に対して膨らみうる倍率
        /// </summary>
        private const int 容量の膨らみ = 2;

        /// <summary>
        /// IO バッファサイズ
        /// </summary>
        private const int IOバッファサイズ = 1 << 20;

        /// <summary>実際の種類数が分かる前に確保する辞書容量の上限</summary>
        private const int 初期容量の上限 = 16_384;

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
        /// <remarks>
        /// メモリ予算を等分するために使う
        /// </remarks>
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
            this._バッファ = new Dictionary<byte[], ulong>(Math.Min(this._フラッシュ閾値, 初期容量の上限), this._等価比較器);
            this._値バッファ = new Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>(Math.Min(this._フラッシュ閾値, 初期容量の上限));
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
        /// <remarks>
        /// バイト列のキーは k-mer ごとに配列を確保し、比較もバイト単位になる
        /// </remarks>
        public void V_登録_値((UInt128 A_上位, UInt128 A_下位) p_値)
        {
            if (this._値バッファ.TryGetValue(p_値, out var l_出現回数))
            {
                this._値バッファ[p_値] = l_出現回数 + 1UL;
                return;
            }

            this._値バッファ[p_値] = 1UL;
            if (this.Is閾値到達())
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
            // メモリ上に残っている未フラッシュ分を書き出す
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
                V_マージ_2ファイル(l_対象ファイル[0], l_対象ファイル[1], l_出力先, this._パック長, this._比較器);
                l_対象ファイル.RemoveRange(0, 2);
                l_対象ファイル.Add(l_出力先);
            }

            // 1 件のみでマージが走らなかった場合、登録したままだと Dispose で消されてしまうため、所有権を呼び出し元へ渡す
            var l_最終ファイル = l_対象ファイル[0];
            _ = this._フラッシュ済みファイル.Remove(l_最終ファイル);
            return l_最終ファイル;
        }

        /// <summary>
        /// 保持している資源を解放する
        /// </summary>
        public void Dispose()
        {
            // 未フラッシュのデータは統合側で処理される想定だが、
            // 統合を呼ばずに破棄された場合に備えて残存ファイルを掃除する
            foreach (var l_ファイル in this._フラッシュ済みファイル)
            {
                if (File.Exists(l_ファイル))
                {
                    File.Delete(l_ファイル);
                }
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
        /// <remarks>
        /// 辞書の実体だけでなく、容量の膨らみと、フラッシュで辞書と同時に生きる整列用の一時配列まで数える<br/>
        /// ここを小さく見積もると -mem の指定より実際の常駐がずっと大きくなる
        /// </remarks>
        /// <returns></returns>
        private static int Get_エントリあたりのバイト数(int p_k長)
        {
            if (p_k長 <= TrustedKmerIndex.パック値のk上限)
            {
                // 鍵 32 + 回数 8 + ハッシュ 4 + 次 4、バケット 4
                const int l_容量に比例する分 = (32 + 8 + 4 + 4) + 4;

                // 整列用の鍵配列 32 と回数配列 8
                const int l_一時配列 = 32 + 8;
                return (l_容量に比例する分 * 容量の膨らみ) + l_一時配列;
            }

            // 参照 8 + 回数 8 + ハッシュ 4 + 次 4、バケット 4
            const int l_辞書の分 = (8 + 8 + 4 + 4) + 4;

            // パック済みバイト列のオブジェクト (ヘッダ 24 + 8 バイト境界へ丸めた本体)
            var l_鍵の実体 = 24 + (((((p_k長 + 3) / 4) + 7) / 8) * 8);

            // 整列用に取り出す KeyValuePair の配列
            const int l_一時配列_大 = 16;
            return (l_辞書の分 * 容量の膨らみ) + l_鍵の実体 + l_一時配列_大;
        }

        /// <summary>
        /// メモリ上の件数がフラッシュ閾値に達したか
        /// </summary>
        /// <returns></returns>
        private bool Is閾値到達()
        {
            return this._バッファ.Count + this._値バッファ.Count >= this._フラッシュ閾値;
        }

        /// <summary>
        /// 書き込み用のストリームを開いて返す
        /// </summary>
        /// <param name="p_ファイル名">開くファイル名</param>
        /// <returns>書き込み用のストリーム</returns>
        private static FileStream Get_書き込みストリーム(string p_ファイル名)
        {
            return new FileStream(p_ファイル名, FileMode.Create, FileAccess.Write, FileShare.None, IOバッファサイズ, FileOptions.SequentialScan);
        }

        /// <summary>
        /// 読み込み用のストリームを開いて返す
        /// </summary>
        /// <param name="p_ファイル名">開くファイル名</param>
        /// <returns>読み込み用のストリーム</returns>
        private static FileStream Get_読み込みストリーム(string p_ファイル名)
        {
            return new FileStream(p_ファイル名, FileMode.Open, FileAccess.Read, FileShare.Read, IOバッファサイズ, FileOptions.SequentialScan);
        }

        /// <summary>
        /// メモリ上の集約済みカウントをキー順にソートしてディスクへ書き出す
        /// </summary>
        /// <remarks>
        /// フラッシュ後のファイルは常にソート済み・集約済みであるため、統合側では再集計 (Dictionary への読み直し) が不要になる<br/>
        /// 2 種類のバッファは同じ並び規則の別ファイルとして書き、同じキーが両方にあってもマージで合算される
        /// </remarks>
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
                this._バッファ = new Dictionary<byte[], ulong>(Math.Min(this._フラッシュ閾値, 初期容量の上限), this._等価比較器);
            }

            if (this._値バッファ.Count > 0)
            {
                var l_ファイル名 = this.Get_次のファイル名();

                // 右詰めのパック値の大小は、ファイル上のバイト列の辞書順と一致する
                var l_キー = new (UInt128 A_上位, UInt128 A_下位)[this._値バッファ.Count];
                var l_回数 = new ulong[l_キー.Length];
                var l_位置 = 0;
                foreach (var (l_値, l_出現回数) in this._値バッファ)
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
                        V_変換_パック済みバイト列(l_キー[i], this._k長, l_バイト列);
                        l_書き込み.Write(l_バイト列);
                        l_書き込み.Write(l_回数[i]);
                    }
                }

                this._フラッシュ済みファイル.Add(l_ファイル名);
                this._値バッファ = new Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>(Math.Min(this._フラッシュ閾値, 初期容量の上限));
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
        /// <param name="p_比較器"></param>
        /// <remarks>
        /// 同じキーが両方に現れた場合はカウントを合算する
        /// </remarks>
        private static void V_マージ_2ファイル(string p_ファイル1, string p_ファイル2, string p_出力先, int p_パック長, ByteArrayComparer p_比較器)
        {
            using (var l_読み込み1 = new BinaryReader(Get_読み込みストリーム(p_ファイル1)))
            {
                using var l_読み込み2 = new BinaryReader(Get_読み込みストリーム(p_ファイル2));
                using var l_書き込み = new BinaryWriter(Get_書き込みストリーム(p_出力先));

                // BinaryReader.ReadBytes は EOF でも長さ 0 の配列を返すため、空ファイルを中身があると誤認しないよう先に確かめる
                var l_キー1 = Util.Has続き(l_読み込み1) ? l_読み込み1.ReadBytes(p_パック長) : null;
                var l_キー2 = Util.Has続き(l_読み込み2) ? l_読み込み2.ReadBytes(p_パック長) : null;

                while (l_キー1 != null && l_キー2 != null)
                {
                    var l_比較結果 = p_比較器.Compare(l_キー1, l_キー2);
                    if (l_比較結果 == 0)
                    {
                        l_書き込み.Write(l_キー1);
                        l_書き込み.Write(l_読み込み1.ReadUInt64() + l_読み込み2.ReadUInt64());
                        l_キー1 = Util.Has続き(l_読み込み1) ? l_読み込み1.ReadBytes(p_パック長) : null;
                        l_キー2 = Util.Has続き(l_読み込み2) ? l_読み込み2.ReadBytes(p_パック長) : null;
                    }
                    else if (l_比較結果 < 0)
                    {
                        l_書き込み.Write(l_キー1);
                        l_書き込み.Write(l_読み込み1.ReadUInt64());
                        l_キー1 = Util.Has続き(l_読み込み1) ? l_読み込み1.ReadBytes(p_パック長) : null;
                    }
                    else
                    {
                        l_書き込み.Write(l_キー2);
                        l_書き込み.Write(l_読み込み2.ReadUInt64());
                        l_キー2 = Util.Has続き(l_読み込み2) ? l_読み込み2.ReadBytes(p_パック長) : null;
                    }
                }

                while (l_キー1 != null)
                {
                    l_書き込み.Write(l_キー1);
                    l_書き込み.Write(l_読み込み1.ReadUInt64());
                    l_キー1 = Util.Has続き(l_読み込み1) ? l_読み込み1.ReadBytes(p_パック長) : null;
                }

                while (l_キー2 != null)
                {
                    l_書き込み.Write(l_キー2);
                    l_書き込み.Write(l_読み込み2.ReadUInt64());
                    l_キー2 = Util.Has続き(l_読み込み2) ? l_読み込み2.ReadBytes(p_パック長) : null;
                }
            }

            File.Delete(p_ファイル1);
            File.Delete(p_ファイル2);
        }

        /// <summary>
        /// 空のソート済みファイルを作って、そのパスを返す
        /// </summary>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_接頭辞"></param>
        /// <remarks>
        /// 登録が 1 件も無かったシャードでも、統合処理に渡せる形を保つために使う
        /// </remarks>
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
