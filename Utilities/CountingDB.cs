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
        /// エントリあたりの推定バイト数
        /// </summary>
        private const int エントリあたりの推定バイト数 = 80;

        /// <summary>
        /// IO バッファサイズ
        /// </summary>
        private const int IOバッファサイズ = 1 << 20;

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
        /// まだディスクへ書き出していない k-mer と出現回数
        /// </summary>
        private Dictionary<byte[], ulong> _バッファ;

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
            this._パック長 = (ConfigurationManager.A_実行時引数.A_k長 + 3) / 4;
            var l_総予算 = ConfigurationManager.A_実行時引数.A_メモリ予算バイト数;
            var l_シャードあたりの予算 = l_総予算 / Math.Max(1, p_シャード数);
            this._フラッシュ閾値 = (int)Math.Max(1_024L, Math.Min(int.MaxValue, l_シャードあたりの予算 / エントリあたりの推定バイト数));
            this._バッファ = new Dictionary<byte[], ulong>(this._フラッシュ閾値, this._等価比較器);
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
            }
            else
            {
                this._バッファ[p_パック済みkmer] = 1UL;
                if (this._バッファ.Count >= this._フラッシュ閾値)
                {
                    this.V_フラッシュ();
                }
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

            // フラッシュ済みファイルが 1 件のみだった場合、マージが一度も走らず
            // その元ファイル (_フラッシュ済みファイル に登録済み) がそのまま返される
            // 登録したままだと、この直後に Dispose () が呼ばれた際
            // _フラッシュ済みファイル を掃除する処理で削除されてしまい、
            // 呼び出し元に返したパスが消える (FileNotFoundException の原因)
            // 呼び出し元へ所有権を渡すため、返す前に登録を外しておく
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
        /// 各シャードの CountingDB が Get_統合ファイル () で出力したソート済み・集約済みファイルを、さらにペアワイズマージして 1 本に統合する
        /// </summary>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_ファイル一覧"></param>
        /// <returns></returns>
        public static (string A_ファイルパス, Dictionary<ulong, long>? A_ヒストグラム) Get_統合結果_シャード間(string p_一時ディレクトリ, List<string> p_ファイル一覧)
        {
            var l_比較器 = new ByteArrayComparer();
            var l_パック長 = (ConfigurationManager.A_実行時引数.A_k長 + 3) / 4;
            var l_対象ファイル = new List<string>(p_ファイル一覧);
            var l_接頭辞 = Guid.NewGuid().ToString("N");
            var l_連番 = 1;

            if (l_対象ファイル.Count == 0)
            {
                return (Get_空ファイル(p_一時ディレクトリ, l_接頭辞), []);
            }

            // 最後のマージで書き出される値だけが最終的な出現回数になる
            // そこで集計しておけば、統合ファイルを読み直さずに済む
            Dictionary<ulong, long>? l_ヒストグラム = null;
            while (l_対象ファイル.Count > 1)
            {
                var l_Is最後のマージ = l_対象ファイル.Count == 2;
                if (l_Is最後のマージ)
                {
                    l_ヒストグラム = [];
                }
                var l_出力先 = Path.Combine(p_一時ディレクトリ, $"{l_接頭辞}_workermerge_{l_連番++}");
                V_マージ_2ファイル(l_対象ファイル[0], l_対象ファイル[1], l_出力先, l_パック長, l_比較器, l_ヒストグラム);
                l_対象ファイル.RemoveRange(0, 2);
                l_対象ファイル.Add(l_出力先);
            }

            // マージが一度も走らなかった場合 (シャードが 1 つ) は集計できていない
            return (l_対象ファイル[0], l_ヒストグラム);
        }

        #endregion

        #region 内部メソッド

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
        /// フラッシュ後のファイルは常にソート済み・集約済みであるため、統合側では再集計 (Dictionary への読み直し) が不要になる
        /// </remarks>
        private void V_フラッシュ()
        {
            if (this._バッファ.Count == 0)
            {
                return;
            }

            this._ファイル連番 += 1;
            var l_ファイル名 = Path.Combine(this._一時ディレクトリ, $"{this._ファイル接頭辞}_{this._ファイル連番}");

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
            this._バッファ = new Dictionary<byte[], ulong>(this._フラッシュ閾値, this._等価比較器);
        }

        /// <summary>
        /// ソート済み・集約済みの 2 ファイルを 1 本にマージする
        /// </summary>
        /// <param name="p_ファイル1"></param>
        /// <param name="p_ファイル2"></param>
        /// <param name="p_出力先"></param>
        /// <param name="p_パック長"></param>
        /// <param name="p_比較器"></param>
        /// <param name="p_ヒストグラム"></param>
        /// <remarks>
        /// 同じキーが両方に現れた場合はカウントを合算する<br/>
        /// シャード内統合とシャード間統合で共有する (二重に持つと片方だけ直したときに静かに食い違う)
        /// </remarks>
        private static void V_マージ_2ファイル(string p_ファイル1, string p_ファイル2, string p_出力先, int p_パック長, ByteArrayComparer p_比較器, Dictionary<ulong, long>? p_ヒストグラム = null)
        {
            using (var l_読み込み1 = new BinaryReader(Get_読み込みストリーム(p_ファイル1)))
            {
                using var l_読み込み2 = new BinaryReader(Get_読み込みストリーム(p_ファイル2));
                using var l_書き込み = new BinaryWriter(Get_書き込みストリーム(p_出力先));

                // BinaryReader.ReadBytes は EOF でも null ではなく長さ 0 の配列を
                // 返すため、初回読み取りを保護しないと空ファイルを
                // 「まだ中身がある」と誤認し、続く ReadUInt64 で破綻する
                // k-mer をハッシュでシャードへ振り分けるようにして以降、
                // 空のシャードが普通に発生するようになったため必須
                var l_キー1 = Util.Has続き(l_読み込み1) ? l_読み込み1.ReadBytes(p_パック長) : null;
                var l_キー2 = Util.Has続き(l_読み込み2) ? l_読み込み2.ReadBytes(p_パック長) : null;

                while (l_キー1 != null && l_キー2 != null)
                {
                    var l_比較結果 = p_比較器.Compare(l_キー1, l_キー2);
                    if (l_比較結果 == 0)
                    {
                        l_書き込み.Write(l_キー1);
                        V_書き込み_出現回数(l_書き込み, l_読み込み1.ReadUInt64() + l_読み込み2.ReadUInt64(), p_ヒストグラム);
                        l_キー1 = Util.Has続き(l_読み込み1) ? l_読み込み1.ReadBytes(p_パック長) : null;
                        l_キー2 = Util.Has続き(l_読み込み2) ? l_読み込み2.ReadBytes(p_パック長) : null;
                    }
                    else if (l_比較結果 < 0)
                    {
                        l_書き込み.Write(l_キー1);
                        V_書き込み_出現回数(l_書き込み, l_読み込み1.ReadUInt64(), p_ヒストグラム);
                        l_キー1 = Util.Has続き(l_読み込み1) ? l_読み込み1.ReadBytes(p_パック長) : null;
                    }
                    else
                    {
                        l_書き込み.Write(l_キー2);
                        V_書き込み_出現回数(l_書き込み, l_読み込み2.ReadUInt64(), p_ヒストグラム);
                        l_キー2 = Util.Has続き(l_読み込み2) ? l_読み込み2.ReadBytes(p_パック長) : null;
                    }
                }

                while (l_キー1 != null)
                {
                    l_書き込み.Write(l_キー1);
                    V_書き込み_出現回数(l_書き込み, l_読み込み1.ReadUInt64(), p_ヒストグラム);
                    l_キー1 = Util.Has続き(l_読み込み1) ? l_読み込み1.ReadBytes(p_パック長) : null;
                }

                while (l_キー2 != null)
                {
                    l_書き込み.Write(l_キー2);
                    V_書き込み_出現回数(l_書き込み, l_読み込み2.ReadUInt64(), p_ヒストグラム);
                    l_キー2 = Util.Has続き(l_読み込み2) ? l_読み込み2.ReadBytes(p_パック長) : null;
                }
            }

            File.Delete(p_ファイル1);
            File.Delete(p_ファイル2);
        }

        /// <summary>
        /// 出現回数を書き出し、ヒストグラムが渡されていれば同時に集計する
        /// </summary>
        /// <param name="p_書き込み"></param>
        /// <param name="p_出現回数"></param>
        /// <param name="p_ヒストグラム"></param>
        /// <remarks>
        /// 最終マージの書き出しで集計しておけば、-kc の自動決定のために統合ファイルをもう一度読む必要がなくなる
        /// </remarks>
        private static void V_書き込み_出現回数(BinaryWriter p_書き込み, ulong p_出現回数, Dictionary<ulong, long>? p_ヒストグラム)
        {
            p_書き込み.Write(p_出現回数);
            p_ヒストグラム?[p_出現回数] = p_ヒストグラム.GetValueOrDefault(p_出現回数, 0L) + 1L;
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
