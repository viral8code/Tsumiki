using Tsumiki.Common;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer の出現回数カウントと、カットオフを通過した信頼できる k-mer の
    /// 厳密な集合を保持する
    /// </summary>
    /// <remarks>
    /// Bloom filter のような近似判定は使わない<br/>
    /// フォールスポジティブによる
    /// グラフ構造の誤判定 (偽の分岐点・偽の隣接) を原理的に排除できないため<br/>
    /// 細菌ゲノム規模なら厳密な集合をメモリに載せられる
    /// </remarks>
    internal class TrustedKmerIndex : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        // k-mer カウント用のシャード
        // k-mer 自身のハッシュで振り分けるため、
        // ある k-mer は必ず 1 つのシャードにしか載らない

        /// <summary>
        /// ワーカーごとの k-mer カウンタ
        /// </summary>
        private CountingDB[]? _カウンタ群;

        // シャードごとのロック
        // k-mer をワーカー単位ではなくハッシュ値で
        // 振り分けるようにしたため、複数スレッドが同じシャードへ書きうる

        /// <summary>
        /// シャードロック
        /// </summary>
        private readonly object[]? _シャードロック;

        // カットオフを通過した k-mer の厳密な集合 (常に正規形)
        // 値はカバレッジ
        //
        // k の範囲で経路を分けるのは、ulong / UInt128 が値型でヒープ確保を
        // 伴わないため
        // KmerKey は毎回 ulong[] を確保し、1 リードあたり
        // 数百〜数千回呼ばれる所属判定では実データ規模で致命的に効く
        // k>64 のときだけ KmerKey へフォールバックする

        /// <summary>
        /// 信頼できる k-mer と出現回数 (k &gt; 64)
        /// </summary>
        private Dictionary<KmerKey, ulong>? _信頼kmer_大;

        /// <summary>
        /// 信頼できる k-mer と出現回数 (k &lt;= 32)
        /// </summary>
        private Dictionary<ulong, ulong>? _信頼kmer_小;

        // 33 <= k <= 64 用
        // 150 bp リードで k=31 のままだと 31 bp 以上の反復配列が
        // すべて潰れてしまい contig N50 が伸びないため、k を 63 前後まで上げられる
        // ことが品質上きわめて重要になる

        /// <summary>
        /// 信頼できる k-mer と出現回数 (33 &lt;= k &lt;= 64)
        /// </summary>
        private Dictionary<UInt128, ulong>? _信頼kmer_中;

        // 全シャードを 1 本にマージしたソート済みファイル
        // 統合は高くつくため
        // 一度だけ行い、ヒストグラムの集計とカットオフの適用で使い回す

        /// <summary>
        /// ワーカーごとのカウントを 1 本へ統合したファイルのパス
        /// </summary>
        private string? _統合ファイルパス;

        // 最終マージの書き出し中に集計したヒストグラム
        // シャードが 1 つで
        // マージが走らなかった場合は null になり、そのときだけ読み直す

        /// <summary>
        /// 統合の際に数えた、出現回数ごとの k-mer 種類数
        /// </summary>
        private Dictionary<ulong, long>? _統合時のヒストグラム;

        /// <summary>
        /// 直近の V_カットオフ で集計した出現回数ヒストグラム
        /// </summary>
        /// <remarks>
        /// カットオフ判定と同じループで作れるため追加のコストはかからない<br/>
        /// ゲノムサイズやカバレッジの推定に使う
        /// </remarks>
        public IReadOnlyDictionary<ulong, long> A_出現回数ヒストグラム { get; private set; }
            = new Dictionary<ulong, long>();

        // k 長は構築時に固定する
        // グローバルから毎回読むと、別の k を使う処理が
        // 走った後にこのインデックスへ問い合わせたとき、内部表現と食い違う経路を
        // 選んで破綻する (multi-k のように k が切り替わる場面で実際に起きる)

        /// <summary>
        /// k 長
        /// </summary>
        private readonly int _k長;

        /// <summary>
        /// 小経路を使うか
        /// </summary>
        private bool 小経路を使うか => this._k長 <= 32;

        /// <summary>
        /// 中経路を使うか
        /// </summary>
        private bool 中経路を使うか => this._k長 is > 32 and <= 64;

        public TrustedKmerIndex(string p_一時ディレクトリ)
        {
            this._k長 = ConfigurationManager.A_実行時引数.A_k長;
            this._一時ディレクトリ = p_一時ディレクトリ;
            var l_シャード数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            this._カウンタ群 = new CountingDB[l_シャード数];
            this._シャードロック = new object[l_シャード数];
            for (var i = 0; i < l_シャード数; i++)
            {
                this._カウンタ群[i] = new CountingDB(p_一時ディレクトリ, l_シャード数);
                this._シャードロック[i] = new object();
            }
        }

        /// <summary>
        /// 曖昧塩基を含む k-mer を、ありうる塩基の組み合わせすべてに展開して登録する
        /// </summary>
        /// <remarks>
        /// 展開は塩基 ID の空間で行い、1 件ずつ通常の登録へ渡す<br/>
        /// パック済みバイト列を
        /// 自前で組み立てると正規化とシャード振り分けが通常経路とずれる
        /// </remarks>
        public void V_登録_曖昧塩基あり(Span<byte[]> p_塩基候補列, int p_ワーカー番号)
        {
            if (this._カウンタ群 is null)
            {
                return;
            }
            var l_kmer = new byte[p_塩基候補列.Length];
            this.V_登録_組み合わせ展開(p_塩基候補列, 0, l_kmer, p_ワーカー番号);
        }

        /// <summary>
        /// 曖昧塩基の候補をすべて展開して登録する
        /// </summary>
        /// <param name="p_塩基候補列">位置ごとの塩基候補</param>
        /// <param name="p_位置">いま決めている位置</param>
        /// <param name="p_kmer">組み立て中の k-mer</param>
        /// <param name="p_ワーカー番号">登録先のワーカー</param>
        private void V_登録_組み合わせ展開(Span<byte[]> p_塩基候補列, int p_位置, byte[] p_kmer, int p_ワーカー番号)
        {
            if (p_位置 == p_塩基候補列.Length)
            {
                this.V_登録(p_kmer.AsSpan(), p_ワーカー番号);
                return;
            }
            foreach (var l_塩基ID in p_塩基候補列[p_位置])
            {
                p_kmer[p_位置] = l_塩基ID;
                this.V_登録_組み合わせ展開(p_塩基候補列, p_位置 + 1, p_kmer, p_ワーカー番号);
            }
        }

        /// <summary>
        /// k-mer を 1 件カウントする
        /// </summary>
        /// <remarks>
        /// 振り分けはワーカー番号ではなく k-mer 自身のハッシュで行う<br/>
        /// ワーカー単位だと同じ k-mer がスレッド数ぶんの辞書に重複して載る<br/>
        /// 数える前に正規形へ寄せるのも同じ理由で、両向きを別キーにすると
        /// エントリ数と書き出し量が倍になる
        /// </remarks>
        public void V_登録(Span<byte> p_kmer, int p_ワーカー番号)
        {
            if (this._カウンタ群 is not { } l_カウンタ群)
            {
                return;
            }

            var l_パック済み = Get_正規化パック(p_kmer);
            var l_シャード = (int)(Get_ハッシュ(l_パック済み) % (uint)l_カウンタ群.Length);
            lock (this._シャードロック![l_シャード])
            {
                l_カウンタ群[l_シャード].V_登録_パック済み(l_パック済み);
            }
        }

        /// <summary>
        /// k-mer を正規形の向きで 2 bit パックする
        /// </summary>
        /// <remarks>
        /// 先頭塩基が上位ビットに来るためバイト列の辞書順が塩基列の辞書順と
        /// 一致し、外部マージソートの順序と整合する
        /// </remarks>
        private static byte[] Get_正規化パック(ReadOnlySpan<byte> p_kmer)
        {
            var l_順鎖を使うか = Get_順鎖が正規形か(p_kmer);
            var l_パック済み = new byte[(p_kmer.Length + 3) / 4];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                // 逆鎖側を採用する場合は、末尾から相補塩基を取り出す
                // 相補は A(1)<->T(4), C(2)<->G(3) なので 5 - x で得られる
                var l_塩基ID = l_順鎖を使うか ? p_kmer[i] : (byte)(5 - p_kmer[p_kmer.Length - 1 - i]);
                l_パック済み[i >> 2] |= (byte)((l_塩基ID - 1) << ((3 - (i & 3)) << 1));
            }
            return l_パック済み;
        }

        /// <summary>
        /// 順鎖側がその逆相補以下 (辞書順) かどうか
        /// </summary>
        /// <remarks>
        /// 確保なしで判定する
        /// </remarks>
        private static bool Get_順鎖が正規形か(ReadOnlySpan<byte> p_kmer)
        {
            int i = 0, j = p_kmer.Length - 1;
            while (i <= j)
            {
                var l_順鎖 = p_kmer[i];
                var l_逆鎖 = (byte)(5 - p_kmer[j]);
                if (l_順鎖 != l_逆鎖)
                {
                    return l_順鎖 < l_逆鎖;
                }
                i++;
                j--;
            }
            // 回文 (自身が逆相補と一致)
            // どちらでも同じなので順鎖扱い
            return true;
        }

        /// <summary>
        /// パック済みキーの FNV-1 a ハッシュ
        /// </summary>
        /// <remarks>
        /// シャードの振り分けに使う
        /// </remarks>
        private static uint Get_ハッシュ(byte[] p_パック済みkmer)
        {
            var l_ハッシュ = 2166136261U;
            foreach (var l_バイト in p_パック済みkmer)
            {
                l_ハッシュ ^= l_バイト;
                l_ハッシュ *= 16777619U;
            }
            return l_ハッシュ;
        }

        /// <summary>
        /// kmer(順鎖・逆鎖いずれの向きでもよい) がカットオフを通過した
        /// 信頼できる k-mer 集合に含まれるかどうかを厳密に判定する
        /// </summary>
        public bool Get_含まれるか(Span<byte> p_kmer)
        {
            return this.小経路を使うか
                ? this._信頼kmer_小!.ContainsKey(Get_正規形_小(p_kmer))
                : this.中経路を使うか ? this._信頼kmer_中!.ContainsKey(Get_正規形_中(p_kmer)) : this._信頼kmer_大!.ContainsKey(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// 正規形にパック済みの値で所属を判定する (k &lt;= 32)
        /// </summary>
        /// <remarks>
        /// walk は 1 塩基ずつ進むためパック値を転がして更新できる<br/>
        /// Span を受ける版だと呼ぶたびに O(k) の詰め直しが要る
        /// </remarks>
        public bool Get_含まれるか_小(ulong p_正規形)
        {
            return this._信頼kmer_小!.ContainsKey(p_正規形);
        }

        /// <summary>
        /// 正規形にパック済みの値で所属を判定する (33 &lt;= k &lt;= 64)
        /// </summary>
        public bool Get_含まれるか_中(UInt128 p_正規形)
        {
            return this._信頼kmer_中!.ContainsKey(p_正規形);
        }

        /// <summary>
        /// kmerの出現回数 (カバレッジ) を返す
        /// </summary>
        /// <remarks>
        /// 信頼できる k-mer 集合に
        /// 含まれない場合は 0 を返す
        /// </remarks>
        public ulong Get_カバレッジ(Span<byte> p_kmer)
        {
            return this.小経路を使うか
                ? this._信頼kmer_小!.GetValueOrDefault(Get_正規形_小(p_kmer), 0UL)
                : this.中経路を使うか
                ? this._信頼kmer_中!.GetValueOrDefault(Get_正規形_中(p_kmer), 0UL)
                : this._信頼kmer_大!.GetValueOrDefault(new KmerKey(p_kmer).Get_正規形(), 0UL);
        }

        /// <summary>
        /// kmer(塩基ID 1-4、長さ 32 以下) を2 bit/塩基で ulong1個にパックする
        /// </summary>
        /// <remarks>
        /// 先頭塩基が最上位側、末尾塩基が最下位側に来る (空きビットは下位側に残る)
        /// </remarks>
        internal static ulong Get_パック_小(ReadOnlySpan<byte> p_kmer)
        {
            var l_値 = 0UL;
            foreach (var l_塩基ID in p_kmer)
            {
                l_値 = (l_値 << 2) | ((ulong)l_塩基ID - 1);
            }
            return l_値;
        }

        /// <summary>
        /// Get_パック_小 でパックした値の逆相補を、ヒープ確保なしで直接計算する
        /// </summary>
        /// <remarks>
        /// 2bitコドンごとに相補を取り (A&lt;-&gt;T, C&lt;-&gt;G)、
        /// 下位から順に取り出しつつ上位へ積み直すことでコドン順序も反転させる
        /// </remarks>
        private static ulong Get_逆相補_小(ulong p_パック済み, int p_長さ)
        {
            var l_残り = p_パック済み;
            var l_結果 = 0UL;
            for (var i = 0; i < p_長さ; i++)
            {
                var l_コドン = l_残り & 0x3UL;
                l_結果 = (l_結果 << 2) | (l_コドン ^ 0x3UL);
                l_残り >>= 2;
            }
            return l_結果;
        }

        /// <summary>
        /// k-mer の正規形を返す (k &lt;= 32)
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        /// <returns>正規形</returns>
        private static ulong Get_正規形_小(ReadOnlySpan<byte> p_kmer)
        {
            var l_パック済み = Get_パック_小(p_kmer);
            var l_逆相補 = Get_逆相補_小(l_パック済み, p_kmer.Length);
            return Math.Min(l_パック済み, l_逆相補);
        }

        /// <summary>
        /// ファイル上のパック済みバイト列を、そのままパック値として読み替える (k &lt;= 32)
        /// </summary>
        /// <remarks>
        /// 並びの規約が同じなので、塩基列へ展開して詰め直す必要はない
        /// </remarks>
        internal static ulong Get_読み替え_小(ReadOnlySpan<byte> p_パック済み, int p_余りビット)
        {
            var l_値 = 0UL;
            foreach (var l_バイト in p_パック済み)
            {
                l_値 = (l_値 << 8) | l_バイト;
            }
            return l_値 >> p_余りビット;
        }

        /// <summary>
        /// Get_読み替え_小 の 128 bit 版 (33 &lt;= k &lt;= 64)
        /// </summary>
        internal static UInt128 Get_読み替え_中(ReadOnlySpan<byte> p_パック済み, int p_余りビット)
        {
            UInt128 l_値 = 0;
            foreach (var l_バイト in p_パック済み)
            {
                l_値 = (l_値 << 8) | l_バイト;
            }
            return l_値 >> p_余りビット;
        }

        /// <summary>
        /// パック済みバイト列から塩基ID列を復元する (k &gt; 64 の経路用)
        /// </summary>
        private static byte[] Get_復元_塩基列(ReadOnlySpan<byte> p_パック済み, int p_k長)
        {
            var l_塩基列 = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                var l_バイト = p_パック済み[i / 4];
                var l_ずらし = 6 - (2 * (i % 4));
                l_塩基列[i] = (byte)(((l_バイト >> l_ずらし) & 3) + 1);
            }
            return l_塩基列;
        }

        /// <summary>
        /// Get_パック_小 の逆変換
        /// </summary>
        /// <remarks>
        /// 末尾塩基が最下位ビット側にあるため、末尾から復元する
        /// </remarks>
        private static byte[] Get_復元_小(ulong p_パック済み, int p_長さ)
        {
            var l_塩基列 = new byte[p_長さ];
            for (var i = p_長さ - 1; i >= 0; i--)
            {
                l_塩基列[i] = (byte)((p_パック済み & 0x3UL) + 1);
                p_パック済み >>= 2;
            }
            return l_塩基列;
        }

        /// <summary>
        /// Get_パック_小 の 128 bit 版 (k は 64 以下)
        /// </summary>
        /// <remarks>
        /// ビット配置の規約は同じで、
        /// kmer の先頭塩基が最上位側、末尾塩基が最下位側に来る
        /// </remarks>
        internal static UInt128 Get_パック_中(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_値 = 0;
            foreach (var l_塩基ID in p_kmer)
            {
                l_値 = (l_値 << 2) | (UInt128)(l_塩基ID - 1);
            }
            return l_値;
        }

        /// <summary>
        /// Get_逆相補_小 の 128 bit 版
        /// </summary>
        private static UInt128 Get_逆相補_中(UInt128 p_パック済み, int p_長さ)
        {
            var l_残り = p_パック済み;
            UInt128 l_結果 = 0;
            for (var i = 0; i < p_長さ; i++)
            {
                var l_コドン = l_残り & 3;
                l_結果 = (l_結果 << 2) | (l_コドン ^ 3);
                l_残り >>= 2;
            }
            return l_結果;
        }

        /// <summary>
        /// k-mer の正規形を返す (33 &lt;= k &lt;= 64)
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        /// <returns>正規形</returns>
        private static UInt128 Get_正規形_中(ReadOnlySpan<byte> p_kmer)
        {
            var l_パック済み = Get_パック_中(p_kmer);
            var l_逆相補 = Get_逆相補_中(l_パック済み, p_kmer.Length);
            return l_パック済み < l_逆相補 ? l_パック済み : l_逆相補;
        }

        /// <summary>
        /// Get_パック_中 の逆変換
        /// </summary>
        private static byte[] Get_復元_中(UInt128 p_パック済み, int p_長さ)
        {
            var l_塩基列 = new byte[p_長さ];
            for (var i = p_長さ - 1; i >= 0; i--)
            {
                l_塩基列[i] = (byte)((ulong)(p_パック済み & 3) + 1);
                p_パック済み >>= 2;
            }
            return l_塩基列;
        }

        /// <summary>
        /// カットオフを通過した信頼できる k-merを (正規化された、いずれかの向きの)
        /// byte 配列として 1 件ずつ列挙する
        /// </summary>
        /// <remarks>
        /// GraphSimplifier の tip 除去等、
        /// 集合全体を舐めて再判定する処理から使う
        /// </remarks>
        public IEnumerable<byte[]> Get_信頼kmer一覧()
        {
            var l_k長 = this._k長;
            if (this.小経路を使うか)
            {
                foreach (var l_パック済み in this._信頼kmer_小!.Keys)
                {
                    yield return Get_復元_小(l_パック済み, l_k長);
                }
            }
            else if (this.中経路を使うか)
            {
                foreach (var l_パック済み in this._信頼kmer_中!.Keys)
                {
                    yield return Get_復元_中(l_パック済み, l_k長);
                }
            }
            else
            {
                foreach (var l_キー in this._信頼kmer_大!.Keys)
                {
                    yield return l_キー.Get_塩基列(l_k長);
                }
            }
        }

        /// <summary>
        /// kmer を信頼できる k-mer集合から除去する (順鎖・逆鎖どちらの向きで
        /// 渡してもよい)
        /// </summary>
        /// <remarks>
        /// GraphSimplifier が tip の構成 k-mer を取り除く際に使う
        /// </remarks>
        public void V_除去(ReadOnlySpan<byte> p_kmer)
        {
            _ = this.小経路を使うか
                ? this._信頼kmer_小!.Remove(Get_正規形_小(p_kmer))
                : this.中経路を使うか ? this._信頼kmer_中!.Remove(Get_正規形_中(p_kmer)) : this._信頼kmer_大!.Remove(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// カットオフ後の信頼できる k-mer 集合へ 1 件足す
        /// </summary>
        /// <remarks>
        /// 既にある場合は何もしない<br/>
        /// 前段の k のアセンブリから配列を引き継ぐときに使う<br/>
        /// 戻り値は実際に足したかどうか<br/>
        /// カバレッジは呼び出し側が与える<br/>
        /// ここを名目値で埋めると、
        /// コピー数推定と低カバレッジ端のトリミングが揃って壊れる
        /// </remarks>
        public bool V_追加_信頼kmer(ReadOnlySpan<byte> p_kmer, ulong p_カバレッジ)
        {
            return this.小経路を使うか
                ? this._信頼kmer_小!.TryAdd(Get_正規形_小(p_kmer), p_カバレッジ)
                : this.中経路を使うか
                ? this._信頼kmer_中!.TryAdd(Get_正規形_中(p_kmer), p_カバレッジ)
                : this._信頼kmer_大!.TryAdd(new KmerKey(p_kmer).Get_正規形(), p_カバレッジ);
        }

        /// <summary>
        /// 信頼できる k-mer 集合を走査し、unitig の開始点をすべて再検出する
        /// </summary>
        /// <remarks>
        /// 各座位について両方の向きを個別に判定する<br/>
        /// 開始点判定は向き依存
        /// (その k-mer 自身の入次数を見る) なので、正規形側だけを調べると
        /// 「順鎖では分岐点の直後だが逆鎖ではそうでない」座位を見逃す
        /// </remarks>
        public List<byte[]> Get_開始kmer一覧()
        {
            // 判定は 1 件あたり最大 8 回のハッシュ引きを要し、それを全 k-mer の
            // 両向きについて行う
            // tip 除去は反復のたびにこれを呼ぶため、
            // 単一スレッドだと実行時間の大半をここが占める
            // 判定は読み取りのみなので並列に行える
            //
            // AsOrdered で走査順を保つ
            // unitig 構築の結果は開始点の順序に
            // 依存するため、順序が変わると出力が実行ごとに変わる
            return [.. this.Get_信頼kmer一覧()
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数))
                .SelectMany(this.Get_開始kmer候補)];
        }

        /// <summary>
        /// その座位が開始点になる向きを列挙する (0〜2 件)
        /// </summary>
        /// <remarks>
        /// 開始点判定は向き依存なので、両方の向きを個別に見る必要がある
        /// </remarks>
        private IEnumerable<byte[]> Get_開始kmer候補(byte[] p_kmer)
        {
            if (this.Get_開始kmerか(p_kmer))
            {
                yield return p_kmer;
            }
            var l_逆相補 = Util.V_逆相補(p_kmer).ToArray();
            if (this.Get_開始kmerか(l_逆相補))
            {
                yield return l_逆相補;
            }
        }

        /// <summary>
        /// 全シャードを 1 本のソート済みファイルへ統合し、そのパスを返す
        /// </summary>
        /// <remarks>
        /// 結果は使い回す<br/>
        /// -kc の自動決定がカットオフ前にヒストグラムを
        /// 読むため、やり直すとディスク I/O が丸ごと二重になる
        /// </remarks>
        private string Get_統合済みファイル()
        {
            if (this._統合ファイルパス != null)
            {
                return this._統合ファイルパス;
            }

            var l_統合済みファイル = new List<string>();
            foreach (var l_カウンタ in this._カウンタ群!)
            {
                l_統合済みファイル.Add(l_カウンタ.Get_統合ファイル());
                l_カウンタ.Dispose();
            }
            this._カウンタ群 = null;

            var (l_パス, l_ヒストグラム) =
                CountingDB.Get_統合結果_シャード間(this._一時ディレクトリ, l_統合済みファイル);
            this._統合ファイルパス = l_パス;
            this._統合時のヒストグラム = l_ヒストグラム;
            return this._統合ファイルパス;
        }

        /// <summary>
        /// 出現回数ヒストグラムだけを作る
        /// </summary>
        /// <remarks>
        /// -kc が未指定のとき、
        /// カットオフを決めるために先に呼ぶ<br/>
        /// k-mer の中身は読み捨てる
        /// </remarks>
        public Dictionary<ulong, long> Get_出現回数ヒストグラム()
        {
            var l_ファイルパス = this.Get_統合済みファイル();
            if (this._統合時のヒストグラム is { } l_集計済み)
            {
                return l_集計済み;
            }

            var l_パック長 = (this._k長 + 3) / 4;
            var l_読み捨てバッファ = new byte[l_パック長];

            Dictionary<ulong, long> l_ヒストグラム = [];
            using var l_読み込み = new BinaryReader(File.Open(l_ファイルパス, FileMode.Open, FileAccess.Read));
            while (Util.Get_続きがあるか(l_読み込み))
            {
                _ = l_読み込み.Read(l_読み捨てバッファ, 0, l_パック長);
                var l_出現回数 = l_読み込み.ReadUInt64();
                l_ヒストグラム[l_出現回数] = l_ヒストグラム.GetValueOrDefault(l_出現回数, 0L) + 1;
            }
            return l_ヒストグラム;
        }

        /// <summary>
        /// 出現回数がカットオフに満たない k-mer を落とし、残ったものを信頼できる集合にする
        /// </summary>
        /// <param name="p_カットオフ">残すために必要な出現回数</param>
        /// <returns>walk の開始点になりうる k-mer</returns>
        public List<byte[]> V_カットオフ(ulong p_カットオフ)
        {
            var l_ファイルパス = this.Get_統合済みファイル();

            var l_パック長 = (this._k長 + 3) / 4;
            var l_小経路 = this.小経路を使うか;
            var l_中経路 = this.中経路を使うか;
            Dictionary<KmerKey, ulong>? l_信頼kmer_大;
            Dictionary<ulong, ulong>? l_信頼kmer_小;
            Dictionary<UInt128, ulong>? l_信頼kmer_中;
            {
                ulong l_採用数 = 0UL;
                ulong l_総種類数 = 0UL;
                // 出現回数 -> その回数を持つユニーク k-mer の種類数
                // エラー由来の低頻度 k-mer と真のゲノム由来 k-mer を分ける「谷」を
                // 推定するために、カットオフ判定と同じこのループで集計する
                // (このファイルはこの後削除されるため、ここでしか見られない)
                Dictionary<ulong, long> l_ヒストグラム = [];

                var l_エントリ長 = l_パック長 + sizeof(ulong);
                var l_総エントリ数 = new FileInfo(l_ファイルパス).Length / l_エントリ長;

                // 数千万件を入れるので、伸ばしながらのリハッシュを避けて先に確保する
                var l_見込み = (int)Math.Min(int.MaxValue / 2, Math.Max(1024, l_総エントリ数 / 4));
                l_信頼kmer_小 = l_小経路 ? new Dictionary<ulong, ulong>(l_見込み) : null;
                l_信頼kmer_中 = l_中経路 ? new Dictionary<UInt128, ulong>(l_見込み) : null;
                l_信頼kmer_大 = l_小経路 || l_中経路 ? null : new Dictionary<KmerKey, ulong>(l_見込み);

                // ファイル上の並びは 2 bit/塩基・先頭塩基が最上位で、パック値の規約と
                // 同じ
                // 塩基列へ展開してから詰め直す必要はなく、余りビット分の
                // シフトだけで正規形の計算へ渡せる
                var l_余りビット = (8 * l_パック長) - (2 * this._k長);

                using var l_流れ = new FileStream(
                    l_ファイルパス, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20, FileOptions.SequentialScan);
                var l_バッファ = new byte[l_エントリ長 * 4096];
                var l_残り = 0;
                while (true)
                {
                    var l_読んだ = l_流れ.Read(l_バッファ, l_残り, l_バッファ.Length - l_残り);
                    var l_有効 = l_残り + l_読んだ;
                    var l_位置 = 0;
                    while (l_位置 + l_エントリ長 <= l_有効)
                    {
                        var l_パック済み = l_バッファ.AsSpan(l_位置, l_パック長);
                        var l_出現回数 = BitConverter.ToUInt64(l_バッファ, l_位置 + l_パック長);
                        l_位置 += l_エントリ長;

                        l_総種類数 += 1;
                        l_ヒストグラム[l_出現回数] = l_ヒストグラム.GetValueOrDefault(l_出現回数, 0L) + 1;
                        if (l_出現回数 < p_カットオフ)
                        {
                            continue;
                        }

                        l_採用数 += 1;
                        // カウント段階で既に正規形へ寄せてあるため、同じ正規形が
                        // 複数エントリとして現れることはない
                        // それでも加算で受けて
                        // おけば、将来カウント側の正規化をやめた場合でも壊れない
                        if (l_小経路)
                        {
                            var l_値 = Get_読み替え_小(l_パック済み, l_余りビット);
                            var l_逆 = Get_逆相補_小(l_値, this._k長);
                            var l_正規形 = Math.Min(l_値, l_逆);
                            l_信頼kmer_小![l_正規形] = l_信頼kmer_小.GetValueOrDefault(l_正規形, 0UL) + l_出現回数;
                        }
                        else if (l_中経路)
                        {
                            var l_値 = Get_読み替え_中(l_パック済み, l_余りビット);
                            var l_逆 = Get_逆相補_中(l_値, this._k長);
                            var l_正規形 = l_値 < l_逆 ? l_値 : l_逆;
                            l_信頼kmer_中![l_正規形] = l_信頼kmer_中.GetValueOrDefault(l_正規形, 0UL) + l_出現回数;
                        }
                        else
                        {
                            var l_正規形 = new KmerKey(Get_復元_塩基列(l_パック済み, this._k長)).Get_正規形();
                            l_信頼kmer_大![l_正規形] = l_信頼kmer_大.GetValueOrDefault(l_正規形, 0UL) + l_出現回数;
                        }
                    }

                    if (l_読んだ == 0)
                    {
                        break;
                    }

                    l_残り = l_有効 - l_位置;
                    if (l_残り > 0)
                    {
                        Array.Copy(l_バッファ, l_位置, l_バッファ, 0, l_残り);
                    }
                }

                Logger.V_出力(メッセージID.kmer種類数, l_総種類数);
                Logger.V_出力(メッセージID.採用kmer数, l_採用数);
                this.A_出現回数ヒストグラム = l_ヒストグラム;
            }
            File.Delete(l_ファイルパス);
            this._統合ファイルパス = null;
            this._信頼kmer_大 = l_信頼kmer_大;
            this._信頼kmer_小 = l_信頼kmer_小;
            this._信頼kmer_中 = l_信頼kmer_中;

            Logger.V_出力(メッセージID.開始kmerの探索);
            // 以前はここで一度カットオフ通過 k-mer をファイルへ書き出し、
            // 読み直して開始点を判定していた
            // 厳密な集合をインメモリで
            // 保持するようになったため、その集合を直接走査すれば同じ結果が
            // 得られ、ディスク I/O を 1 往復省略できる
            return this.Get_開始kmer一覧();
        }

        /// <summary>
        /// k-mer が unitig の開始点かどうか
        /// </summary>
        /// <remarks>
        /// 入次数が 1 でも、唯一の予測元が分岐点なら開始点として扱う<br/>
        /// 前進 walk は予測元の時点で停止するため、この k-mer は誰からも
        /// 訪れてもらえず、扱わないと配列が丸ごと欠落する
        /// </remarks>
        public bool Get_開始kmerか(Span<byte> p_kmer)
        {
            var l_入次数 = this.Get_入次数(p_kmer, out var l_唯一の予測元);
            return l_入次数 != 1 || this.Get_出次数(l_唯一の予測元!) != 1;
        }

        /// <summary>
        /// kmer への入次数 (前方に接続しうる異なる 1 塩基拡張の数)
        /// </summary>
        /// <remarks>
        /// 前進伸長が 後続 = kmer[1..] + c である以上、その逆を解くと
        /// 予測元は P = c + kmer[..^1] になる<br/>
        /// kmer[1..] を使うと
        /// c = kmer[0] のとき候補が kmer 自身と一致して常に自己ヒットし、
        /// 入次数 0 が検出できなくなる
        /// </remarks>
        public int Get_入次数(Span<byte> p_kmer)
        {
            return this.Get_入次数(p_kmer, out _);
        }

        /// <summary>
        /// 入次数計算の本体
        /// </summary>
        /// <remarks>
        /// 入次数がちょうど 1 だった場合、その唯一の
        /// 予測元 (kmer 長の byte 配列) も同時に返す (開始点判定が使う)
        /// </remarks>
        private int Get_入次数(Span<byte> p_kmer, out byte[]? p_唯一の予測元)
        {
            var l_候補 = new byte[p_kmer.Length];
            p_kmer[..^1].CopyTo(l_候補.AsSpan(1));
            var l_件数 = 0;
            byte[]? l_一致 = null;
            for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                l_候補[0] = i;
                if (this.Get_含まれるか(l_候補))
                {
                    l_件数++;
                    l_一致 = l_件数 == 1 ? (byte[])l_候補.Clone() : null;
                }
            }
            p_唯一の予測元 = l_件数 == 1 ? l_一致 : null;
            return l_件数;
        }

        /// <summary>
        /// kmer からの出次数 (後方に接続しうる異なる 1 塩基拡張の数) を数える
        /// </summary>
        /// <remarks>
        /// UnitigMaker の前進伸長規則 (kmer[1..] + c) そのものを試す<br/>
        /// GraphSimplifier がunitigの末尾端の次数 (=tip 判定) を見る際に使う
        /// </remarks>
        public int Get_出次数(Span<byte> p_kmer)
        {
            var l_候補 = new byte[p_kmer.Length];
            p_kmer[1..].CopyTo(l_候補.AsSpan(0, p_kmer.Length - 1));
            var l_件数 = 0;
            for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                l_候補[^1] = i;
                if (this.Get_含まれるか(l_候補))
                {
                    l_件数++;
                }
            }
            return l_件数;
        }

        /// <summary>
        /// 保持している資源を解放する
        /// </summary>
        public void Dispose()
        {
            if (this._カウンタ群 != null)
            {
                foreach (var l_カウンタ in this._カウンタ群)
                {
                    l_カウンタ.Dispose();
                }
            }
        }
    }
}
