using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// k-mer の出現回数カウントと、カットオフを通過した信頼できる k-mer の厳密な集合を保持する
    /// </summary>
    internal class TrustedKmerIndex : IDisposable, IKmerLookup
    {
        #region 定数

        /// <summary>
        /// 1 窓で許す曖昧塩基の組み合わせ数
        /// </summary>
        private const int C_曖昧塩基の展開上限 = 64;

        /// <summary>
        /// パック値で扱える k の上限
        /// </summary>
        public const int C_パック値のk上限 = 128;

        /// <summary>
        /// ヒストグラムを配列で数える出現回数の上限
        /// </summary>
        private const int C_配列で数える出現回数の上限 = 1 << 16;

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        private readonly string _一時ディレクトリ;

        /// <summary>
        /// シャードごとの k-mer カウンタ
        /// </summary>
        private CountingDB[]? _カウンタ群;

        /// <summary>
        /// シャードロック
        /// </summary>
        private readonly Lock[]? _シャードロック;

        /// <summary>
        /// 信頼できる k-mer と出現回数 (k &gt; 128)
        /// </summary>
        private Dictionary<KmerKey, ulong>? _信頼kmer_大;

        /// <summary>
        /// 信頼できる k-mer と出現回数 (k &lt;= 32)
        /// </summary>
        private Dictionary<ulong, ulong>? _信頼kmer_小;

        /// <summary>
        /// 信頼できる k-mer と出現回数 (33 &lt;= k &lt;= 64)
        /// </summary>
        private Dictionary<UInt128, ulong>? _信頼kmer_中;

        /// <summary>
        /// 信頼できる k-mer と出現回数 (65 &lt;= k &lt;= 128)
        /// </summary>
        private Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>? _信頼kmer_長;

        /// <summary>
        /// カットオフ未満で控えた k-mer と出現回数 (k &gt; 128)
        /// </summary>
        private Dictionary<KmerKey, ulong>? _控えkmer_大;

        /// <summary>
        /// カットオフ未満で控えた k-mer と出現回数 (k &lt;= 32)
        /// </summary>
        private Dictionary<ulong, ulong>? _控えkmer_小;

        /// <summary>
        /// カットオフ未満で控えた k-mer と出現回数 (33 &lt;= k &lt;= 64)
        /// </summary>
        private Dictionary<UInt128, ulong>? _控えkmer_中;

        /// <summary>
        /// カットオフ未満で控えた k-mer と出現回数 (65 &lt;= k &lt;= 128)
        /// </summary>
        private Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>? _控えkmer_長;

        /// <summary>
        /// シャードごとに統合したファイルのパス
        /// </summary>
        private List<string>? _統合ファイル群;

        /// <summary>
        /// 統合ファイルから集計した、出現回数ごとの k-mer 種類数
        /// </summary>
        private Dictionary<ulong, long>? _統合時のヒストグラム;

        /// <summary>
        /// k 長
        /// </summary>
        private readonly int _k長;

        /// <summary>
        /// 64 bit のキーで所属を調べるか
        /// </summary>
        private readonly bool _Is小経路使用;

        /// <summary>
        /// 128 bit のキーで所属を調べるか
        /// </summary>
        private readonly bool _Is中経路使用;

        /// <summary>
        /// 128 bit 2 語のキーで所属を調べるか
        /// </summary>
        private readonly bool _Is長経路使用;

        #endregion

        #region プロパティ

        /// <summary>
        /// 直近の V_カットオフ で集計した出現回数ヒストグラム
        /// </summary>
        public IReadOnlyDictionary<ulong, long> A_出現回数ヒストグラム { get; private set; } = new Dictionary<ulong, long>();

        /// <summary>
        /// カットオフ未満で控えている k-mer の数
        /// </summary>
        public int A_控えkmer数 => this._控えkmer_小?.Count ?? this._控えkmer_中?.Count ?? this._控えkmer_長?.Count ?? this._控えkmer_大?.Count ?? 0;

        /// <summary>
        /// カウンタのシャード数
        /// </summary>
        public int A_シャード数 => this._シャードロック?.Length ?? 0;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// シャードごとのカウンタを一時ディレクトリの下に用意する
        /// </summary>
        /// <param name="p_一時ディレクトリ">カウンタの置き場</param>
        public TrustedKmerIndex(string p_一時ディレクトリ)
        {
            this._k長 = ConfigurationManager.A_実行時引数.A_k長;
            this._Is小経路使用 = this._k長 <= 32;
            this._Is中経路使用 = this._k長 is > 32 and <= 64;
            this._Is長経路使用 = this._k長 is > 64 and <= C_パック値のk上限;
            this._一時ディレクトリ = p_一時ディレクトリ;
            var l_シャード数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            this._カウンタ群 = new CountingDB[l_シャード数];
            this._シャードロック = new Lock[l_シャード数];
            for (var i = 0; i < l_シャード数; i++)
            {
                this._カウンタ群[i] = new CountingDB(p_一時ディレクトリ, l_シャード数);
                this._シャードロック[i] = new Lock();
            }
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 曖昧塩基を含む k-mer を、ありうる塩基の組み合わせすべてに展開して登録する
        /// </summary>
        /// <param name="p_塩基候補列"></param>
        /// <param name="p_ワーカー番号"></param>
        public void V_登録_曖昧塩基あり(Span<byte[]> p_塩基候補列, int p_ワーカー番号)
        {
            if (this._カウンタ群 is null)
            {
                return;
            }

            var l_組み合わせ数 = 1;
            foreach (var l_候補 in p_塩基候補列)
            {
                if (l_候補.Length == 0 || l_候補.Length > C_曖昧塩基の展開上限 / l_組み合わせ数)
                {
                    return;
                }

                l_組み合わせ数 *= l_候補.Length;
            }

            var l_kmer = new byte[p_塩基候補列.Length];
            this.V_登録_組み合わせ展開(p_塩基候補列, 0, l_kmer, p_ワーカー番号);
        }

        /// <summary>
        /// k-mer を 1 件カウントする
        /// </summary>
        /// <param name="p_kmer"></param>
        public void V_登録(Span<byte> p_kmer)
        {
            if (this._カウンタ群 is not { } l_カウンタ群)
            {
                return;
            }

            if (this._k長 <= C_パック値のk上限)
            {
                var l_値 = Get_正規形_値(p_kmer);
                var l_値のシャード = Get_シャード番号(l_値, l_カウンタ群.Length);
                lock (this._シャードロック![l_値のシャード])
                {
                    l_カウンタ群[l_値のシャード].V_登録_値(l_値);
                }

                return;
            }

            var l_パック済み = TryGet_正規化パック(p_kmer);
            var l_シャード = (int)(Get_ハッシュ(l_パック済み) % (uint)l_カウンタ群.Length);
            lock (this._シャードロック![l_シャード])
            {
                l_カウンタ群[l_シャード].V_登録_パック済み(l_パック済み);
            }
        }

        /// <summary>
        /// 同じシャードへ振り分けた正規形のパック値をまとめて数える (k &lt;= 128)
        /// </summary>
        /// <param name="p_シャード"></param>
        /// <param name="p_値群"></param>
        public void V_登録_値群(int p_シャード, ReadOnlySpan<(UInt128 A_上位, UInt128 A_下位)> p_値群)
        {
            if (this._カウンタ群 is not { } l_カウンタ群)
            {
                return;
            }

            lock (this._シャードロック![p_シャード])
            {
                foreach (var l_値 in p_値群)
                {
                    l_カウンタ群[p_シャード].V_登録_値(l_値);
                }
            }
        }

        /// <summary>
        /// 正規形のパック値を振り分けるシャード
        /// </summary>
        /// <param name="p_値"></param>
        /// <param name="p_シャード数"></param>
        /// <returns></returns>
        public static int Get_シャード番号((UInt128 A_上位, UInt128 A_下位) p_値, int p_シャード数)
        {
            var l_混合 = (ulong)p_値.A_下位 ^ (ulong)(p_値.A_下位 >> 64) ^ (ulong)p_値.A_上位 ^ (ulong)(p_値.A_上位 >> 64);
            l_混合 *= 0x9E37_79B9_7F4A_7C15UL;
            return (int)(((l_混合 >> 32) * (ulong)p_シャード数) >> 32);
        }

        /// <summary>
        /// kmer (順鎖・逆鎖いずれの向きでもよい) がカットオフを通過した信頼できる k-mer 集合に含まれるかどうかを厳密に判定する
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public bool Haskmer(Span<byte> p_kmer)
        {
            return this._Is小経路使用
                ? this._信頼kmer_小!.ContainsKey(Get_正規形_小(p_kmer))
                : this._Is中経路使用
                ? this._信頼kmer_中!.ContainsKey(Get_正規形_中(p_kmer))
                : this._Is長経路使用
                ? this._信頼kmer_長!.ContainsKey(Get_正規形_長(p_kmer))
                : this._信頼kmer_大!.ContainsKey(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// 正規形にパック済みの値で所属を判定する (k &lt;= 32)
        /// </summary>
        /// <param name="p_正規形"></param>
        /// <returns></returns>
        public bool Haskmer_小(ulong p_正規形)
        {
            return this._信頼kmer_小!.ContainsKey(p_正規形);
        }

        /// <summary>
        /// 正規形にパック済みの値で所属を判定する (33 &lt;= k &lt;= 64)
        /// </summary>
        /// <param name="p_正規形"></param>
        /// <returns></returns>
        public bool Haskmer_中(UInt128 p_正規形)
        {
            return this._信頼kmer_中!.ContainsKey(p_正規形);
        }

        /// <summary>
        /// 正規形の右詰めパック値で所属を判定する (k &lt;= 128)
        /// </summary>
        /// <param name="p_上位">128 bit を超える側、k &lt;= 64 なら 0</param>
        /// <param name="p_下位"></param>
        /// <returns></returns>
        public bool Haskmer_正規形(UInt128 p_上位, UInt128 p_下位)
        {
            return this._Is小経路使用
                ? this._信頼kmer_小!.ContainsKey((ulong)p_下位)
                : this._Is中経路使用 ? this._信頼kmer_中!.ContainsKey(p_下位) : this._信頼kmer_長!.ContainsKey((p_上位, p_下位));
        }

        /// <summary>
        /// 正規形のキーで所属を判定する (k &gt; 128)
        /// </summary>
        /// <param name="p_正規形"></param>
        /// <returns></returns>
        public bool Haskmer_正規形(KmerKey p_正規形)
        {
            return this._信頼kmer_大!.ContainsKey(p_正規形);
        }

        /// <summary>
        /// kmer の出現回数 (カバレッジ) を返す
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public ulong Get_カバレッジ(Span<byte> p_kmer)
        {
            return this._Is小経路使用
                ? this._信頼kmer_小!.GetValueOrDefault(Get_正規形_小(p_kmer), 0UL)
                : this._Is中経路使用
                ? this._信頼kmer_中!.GetValueOrDefault(Get_正規形_中(p_kmer), 0UL)
                : this._Is長経路使用
                ? this._信頼kmer_長!.GetValueOrDefault(Get_正規形_長(p_kmer), 0UL)
                : this._信頼kmer_大!.GetValueOrDefault(new KmerKey(p_kmer).Get_正規形(), 0UL);
        }

        /// <summary>
        /// カットオフ未満で控えた k-mer の出現回数を返す
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns>控えに無ければ 0</returns>
        public ulong Get_控えカバレッジ(Span<byte> p_kmer)
        {
            return this._Is小経路使用
                ? this._控えkmer_小?.GetValueOrDefault(Get_正規形_小(p_kmer), 0UL) ?? 0UL
                : this._Is中経路使用
                ? this._控えkmer_中?.GetValueOrDefault(Get_正規形_中(p_kmer), 0UL) ?? 0UL
                : this._Is長経路使用
                ? this._控えkmer_長?.GetValueOrDefault(Get_正規形_長(p_kmer), 0UL) ?? 0UL
                : this._控えkmer_大?.GetValueOrDefault(new KmerKey(p_kmer).Get_正規形(), 0UL) ?? 0UL;
        }

        /// <summary>
        /// 控えた k-mer を手放す
        /// </summary>
        public void V_解放_控え()
        {
            this._控えkmer_小 = null;
            this._控えkmer_中 = null;
            this._控えkmer_長 = null;
            this._控えkmer_大 = null;
        }

        /// <summary>
        /// kmer (塩基 ID 1-4、長さ 32 以下) を 2 bit/塩基で ulong1 個にパックする
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public static ulong TryGet_パック_小(ReadOnlySpan<byte> p_kmer)
        {
            var l_値 = 0UL;
            foreach (var l_塩基ID in p_kmer)
            {
                l_値 = (l_値 << 2) | (l_塩基ID - 1UL);
            }

            return l_値;
        }

        /// <summary>
        /// ファイル上のパック済みバイト列を、そのままパック値として読み替える (k &lt;= 32)
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_余りビット"></param>
        /// <returns></returns>
        public static ulong Get_読み替え_小(ReadOnlySpan<byte> p_パック済み, int p_余りビット)
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
        /// <param name="p_パック済み"></param>
        /// <param name="p_余りビット"></param>
        /// <returns></returns>
        public static UInt128 Get_読み替え_中(ReadOnlySpan<byte> p_パック済み, int p_余りビット)
        {
            UInt128 l_値 = 0;
            foreach (var l_バイト in p_パック済み)
            {
                l_値 = (l_値 << 8) | l_バイト;
            }

            return l_値 >> p_余りビット;
        }

        /// <summary>
        /// Get_読み替え_小 の 128 bit 2 語版 (k &lt;= 128)
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_余りビット"></param>
        /// <returns></returns>
        public static (UInt128 A_上位, UInt128 A_下位) Get_読み替え_長(ReadOnlySpan<byte> p_パック済み, int p_余りビット)
        {
            UInt128 l_上位 = 0;
            UInt128 l_下位 = 0;
            foreach (var l_バイト in p_パック済み)
            {
                l_上位 = (l_上位 << 8) | (l_下位 >> 120);
                l_下位 = (l_下位 << 8) | l_バイト;
            }

            if (p_余りビット > 0)
            {
                l_下位 = (l_下位 >> p_余りビット) | (l_上位 << (128 - p_余りビット));
                l_上位 >>= p_余りビット;
            }

            return (l_上位, l_下位);
        }

        /// <summary>
        /// TryGet_パック_小 の 128 bit 版 (k は 64 以下)
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public static UInt128 TryGet_パック_中(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_値 = 0;
            foreach (var l_塩基ID in p_kmer)
            {
                l_値 = (l_値 << 2) | (UInt128)(l_塩基ID - 1);
            }

            return l_値;
        }

        /// <summary>
        /// TryGet_パック_小 の 128 bit 2 語版 (k は 128 以下)
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public static (UInt128 A_上位, UInt128 A_下位) TryGet_パック_長(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_上位 = 0;
            UInt128 l_下位 = 0;
            foreach (var l_塩基ID in p_kmer)
            {
                l_上位 = (l_上位 << 2) | (l_下位 >> 126);
                l_下位 = (l_下位 << 2) | (UInt128)(l_塩基ID - 1);
            }

            return (l_上位, l_下位);
        }

        /// <summary>
        /// カットオフを通過した信頼できる k-mer を (正規化された、いずれかの向きの) byte 配列として 1 件ずつ列挙する
        /// </summary>
        /// <returns></returns>
        public IEnumerable<byte[]> Get_信頼kmer一覧()
        {
            var l_k長 = this._k長;
            if (this._Is小経路使用)
            {
                foreach (var l_パック済み in this._信頼kmer_小!.Keys)
                {
                    yield return Get_復元_小(l_パック済み, l_k長);
                }
            }
            else if (this._Is中経路使用)
            {
                foreach (var l_パック済み in this._信頼kmer_中!.Keys)
                {
                    yield return Get_復元_中(l_パック済み, l_k長);
                }
            }
            else if (this._Is長経路使用)
            {
                foreach (var l_パック済み in this._信頼kmer_長!.Keys)
                {
                    yield return Get_復元_長(l_パック済み, l_k長);
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
        /// kmer を信頼できる k-mer 集合から除去する (順鎖・逆鎖どちらの向きで渡してもよい)
        /// </summary>
        /// <param name="p_kmer"></param>
        public void V_除去(ReadOnlySpan<byte> p_kmer)
        {
            _ = this._Is小経路使用
                ? this._信頼kmer_小!.Remove(Get_正規形_小(p_kmer))
                : this._Is中経路使用
                ? this._信頼kmer_中!.Remove(Get_正規形_中(p_kmer))
                : this._Is長経路使用 ? this._信頼kmer_長!.Remove(Get_正規形_長(p_kmer)) : this._信頼kmer_大!.Remove(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// カットオフ後の信頼できる k-mer 集合へ 1 件足す
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <param name="p_カバレッジ"></param>
        /// <returns></returns>
        public bool Try追加_信頼kmer(ReadOnlySpan<byte> p_kmer, ulong p_カバレッジ)
        {
            return this._Is小経路使用
                ? this._信頼kmer_小!.TryAdd(Get_正規形_小(p_kmer), p_カバレッジ)
                : this._Is中経路使用
                ? this._信頼kmer_中!.TryAdd(Get_正規形_中(p_kmer), p_カバレッジ)
                : this._Is長経路使用
                ? this._信頼kmer_長!.TryAdd(Get_正規形_長(p_kmer), p_カバレッジ)
                : this._信頼kmer_大!.TryAdd(new KmerKey(p_kmer).Get_正規形(), p_カバレッジ);
        }

        /// <summary>
        /// 正規形の右詰めパック値で信頼できる k-mer 集合へ 1 件足す (k &lt;= 128)
        /// </summary>
        /// <param name="p_上位">128 bit を超える側、k &lt;= 64 なら 0</param>
        /// <param name="p_下位"></param>
        /// <param name="p_カバレッジ">未登録の場合に設定する観測回数</param>
        /// <returns>未登録のキーを追加した場合は true</returns>
        public bool Try追加_信頼kmer_正規形(UInt128 p_上位, UInt128 p_下位, ulong p_カバレッジ)
        {
            return this._Is小経路使用
                ? this._信頼kmer_小!.TryAdd((ulong)p_下位, p_カバレッジ)
                : this._Is中経路使用 ? this._信頼kmer_中!.TryAdd(p_下位, p_カバレッジ) : this._信頼kmer_長!.TryAdd((p_上位, p_下位), p_カバレッジ);
        }

        /// <summary>
        /// 正規形のキーで信頼できる k-mer 集合へ 1 件足す (k &gt; 128)
        /// </summary>
        /// <param name="p_正規形">作業領域から切り離したキー</param>
        /// <param name="p_カバレッジ">未登録の場合に設定する観測回数</param>
        /// <returns>未登録のキーを追加した場合は true</returns>
        public bool Try追加_信頼kmer_正規形(KmerKey p_正規形, ulong p_カバレッジ)
        {
            return this._信頼kmer_大!.TryAdd(p_正規形, p_カバレッジ);
        }

        /// <summary>
        /// 信頼できる k-mer 集合を走査し、unitig の開始点をすべて再検出する
        /// </summary>
        /// <returns></returns>
        public List<byte[]> Get_開始kmer一覧()
        {
            return [.. this.Get_信頼kmer一覧()
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数))
                .SelectMany(this.Get_開始kmer候補)];
        }

        /// <summary>
        /// 出現回数ヒストグラムだけを作る
        /// </summary>
        /// <returns></returns>
        public Dictionary<ulong, long> Get_出現回数ヒストグラム()
        {
            var l_ファイル群 = this.Get_統合ファイル群();
            if (this._統合時のヒストグラム is { } l_集計済み)
            {
                return l_集計済み;
            }

            using var l_計測 = new StageTimer($"kmer-histogram k={this._k長}");
            var l_パック長 = (this._k長 + 3) / 4;
            var l_シャード別 = new (long[] A_配列, Dictionary<ulong, long> A_大きい回数)[l_ファイル群.Count];
            _ = Parallel.For(0, l_ファイル群.Count, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数) }, s =>
            {
                var l_配列 = new long[C_配列で数える出現回数の上限];
                Dictionary<ulong, long> l_大きい回数 = [];
                V_走査_エントリ(l_ファイル群[s], l_パック長, (_, l_出現回数) => V_加算_ヒストグラム(l_配列, l_大きい回数, l_出現回数));
                l_シャード別[s] = (l_配列, l_大きい回数);
            });

            this._統合時のヒストグラム = Get_合算ヒストグラム(l_シャード別);
            return this._統合時のヒストグラム;
        }

        /// <summary>
        /// 出現回数がカットオフに満たない k-mer を落とし、残ったものを信頼できる集合にして開始点を返す
        /// </summary>
        /// <param name="p_カットオフ">残すために必要な出現回数</param>
        /// <returns>walk の開始点になりうる k-mer</returns>
        public List<byte[]> V_カットオフ(ulong p_カットオフ)
        {
            this.V_適用_カットオフ(p_カットオフ);
            Logger.V_出力(メッセージID.開始kmerの探索);
            return this.Get_開始kmer一覧();
        }

        /// <summary>
        /// 出現回数がカットオフに満たない k-mer を落とし、残ったものを信頼できる集合にする
        /// </summary>
        /// <param name="p_カットオフ">残すために必要な出現回数</param>
        /// <param name="p_控え下限">カットオフ未満でも控えておく出現回数の下限、0 なら控えない</param>
        public void V_適用_カットオフ(ulong p_カットオフ, ulong p_控え下限 = 0UL)
        {
            var l_ファイル群 = this.Get_統合ファイル群();
            var l_パック長 = (this._k長 + 3) / 4;
            var l_余りビット = (8 * l_パック長) - (2 * this._k長);
            var l_Is控え使用 = p_控え下限 > 0UL && p_控え下限 < p_カットオフ;
            var l_k長 = this._k長;
            var l_Isパック値 = l_k長 <= C_パック値のk上限;

            var l_シャード別ヒストグラム = new (long[] A_配列, Dictionary<ulong, long> A_大きい回数)[l_ファイル群.Count];
            var l_シャード別採用 = new List<(UInt128 A_上位, UInt128 A_下位, ulong A_出現回数)>[l_ファイル群.Count];
            var l_シャード別控え = new List<(UInt128 A_上位, UInt128 A_下位, ulong A_出現回数)>[l_ファイル群.Count];
            var l_シャード別採用_大 = new List<(KmerKey A_キー, ulong A_出現回数)>[l_ファイル群.Count];
            var l_シャード別控え_大 = new List<(KmerKey A_キー, ulong A_出現回数)>[l_ファイル群.Count];
            var l_シャード別種類数 = new long[l_ファイル群.Count];

            _ = Parallel.For(0, l_ファイル群.Count, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数) }, s =>
            {
                var l_配列 = new long[C_配列で数える出現回数の上限];
                Dictionary<ulong, long> l_大きい回数 = [];
                List<(UInt128 A_上位, UInt128 A_下位, ulong A_出現回数)> l_採用 = [];
                List<(UInt128 A_上位, UInt128 A_下位, ulong A_出現回数)> l_控え = [];
                List<(KmerKey A_キー, ulong A_出現回数)> l_採用_大 = [];
                List<(KmerKey A_キー, ulong A_出現回数)> l_控え_大 = [];
                var l_種類数 = 0L;

                V_走査_エントリ(l_ファイル群[s], l_パック長, (l_パック済み, l_出現回数) =>
                {
                    l_種類数++;
                    V_加算_ヒストグラム(l_配列, l_大きい回数, l_出現回数);
                    var l_Is採用 = l_出現回数 >= p_カットオフ;
                    if (!l_Is採用 && (!l_Is控え使用 || l_出現回数 < p_控え下限))
                    {
                        return;
                    }

                    if (l_Isパック値)
                    {
                        var (l_上位, l_下位) = Get_正規形_読み替え(l_パック済み, l_余りビット, l_k長);
                        (l_Is採用 ? l_採用 : l_控え).Add((l_上位, l_下位, l_出現回数));
                    }
                    else
                    {
                        var l_キー = new KmerKey(Get_復元_塩基列(l_パック済み, l_k長)).Get_正規形();
                        (l_Is採用 ? l_採用_大 : l_控え_大).Add((l_キー, l_出現回数));
                    }
                });

                l_シャード別ヒストグラム[s] = (l_配列, l_大きい回数);
                l_シャード別採用[s] = l_採用;
                l_シャード別控え[s] = l_控え;
                l_シャード別採用_大[s] = l_採用_大;
                l_シャード別控え_大[s] = l_控え_大;
                l_シャード別種類数[s] = l_種類数;
            });

            var l_採用数 = l_Isパック値 ? l_シャード別採用.Sum(x => (long)x.Count) : l_シャード別採用_大.Sum(x => (long)x.Count);
            var l_控え数 = l_Isパック値 ? l_シャード別控え.Sum(x => (long)x.Count) : l_シャード別控え_大.Sum(x => (long)x.Count);
            var l_採用容量 = (int)Math.Min(int.MaxValue / 2, Math.Max(1_024L, l_採用数));
            var l_控え容量 = (int)Math.Min(int.MaxValue / 2, l_控え数);

            this._信頼kmer_小 = this._Is小経路使用 ? new Dictionary<ulong, ulong>(l_採用容量) : null;
            this._信頼kmer_中 = this._Is中経路使用 ? new Dictionary<UInt128, ulong>(l_採用容量) : null;
            this._信頼kmer_長 = this._Is長経路使用 ? new Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>(l_採用容量) : null;
            this._信頼kmer_大 = l_Isパック値 ? null : new Dictionary<KmerKey, ulong>(l_採用容量);
            this._控えkmer_小 = l_Is控え使用 && this._Is小経路使用 ? new Dictionary<ulong, ulong>(l_控え容量) : null;
            this._控えkmer_中 = l_Is控え使用 && this._Is中経路使用 ? new Dictionary<UInt128, ulong>(l_控え容量) : null;
            this._控えkmer_長 = l_Is控え使用 && this._Is長経路使用 ? new Dictionary<(UInt128 A_上位, UInt128 A_下位), ulong>(l_控え容量) : null;
            this._控えkmer_大 = l_Is控え使用 && !l_Isパック値 ? new Dictionary<KmerKey, ulong>(l_控え容量) : null;

            for (var s = 0; s < l_ファイル群.Count; s++)
            {
                foreach (var (l_上位, l_下位, l_出現回数) in l_シャード別採用[s])
                {
                    this.V_加算_正規形(l_上位, l_下位, l_出現回数, p_Is控え: false);
                }

                foreach (var (l_上位, l_下位, l_出現回数) in l_シャード別控え[s])
                {
                    this.V_加算_正規形(l_上位, l_下位, l_出現回数, p_Is控え: true);
                }

                foreach (var (l_キー, l_出現回数) in l_シャード別採用_大[s])
                {
                    this._信頼kmer_大![l_キー] = this._信頼kmer_大.GetValueOrDefault(l_キー, 0UL) + l_出現回数;
                }

                foreach (var (l_キー, l_出現回数) in l_シャード別控え_大[s])
                {
                    this._控えkmer_大![l_キー] = this._控えkmer_大.GetValueOrDefault(l_キー, 0UL) + l_出現回数;
                }
            }

            Logger.V_出力(メッセージID.kmer種類数, (ulong)l_シャード別種類数.Sum());
            Logger.V_出力(メッセージID.採用kmer数, (ulong)l_採用数);
            this.A_出現回数ヒストグラム = Get_合算ヒストグラム(l_シャード別ヒストグラム);

            foreach (var l_ファイル in l_ファイル群)
            {
                中間データ置き場.V_削除(l_ファイル);
            }

            this._統合ファイル群 = null;
        }

        /// <summary>
        /// k-mer が unitig の開始点かどうか
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public bool Is開始kmer(Span<byte> p_kmer)
        {
            var l_入次数 = this.Get_入次数(p_kmer, out var l_唯一の予測元);
            return l_入次数 != 1 || this.Get_出次数(l_唯一の予測元!) != 1;
        }

        /// <summary>
        /// kmer への入次数 (前方に接続しうる異なる 1 塩基拡張の数)
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public int Get_入次数(Span<byte> p_kmer)
        {
            return this.Get_入次数(p_kmer, out _);
        }

        /// <summary>
        /// kmer からの出次数 (後方に接続しうる異なる 1 塩基拡張の数) を数える
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public int Get_出次数(Span<byte> p_kmer)
        {
            var l_候補 = p_kmer.Length <= 256 ? stackalloc byte[p_kmer.Length] : new byte[p_kmer.Length];
            p_kmer[1..].CopyTo(l_候補);
            var l_件数 = 0;
            for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                l_候補[^1] = i;
                if (this.Haskmer(l_候補))
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

            if (this._統合ファイル群 != null)
            {
                foreach (var l_ファイル in this._統合ファイル群)
                {
                    中間データ置き場.V_削除(l_ファイル);
                }
            }
        }

        /// <summary>
        /// k-mer の正規形を返す (k &lt;= 32)
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        /// <returns>正規形</returns>
        internal static ulong Get_正規形_小(ReadOnlySpan<byte> p_kmer)
        {
            var l_パック済み = TryGet_パック_小(p_kmer);
            var l_逆相補 = Get_逆相補_小(l_パック済み, p_kmer.Length);
            return Math.Min(l_パック済み, l_逆相補);
        }

        /// <summary>
        /// k-mer の正規形を返す (33 &lt;= k &lt;= 64)
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        /// <returns>正規形</returns>
        internal static UInt128 Get_正規形_中(ReadOnlySpan<byte> p_kmer)
        {
            var l_パック済み = TryGet_パック_中(p_kmer);
            var l_逆相補 = Get_逆相補_中(l_パック済み, p_kmer.Length);
            return l_パック済み < l_逆相補 ? l_パック済み : l_逆相補;
        }

        /// <summary>
        /// k-mer の正規形を返す (k &lt;= 128)
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        /// <returns>正規形の右詰めパック値</returns>
        internal static (UInt128 A_上位, UInt128 A_下位) Get_正規形_長(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_順上 = 0;
            UInt128 l_順下 = 0;
            UInt128 l_逆上 = 0;
            UInt128 l_逆下 = 0;
            var l_末尾 = p_kmer.Length - 1;
            for (var i = 0; i <= l_末尾; i++)
            {
                l_順上 = (l_順上 << 2) | (l_順下 >> 126);
                l_順下 = (l_順下 << 2) | (UInt128)(p_kmer[i] - 1);
                l_逆上 = (l_逆上 << 2) | (l_逆下 >> 126);
                l_逆下 = (l_逆下 << 2) | (UInt128)(4 - p_kmer[l_末尾 - i]);
            }

            return l_順上 < l_逆上 || (l_順上 == l_逆上 && l_順下 <= l_逆下) ? (l_順上, l_順下) : (l_逆上, l_逆下);
        }

        /// <summary>
        /// k-mer の正規形を、k によらず右詰めのパック値で返す (k &lt;= 128)
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        /// <returns>正規形の右詰めパック値、k &lt;= 64 なら上位は 0</returns>
        internal static (UInt128 A_上位, UInt128 A_下位) Get_正規形_値(ReadOnlySpan<byte> p_kmer)
        {
            return p_kmer.Length <= 64 ? ((UInt128)0, Get_正規形_中(p_kmer)) : Get_正規形_長(p_kmer);
        }

        /// <summary>
        /// TryGet_パック_小 でパックした値の逆相補を、ヒープ確保なしで直接計算する
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_長さ"></param>
        /// <returns></returns>
        internal static ulong Get_逆相補_小(ulong p_パック済み, int p_長さ)
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
        /// Get_逆相補_小 の 128 bit 版
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_長さ"></param>
        /// <returns></returns>
        internal static UInt128 Get_逆相補_中(UInt128 p_パック済み, int p_長さ)
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

        #endregion

        #region 内部メソッド

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
                this.V_登録(p_kmer.AsSpan());
                return;
            }

            foreach (var l_塩基ID in p_塩基候補列[p_位置])
            {
                p_kmer[p_位置] = l_塩基ID;
                this.V_登録_組み合わせ展開(p_塩基候補列, p_位置 + 1, p_kmer, p_ワーカー番号);
            }
        }

        /// <summary>
        /// 右詰めのパック値を集合へ加算する
        /// </summary>
        /// <param name="p_上位"></param>
        /// <param name="p_下位"></param>
        /// <param name="p_出現回数"></param>
        /// <param name="p_Is控え">控えの集合へ入れるか</param>
        private void V_加算_正規形(UInt128 p_上位, UInt128 p_下位, ulong p_出現回数, bool p_Is控え)
        {
            if (this._Is小経路使用)
            {
                var l_集合 = p_Is控え ? this._控えkmer_小! : this._信頼kmer_小!;
                l_集合[(ulong)p_下位] = l_集合.GetValueOrDefault((ulong)p_下位, 0UL) + p_出現回数;
            }
            else if (this._Is中経路使用)
            {
                var l_集合 = p_Is控え ? this._控えkmer_中! : this._信頼kmer_中!;
                l_集合[p_下位] = l_集合.GetValueOrDefault(p_下位, 0UL) + p_出現回数;
            }
            else
            {
                var l_集合 = p_Is控え ? this._控えkmer_長! : this._信頼kmer_長!;
                l_集合[(p_上位, p_下位)] = l_集合.GetValueOrDefault((p_上位, p_下位), 0UL) + p_出現回数;
            }
        }

        /// <summary>
        /// ファイル上のパック済みバイト列を読み替え、逆相補と比べた正規形にする (k &lt;= 128)
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_余りビット"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        private static (UInt128 A_上位, UInt128 A_下位) Get_正規形_読み替え(ReadOnlySpan<byte> p_パック済み, int p_余りビット, int p_k長)
        {
            if (p_k長 <= 64)
            {
                var l_値 = Get_読み替え_中(p_パック済み, p_余りビット);
                var l_逆 = Get_逆相補_中(l_値, p_k長);
                return (0, l_値 < l_逆 ? l_値 : l_逆);
            }

            return Get_正規形_長(Get_復元_長(Get_読み替え_長(p_パック済み, p_余りビット), p_k長));
        }

        /// <summary>
        /// 統合ファイルのエントリを先頭から順に渡す
        /// </summary>
        /// <param name="p_パス"></param>
        /// <param name="p_パック長"></param>
        /// <param name="p_処理">パック済みキーと出現回数を受け取る</param>
        private static void V_走査_エントリ(string p_パス, int p_パック長, Action<ReadOnlySpan<byte>, ulong> p_処理)
        {
            var l_エントリ長 = p_パック長 + sizeof(ulong);
            using var l_流れ = 中間データ置き場.Get_読込ストリーム(p_パス);
            var l_バッファ = new byte[l_エントリ長 * 4_096];
            var l_残り = 0;
            while (true)
            {
                var l_読んだ = l_流れ.Read(l_バッファ, l_残り, l_バッファ.Length - l_残り);
                var l_有効 = l_残り + l_読んだ;
                var l_位置 = 0;
                while (l_位置 + l_エントリ長 <= l_有効)
                {
                    p_処理(l_バッファ.AsSpan(l_位置, p_パック長), BitConverter.ToUInt64(l_バッファ, l_位置 + p_パック長));
                    l_位置 += l_エントリ長;
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
        }

        /// <summary>
        /// 出現回数を 1 件ヒストグラムへ足す
        /// </summary>
        /// <param name="p_配列">小さい出現回数を数える配列</param>
        /// <param name="p_大きい回数">配列に収まらない出現回数</param>
        /// <param name="p_出現回数"></param>
        private static void V_加算_ヒストグラム(long[] p_配列, Dictionary<ulong, long> p_大きい回数, ulong p_出現回数)
        {
            if (p_出現回数 < (ulong)p_配列.Length)
            {
                p_配列[p_出現回数]++;
                return;
            }

            p_大きい回数[p_出現回数] = p_大きい回数.GetValueOrDefault(p_出現回数, 0L) + 1L;
        }

        /// <summary>
        /// シャードごとのヒストグラムを 1 つにまとめる
        /// </summary>
        /// <param name="p_シャード別"></param>
        /// <returns></returns>
        private static Dictionary<ulong, long> Get_合算ヒストグラム((long[] A_配列, Dictionary<ulong, long> A_大きい回数)[] p_シャード別)
        {
            Dictionary<ulong, long> l_結果 = [];
            foreach (var (l_配列, l_大きい回数) in p_シャード別)
            {
                for (var i = 0; i < l_配列.Length; i++)
                {
                    if (l_配列[i] > 0L)
                    {
                        l_結果[(ulong)i] = l_結果.GetValueOrDefault((ulong)i, 0L) + l_配列[i];
                    }
                }

                foreach (var (l_出現回数, l_種類数) in l_大きい回数)
                {
                    l_結果[l_出現回数] = l_結果.GetValueOrDefault(l_出現回数, 0L) + l_種類数;
                }
            }

            return l_結果;
        }

        /// <summary>
        /// k-mer を正規形の向きで 2 bit パックする
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        private static byte[] TryGet_正規化パック(ReadOnlySpan<byte> p_kmer)
        {
            var l_Is順鎖使用 = Is順鎖正規形(p_kmer);
            var l_パック済み = new byte[(p_kmer.Length + 3) / 4];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                var l_塩基ID = l_Is順鎖使用 ? p_kmer[i] : (byte)(5 - p_kmer[p_kmer.Length - 1 - i]);
                l_パック済み[i >> 2] |= (byte)((l_塩基ID - 1) << ((3 - (i & 3)) << 1));
            }

            return l_パック済み;
        }

        /// <summary>
        /// 順鎖側がその逆相補以下 (辞書順) かどうか
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        private static bool Is順鎖正規形(ReadOnlySpan<byte> p_kmer)
        {
            int l_i = 0, l_j = p_kmer.Length - 1;
            while (l_i <= l_j)
            {
                var l_順鎖 = p_kmer[l_i];
                var l_逆鎖 = (byte)(5 - p_kmer[l_j]);
                if (l_順鎖 != l_逆鎖)
                {
                    return l_順鎖 < l_逆鎖;
                }

                l_i++;
                l_j--;
            }

            return true;
        }

        /// <summary>
        /// パック済みキーの FNV-1 a ハッシュ
        /// </summary>
        /// <param name="p_パック済みkmer"></param>
        /// <returns></returns>
        private static uint Get_ハッシュ(byte[] p_パック済みkmer)
        {
            var l_ハッシュ = 2_166_136_261U;
            foreach (var l_バイト in p_パック済みkmer)
            {
                l_ハッシュ ^= l_バイト;
                l_ハッシュ *= 16_777_619U;
            }

            return l_ハッシュ;
        }

        /// <summary>
        /// パック済みバイト列から塩基 ID 列を復元する (k &gt; 128 の経路用)
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
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
        /// TryGet_パック_小 の逆変換
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_長さ"></param>
        /// <returns></returns>
        private static byte[] Get_復元_小(ulong p_パック済み, int p_長さ)
        {
            var l_塩基列 = new byte[p_長さ];
            for (var i = p_長さ - 1; i >= 0; i--)
            {
                l_塩基列[i] = (byte)((p_パック済み & 0x3UL) + 1UL);
                p_パック済み >>= 2;
            }

            return l_塩基列;
        }

        /// <summary>
        /// TryGet_パック_中 の逆変換
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_長さ"></param>
        /// <returns></returns>
        private static byte[] Get_復元_中(UInt128 p_パック済み, int p_長さ)
        {
            var l_塩基列 = new byte[p_長さ];
            for (var i = p_長さ - 1; i >= 0; i--)
            {
                l_塩基列[i] = (byte)((ulong)(p_パック済み & 3) + 1UL);
                p_パック済み >>= 2;
            }

            return l_塩基列;
        }

        /// <summary>
        /// TryGet_パック_長 の逆変換
        /// </summary>
        /// <param name="p_パック済み"></param>
        /// <param name="p_長さ"></param>
        /// <returns></returns>
        private static byte[] Get_復元_長((UInt128 A_上位, UInt128 A_下位) p_パック済み, int p_長さ)
        {
            var (l_上位, l_下位) = p_パック済み;
            var l_塩基列 = new byte[p_長さ];
            for (var i = p_長さ - 1; i >= 0; i--)
            {
                l_塩基列[i] = (byte)((ulong)(l_下位 & 3) + 1UL);
                l_下位 = (l_下位 >> 2) | (l_上位 << 126);
                l_上位 >>= 2;
            }

            return l_塩基列;
        }

        /// <summary>
        /// その座位が開始点になる向きを列挙する (0〜2 件)
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        private IEnumerable<byte[]> Get_開始kmer候補(byte[] p_kmer)
        {
            if (this.Is開始kmer(p_kmer))
            {
                yield return p_kmer;
            }

            var l_逆相補 = Util.V_逆相補(p_kmer).ToArray();
            if (this.Is開始kmer(l_逆相補))
            {
                yield return l_逆相補;
            }
        }

        /// <summary>
        /// 全シャードをそれぞれ 1 本のソート済みファイルへ統合し、そのパスを返す
        /// </summary>
        /// <returns></returns>
        private List<string> Get_統合ファイル群()
        {
            if (this._統合ファイル群 != null)
            {
                return this._統合ファイル群;
            }

            using var l_計測 = new StageTimer($"kmer-merge k={this._k長}");
            var l_カウンタ群 = this._カウンタ群!;
            var l_ファイル群 = new string[l_カウンタ群.Length];
            _ = Parallel.For(0, l_カウンタ群.Length, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数) }, s =>
            {
                l_ファイル群[s] = l_カウンタ群[s].Get_統合ファイル();
                l_カウンタ群[s].Dispose();
            });
            this._カウンタ群 = null;
            this._統合ファイル群 = [.. l_ファイル群];
            return this._統合ファイル群;
        }

        /// <summary>
        /// 入次数計算の本体
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <param name="p_唯一の予測元">入次数がちょうど 1 だった場合の、その唯一の予測元</param>
        /// <returns></returns>
        private int Get_入次数(Span<byte> p_kmer, out byte[]? p_唯一の予測元)
        {
            var l_候補 = p_kmer.Length <= 256 ? stackalloc byte[p_kmer.Length] : new byte[p_kmer.Length];
            p_kmer[..^1].CopyTo(l_候補[1..]);
            var l_件数 = 0;
            var l_一致した塩基 = (byte)0;
            for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                l_候補[0] = i;
                if (this.Haskmer(l_候補))
                {
                    l_件数++;
                    l_一致した塩基 = i;
                }
            }

            if (l_件数 == 1)
            {
                l_候補[0] = l_一致した塩基;
                p_唯一の予測元 = l_候補.ToArray();
            }
            else
            {
                p_唯一の予測元 = null;
            }

            return l_件数;
        }

        #endregion
    }
}
