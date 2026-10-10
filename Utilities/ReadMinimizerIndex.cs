using Tsumiki.Commons;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// リードを minimizer で引けるようにした索引 (長い配列がリードかその逆相補に出てくるかを答える)
    /// </summary>
    internal sealed class ReadMinimizerIndex
    {
        #region 定数

        /// <summary>
        /// minimizer にする k-mer の長さ
        /// </summary>
        public const int C_種長 = 15;


        /// <summary>
        /// minimizer を選ぶ窓に並ぶ k-mer の数
        /// </summary>
        public const int C_窓の種数 = 27;

        /// <summary>
        /// 問い合わせられる配列の最短の長さ (これより短いと、どの minimizer も丸ごと含むとは限らない)
        /// </summary>
        public const int C_最短の問い合わせ長 = C_種長 + C_窓の種数 - 1;

        /// <summary>
        /// 短い問い合わせ用の索引の窓の種数 (既定より小さく、最短の問い合わせ長を短くする)
        /// </summary>
        public const int C_短い問い合わせ用の窓の種数 = 17;

        /// <summary>
        /// 短い問い合わせ用の索引に問い合わせられる配列の最短の長さ
        /// </summary>
        public const int C_短い問い合わせの最短長 = C_種長 + C_短い問い合わせ用の窓の種数 - 1;

        /// <summary>
        /// 1 語に詰める文字数
        /// </summary>
        private const int C_語あたりの文字数 = 32;

        /// <summary>
        /// 2 bit の値の順の塩基
        /// </summary>
        private const string C_塩基の並び = "ACGT";

        /// <summary>
        /// minimizer を並列に集めるときの 1 束の区間数
        /// </summary>
        private const int C_区間の束の大きさ = 4_096;

        /// <summary>
        /// 一度に詰める区間の数
        /// </summary>
        private const int C_詰める束の大きさ = 8_192;

        /// <summary>
        /// 詰める先の最初の語数 (足りなければ倍に広げる)
        /// </summary>
        private const int C_語の初期数 = 1 << 20;

        /// <summary>
        /// 並べるときに振り分ける種の上位ビットの数
        /// </summary>
        private const int C_振り分けのビット数 = 10;

        /// <summary>
        /// 種の値のマスク
        /// </summary>
        private const ulong C_種のマスク = (1UL << (2 * C_種長)) - 1;

        #endregion

        #region 内部変数

        /// <summary>
        /// 2 bit で詰めたリード (A・C・G・T が続く区間ごと)
        /// </summary>
        private readonly ulong[] _語;

        /// <summary>
        /// 区間の最後の文字の位置の印
        /// </summary>
        private readonly ulong[] _末尾印;

        /// <summary>
        /// minimizer の正準値 (昇順)
        /// </summary>
        private readonly uint[] _種;

        /// <summary>
        /// minimizer の位置 (_種 と同じ並び、塩基数が 32 bit に収まるとき)
        /// </summary>
        private readonly uint[]? _位置32;

        /// <summary>
        /// minimizer の位置 (_種 と同じ並び、塩基数が 32 bit に収まらないとき)
        /// </summary>
        private readonly long[]? _位置64;

        /// <summary>
        /// 問い合わせに使える長さの区間 (開始位置と長さ、位置の順、窓の種数によらず共通)
        /// </summary>
        private readonly (long A_開始, int A_長さ)[] _区間群;

        #endregion

        #region プロパティ

        /// <summary>
        /// 索引に入れた塩基数
        /// </summary>
        public long A_塩基数 { get; }

        /// <summary>
        /// minimizer を選ぶ窓に並ぶ k-mer の数
        /// </summary>
        public int A_窓の種数 { get; }

        /// <summary>
        /// 問い合わせられる配列の最短の長さ (これより短いと、どの minimizer も丸ごと含むとは限らない)
        /// </summary>
        public int A_最短の問い合わせ長 => C_種長 + this.A_窓の種数 - 1;

        /// <summary>
        /// 索引が使うおおよそのバイト数
        /// </summary>
        public long A_使用量 => (8L * (this._語.LongLength + this._末尾印.LongLength + (this._位置64?.LongLength ?? 0))) + (4L * (this._種.LongLength + (this._位置32?.LongLength ?? 0)));

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 組み立てた中身で作る
        /// </summary>
        /// <param name="p_語">詰めたリード</param>
        /// <param name="p_末尾印">区間の末尾の印</param>
        /// <param name="p_種">minimizer の正準値 (昇順)</param>
        /// <param name="p_位置32">minimizer の位置 (32 bit)</param>
        /// <param name="p_位置64">minimizer の位置 (64 bit)</param>
        /// <param name="p_塩基数">塩基数</param>
        /// <param name="p_窓の種数">minimizer を選ぶ窓に並ぶ k-mer の数</param>
        /// <param name="p_区間群">問い合わせに使える長さの区間</param>
        private ReadMinimizerIndex(ulong[] p_語, ulong[] p_末尾印, uint[] p_種, uint[]? p_位置32, long[]? p_位置64, long p_塩基数, int p_窓の種数, (long A_開始, int A_長さ)[] p_区間群)
        {
            this._語 = p_語;
            this._末尾印 = p_末尾印;
            this._種 = p_種;
            this._位置32 = p_位置32;
            this._位置64 = p_位置64;
            this.A_塩基数 = p_塩基数;
            this.A_窓の種数 = p_窓の種数;
            this._区間群 = p_区間群;
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群から索引を作る
        /// </summary>
        /// <param name="p_配列列">配列を先頭から返すもの (2 回呼ぶ)</param>
        /// <param name="p_窓の種数">minimizer を選ぶ窓に並ぶ k-mer の数 (17 以上、既定は 27)</param>
        /// <returns>作った索引</returns>
        public static ReadMinimizerIndex V_構築(Func<IEnumerable<string>> p_配列列, int p_窓の種数 = C_窓の種数)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(p_窓の種数, C_短い問い合わせ用の窓の種数);
            var l_語 = new ulong[C_語の初期数];
            var l_末尾印 = new ulong[C_語の初期数];
            List<(long A_開始, int A_長さ)> l_区間群 = [];
            List<string> l_束 = new(C_詰める束の大きさ);
            long l_位置 = 0;
            foreach (var l_配列 in p_配列列())
            {
                l_束.Add(l_配列);
                if (l_束.Count >= C_詰める束の大きさ)
                {
                    l_位置 = V_詰める_束(ref l_語, ref l_末尾印, l_区間群, l_束, l_位置);
                    l_束.Clear();
                }
            }

            var l_長さ = V_詰める_束(ref l_語, ref l_末尾印, l_区間群, l_束, l_位置);
            Array.Resize(ref l_語, (int)((l_長さ / C_語あたりの文字数) + 2));
            Array.Resize(ref l_末尾印, (int)((l_長さ / 64) + 2));

            return V_組み立てる(l_語, l_末尾印, l_区間群.ToArray(), l_長さ, p_窓の種数);
        }

        /// <summary>
        /// 配列がリードかその逆相補に 1 回以上出てくるか
        /// </summary>
        /// <param name="p_配列">調べる配列 (最短の問い合わせ長以上、A・C・G・T だけ)</param>
        /// <returns>出てくれば true</returns>
        public bool Has出現(ReadOnlySpan<char> p_配列)
        {
            return this.Get_出現数(p_配列, 1) > 0;
        }

        /// <summary>
        /// 長さ p_窓長 の全ての窓について、Has出現 と同じ答えを一度に求める (配列全体の 15-mer を一度だけ計算する)
        /// </summary>
        /// <param name="p_配列">調べる配列 (A・C・G・T 以外を含みうる)</param>
        /// <param name="p_窓長">窓の長さ (最短の問い合わせ長以上)</param>
        /// <param name="p_結果">窓の開始位置ごとの答え (長さは 配列の長さ - 窓長 + 1)</param>
        public void V_判定_出現_全窓(ReadOnlySpan<char> p_配列, int p_窓長, Span<bool> p_結果)
        {
            if (p_窓長 < this.A_最短の問い合わせ長)
            {
                throw new ArgumentException($"窓は {this.A_最短の問い合わせ長} 塩基以上が要る", nameof(p_窓長));
            }

            var l_窓数 = p_配列.Length - p_窓長 + 1;
            if (p_結果.Length != l_窓数)
            {
                throw new ArgumentException("結果の長さは 配列の長さ - 窓長 + 1 でなければならない", nameof(p_結果));
            }

            if (l_窓数 <= 0)
            {
                return;
            }

            // 配列全体の 15-mer の正準値とハッシュを一度だけ転がして求める (A・C・G・T 以外を含む 15-mer の値は使われない)
            var l_種数 = p_配列.Length - C_種長 + 1;
            var l_正準値 = new uint[l_種数];
            var l_ハッシュ = new ulong[l_種数];
            ulong l_順 = 0;
            ulong l_逆 = 0;
            var l_連続 = 0;
            for (var p = 0; p < p_配列.Length; p++)
            {
                var l_文字 = Get_2bit値(p_配列[p]);
                if (l_文字 < 0)
                {
                    l_連続 = 0;
                    l_順 = 0;
                    l_逆 = 0;
                    continue;
                }

                l_順 = ((l_順 << 2) | (uint)l_文字) & C_種のマスク;
                l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_文字) << (2 * (C_種長 - 1)));
                l_連続 = Math.Min(l_連続 + 1, C_種長);
                if (l_連続 == C_種長)
                {
                    var l_開始 = p - C_種長 + 1;
                    var l_値 = (uint)Math.Min(l_順, l_逆);
                    l_正準値[l_開始] = l_値;
                    l_ハッシュ[l_開始] = Get_ハッシュ(l_値);
                }
            }

            // 窓の中に A・C・G・T 以外があるかを、窓の始まりごとに次の不正な位置で確かめる
            var l_次の不正 = 0;
            for (var i = 0; i < l_窓数; i++)
            {
                l_次の不正 = Math.Max(l_次の不正, i);
                while (l_次の不正 < p_配列.Length && Get_2bit値(p_配列[l_次の不正]) >= 0)
                {
                    l_次の不正++;
                }

                if (l_次の不正 < i + p_窓長)
                {
                    p_結果[i] = false;
                    continue;
                }

                p_結果[i] = this.Get_出現場所群_種から(p_配列.Slice(i, p_窓長), l_ハッシュ.AsSpan(i, this.A_窓の種数), l_正準値.AsSpan(i, this.A_窓の種数), 1).Count > 0;
            }
        }

        /// <summary>
        /// 配列がリードかその逆相補に出てくる箇所の数を、上限まで数える
        /// </summary>
        /// <param name="p_配列">調べる配列 (最短の問い合わせ長以上、A・C・G・T だけ)</param>
        /// <param name="p_上限">ここまで数えたら打ち切る</param>
        /// <returns>出てくる箇所の数 (上限で頭打ち)</returns>
        public int Get_出現数(ReadOnlySpan<char> p_配列, int p_上限)
        {
            return this.Get_出現場所群(p_配列, p_上限).Count;
        }

        /// <summary>
        /// 配列が出てくるリードごとに、配列の直前と直後の塩基を配列の向きで取り出す (リードの端で止まる)
        /// </summary>
        /// <param name="p_配列">錨にする配列 (最短の問い合わせ長以上、A・C・G・T だけ)</param>
        /// <param name="p_長さ">前後それぞれ取り出す長さの上限</param>
        /// <param name="p_上限">取り出すリードの数の上限</param>
        /// <returns>(直前の塩基列, 直後の塩基列)、無ければ空文字</returns>
        public List<(string A_前, string A_続き)> Get_前後群(ReadOnlySpan<char> p_配列, int p_長さ, int p_上限)
        {
            List<(string A_前, string A_続き)> l_前後群 = [];
            foreach (var (l_開始, l_Is逆鎖) in this.Get_出現場所群(p_配列, p_上限))
            {
                var l_終わり = l_開始 + p_配列.Length;
                var l_右側 = this.Get_外側(l_終わり, 1, p_長さ);
                var l_左側 = this.Get_外側(l_開始 - 1, -1, p_長さ);
                l_前後群.Add(l_Is逆鎖 ? (Get_逆相補(l_右側, 0), Get_逆相補(l_左側, 1)) : (Get_逆相補(l_左側, 2), l_右側));
            }

            return l_前後群;
        }

        /// <summary>
        /// 詰めたリードと区間群は共有したまま、窓の種数だけ変えた索引を作る (minimizer だけを集め直す)
        /// </summary>
        /// <param name="p_窓の種数">minimizer を選ぶ窓に並ぶ k-mer の数 (17 以上)</param>
        /// <returns>窓の種数を変えた索引</returns>
        public ReadMinimizerIndex Get_窓違い(int p_窓の種数)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(p_窓の種数, C_短い問い合わせ用の窓の種数);
            return V_組み立てる(this._語, this._末尾印, this._区間群, this.A_塩基数, p_窓の種数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 詰めたリードと区間群から、指定の窓の種数で minimizer を集めて並べ、索引を組み立てる
        /// </summary>
        /// <param name="p_語">詰めたリード</param>
        /// <param name="p_末尾印">区間の末尾の印</param>
        /// <param name="p_区間群">問い合わせに使える長さの区間 (位置の順)</param>
        /// <param name="p_塩基数">塩基数</param>
        /// <param name="p_窓の種数">minimizer を選ぶ窓に並ぶ k-mer の数</param>
        /// <returns>作った索引</returns>
        private static ReadMinimizerIndex V_組み立てる(ulong[] p_語, ulong[] p_末尾印, (long A_開始, int A_長さ)[] p_区間群, long p_塩基数, int p_窓の種数)
        {
            var l_束数 = (p_区間群.Length + C_区間の束の大きさ - 1) / C_区間の束の大きさ;
            var l_束の先頭 = new long[l_束数 + 1];
            _ = Parallel.For(0, l_束数, l_束 =>
            {
                l_束の先頭[l_束 + 1] = V_集める_束(p_語, p_区間群, l_束, p_窓の種数, null, null, null, 0);
            });
            for (var i = 0; i < l_束数; i++)
            {
                l_束の先頭[i + 1] += l_束の先頭[i];
            }

            var l_種 = new uint[l_束の先頭[l_束数]];
            var l_Is32bit = p_塩基数 <= uint.MaxValue;
            var l_位置32 = l_Is32bit ? new uint[l_種.LongLength] : null;
            var l_位置64 = l_Is32bit ? null : new long[l_種.LongLength];
            _ = Parallel.For(0, l_束数, l_束 =>
            {
                _ = V_集める_束(p_語, p_区間群, l_束, p_窓の種数, l_種, l_位置32, l_位置64, l_束の先頭[l_束]);
            });
            if (l_位置32 is not null)
            {
                V_並べる(ref l_種, ref l_位置32);
            }
            else
            {
                V_並べる(ref l_種, ref l_位置64!);
            }

            return new ReadMinimizerIndex(p_語, p_末尾印, l_種, l_位置32, l_位置64, p_塩基数, p_窓の種数, p_区間群);
        }

        /// <summary>
        /// 配列がリードかその逆相補に出てくる場所を、上限まで集める
        /// </summary>
        /// <param name="p_配列">調べる配列 (最短の問い合わせ長以上、A・C・G・T だけ)</param>
        /// <param name="p_上限">ここまで集めたら打ち切る</param>
        /// <returns>(詰めた配列上の始まり, 逆相補で出てきたか)</returns>
        private List<(long A_開始, bool A_Is逆鎖)> Get_出現場所群(ReadOnlySpan<char> p_配列, int p_上限)
        {
            if (p_配列.Length < this.A_最短の問い合わせ長)
            {
                throw new ArgumentException($"配列は {this.A_最短の問い合わせ長} 塩基以上が要る");
            }

            Span<ulong> l_ハッシュ = stackalloc ulong[this.A_窓の種数];
            Span<uint> l_正準値 = stackalloc uint[this.A_窓の種数];
            for (var j = 0; j < this.A_窓の種数; j++)
            {
                var l_値 = Get_種の値(p_配列.Slice(j, C_種長));
                if (l_値 < 0)
                {
                    return [];
                }

                l_正準値[j] = (uint)l_値;
                l_ハッシュ[j] = Get_ハッシュ((uint)l_値);
            }

            return this.Get_出現場所群_種から(p_配列, l_ハッシュ, l_正準値, p_上限);
        }

        /// <summary>
        /// 窓の先頭の種 (15-mer) の正準値とハッシュから、最小ハッシュの種を選んで配列がリードかその逆相補に出てくる場所を、上限まで集める
        /// </summary>
        /// <param name="p_配列">調べる配列 (最短の問い合わせ長以上、A・C・G・T だけ)</param>
        /// <param name="p_ハッシュ">窓の先頭の A_窓の種数 個の種のハッシュ</param>
        /// <param name="p_正準値">窓の先頭の A_窓の種数 個の種の正準値</param>
        /// <param name="p_上限">ここまで集めたら打ち切る</param>
        /// <returns>(詰めた配列上の始まり, 逆相補で出てきたか)</returns>
        private List<(long A_開始, bool A_Is逆鎖)> Get_出現場所群_種から(ReadOnlySpan<char> p_配列, ReadOnlySpan<ulong> p_ハッシュ, ReadOnlySpan<uint> p_正準値, int p_上限)
        {
            var l_最小 = ulong.MaxValue;
            for (var j = 0; j < this.A_窓の種数; j++)
            {
                l_最小 = Math.Min(l_最小, p_ハッシュ[j]);
            }

            var l_種 = uint.MaxValue;
            for (var j = 0; j < this.A_窓の種数; j++)
            {
                if (p_ハッシュ[j] == l_最小)
                {
                    l_種 = p_正準値[j];
                    break;
                }
            }

            List<(long A_開始, bool A_Is逆鎖)> l_見つけた場所 = [];
            HashSet<(long, bool)> l_見つけた印 = [];
            var l_下 = Get_下限(this._種, l_種);
            for (var l_項 = l_下; l_項 < this._種.LongLength && this._種[l_項] == l_種; l_項++)
            {
                var l_場所 = this._位置32?[l_項] ?? this._位置64![l_項];
                for (var j = 0; j < this.A_窓の種数; j++)
                {
                    if (p_ハッシュ[j] != l_最小)
                    {
                        continue;
                    }

                    var l_逆鎖の開始 = l_場所 - (p_配列.Length - j - C_種長);
                    foreach (var (l_開始, l_Is逆鎖) in (ReadOnlySpan<(long, bool)>)[(l_場所 - j, false), (l_逆鎖の開始, true)])
                    {
                        if (!this.Is一致(l_開始, p_配列, l_Is逆鎖) || !l_見つけた印.Add((l_開始, l_Is逆鎖)))
                        {
                            continue;
                        }

                        l_見つけた場所.Add((l_開始, l_Is逆鎖));
                        if (l_見つけた場所.Count >= p_上限)
                        {
                            return l_見つけた場所;
                        }
                    }
                }
            }

            return l_見つけた場所;
        }

        /// <summary>
        /// 位置から向きへ、区間 (リード) の端まで塩基を読む
        /// </summary>
        /// <param name="p_位置">読み始める位置</param>
        /// <param name="p_向き">+1 なら右へ、-1 なら左へ</param>
        /// <param name="p_長さ">読む長さの上限</param>
        /// <returns>読んだ順の塩基列</returns>
        private string Get_外側(long p_位置, int p_向き, int p_長さ)
        {
            var l_塩基 = new char[p_長さ];
            var l_数 = 0;
            for (var l_位置 = p_位置; l_数 < p_長さ && l_位置 >= 0 && l_位置 < this.A_塩基数 && !(p_向き > 0 ? this.Is区間の末尾(l_位置 - 1) : this.Is区間の末尾(l_位置)); l_位置 += p_向き)
            {
                l_塩基[l_数++] = C_塩基の並び[Get_文字(this._語, l_位置)];
            }

            return new string(l_塩基, 0, l_数);
        }

        /// <summary>
        /// 読んだ塩基列を、錨の向きの並びに直す
        /// </summary>
        /// <param name="p_塩基列">読んだ順の塩基列</param>
        /// <param name="p_直し方">0: 右へ読んだものを逆相補 (逆鎖の直前)、1: 左へ読んだものを相補 (逆鎖の直後)、2: 左へ読んだものを逆順 (順鎖の直前)</param>
        /// <returns></returns>
        private static string Get_逆相補(string p_塩基列, int p_直し方)
        {
            var l_結果 = new char[p_塩基列.Length];
            for (var i = 0; i < p_塩基列.Length; i++)
            {
                l_結果[i] = p_直し方 switch
                {
                    0 => Util.Get_相補塩基(p_塩基列[p_塩基列.Length - 1 - i]),
                    1 => Util.Get_相補塩基(p_塩基列[i]),
                    _ => p_塩基列[p_塩基列.Length - 1 - i],
                };
            }

            return new string(l_結果);
        }

        /// <summary>
        /// 位置が区間 (リードの A・C・G・T の続き) の最後の文字か
        /// </summary>
        /// <param name="p_位置"></param>
        /// <returns></returns>
        private bool Is区間の末尾(long p_位置)
        {
            return (this._末尾印[p_位置 / 64] & (1UL << (int)(p_位置 % 64))) != 0;
        }

        /// <summary>
        /// 束ねた配列を並列に 2 bit へ詰める (A・C・G・T が続く区間ごとに詰め、区間の境目で語を分け合うので、語へは OR で書く)
        /// </summary>
        /// <param name="p_語">詰める先 (足りなければ広げる)</param>
        /// <param name="p_末尾印">区間の末尾の印 (足りなければ広げる)</param>
        /// <param name="p_区間群">問い合わせに使える長さの区間を、位置の順に足す先</param>
        /// <param name="p_束">配列</param>
        /// <param name="p_位置">束の最初の配列を詰め始める位置</param>
        /// <returns>束の最後の配列の直後の位置</returns>
        private static long V_詰める_束(ref ulong[] p_語, ref ulong[] p_末尾印, List<(long A_開始, int A_長さ)> p_区間群, List<string> p_束, long p_位置)
        {
            var l_件数 = p_束.Count;
            var l_開始位置 = new long[l_件数 + 1];
            var l_長い区間の先頭 = new int[l_件数 + 1];
            _ = Parallel.For(0, l_件数, i =>
            {
                var (l_塩基数, l_長い区間数) = Get_塩基数と長い区間数(p_束[i]);
                l_開始位置[i + 1] = l_塩基数;
                l_長い区間の先頭[i + 1] = l_長い区間数;
            });
            l_開始位置[0] = p_位置;
            for (var i = 0; i < l_件数; i++)
            {
                l_開始位置[i + 1] += l_開始位置[i];
                l_長い区間の先頭[i + 1] += l_長い区間の先頭[i];
            }

            var l_終端 = l_開始位置[l_件数];
            var l_要る語数 = (l_終端 / C_語あたりの文字数) + 2;
            if (p_語.LongLength < l_要る語数)
            {
                Array.Resize(ref p_語, (int)Math.Max(l_要る語数, Math.Min(Array.MaxLength, p_語.LongLength * 2)));
            }

            var l_要る印数 = (l_終端 / 64) + 2;
            if (p_末尾印.LongLength < l_要る印数)
            {
                Array.Resize(ref p_末尾印, (int)Math.Max(l_要る印数, Math.Min(Array.MaxLength, p_末尾印.LongLength * 2)));
            }

            var l_語 = p_語;
            var l_末尾印 = p_末尾印;
            var l_長い区間 = new (long A_開始, int A_長さ)[l_長い区間の先頭[l_件数]];
            _ = Parallel.For(0, l_件数, i =>
            {
                var l_配列 = p_束[i];
                var l_今 = l_開始位置[i];
                var l_書く = l_長い区間の先頭[i];
                var l_区間長 = 0;
                var l_語番号 = l_今 / C_語あたりの文字数;
                var l_溜め = 0UL;
                for (var j = 0; j <= l_配列.Length; j++)
                {
                    var l_値 = j < l_配列.Length ? Get_2bit値(l_配列[j]) : -1;
                    if (l_値 < 0)
                    {
                        if (l_区間長 > 0)
                        {
                            var l_末尾 = l_今 - 1;
                            _ = Interlocked.Or(ref l_末尾印[l_末尾 / 64], 1UL << (int)(l_末尾 % 64));
                            if (l_区間長 >= C_短い問い合わせの最短長)
                            {
                                l_長い区間[l_書く++] = (l_今 - l_区間長, l_区間長);
                            }

                            l_区間長 = 0;
                        }

                        continue;
                    }

                    if (l_今 / C_語あたりの文字数 != l_語番号)
                    {
                        _ = Interlocked.Or(ref l_語[l_語番号], l_溜め);
                        l_溜め = 0UL;
                        l_語番号 = l_今 / C_語あたりの文字数;
                    }

                    l_溜め |= (ulong)l_値 << (62 - (int)(2 * (l_今 % C_語あたりの文字数)));
                    l_今++;
                    l_区間長++;
                }

                if (l_溜め != 0UL)
                {
                    _ = Interlocked.Or(ref l_語[l_語番号], l_溜め);
                }
            });
            p_区間群.AddRange(l_長い区間);
            return l_終端;
        }

        /// <summary>
        /// 配列の A・C・G・T の数と、問い合わせに使える長さの区間の数
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <returns>塩基数と長い区間の数</returns>
        private static (long A_塩基数, int A_長い区間数) Get_塩基数と長い区間数(string p_配列)
        {
            var l_塩基数 = 0L;
            var l_長い区間数 = 0;
            var l_区間長 = 0;
            foreach (var l_文字 in p_配列)
            {
                if (Get_2bit値(l_文字) >= 0)
                {
                    l_塩基数++;
                    l_区間長++;
                    continue;
                }

                l_長い区間数 += l_区間長 >= C_短い問い合わせの最短長 ? 1 : 0;
                l_区間長 = 0;
            }

            l_長い区間数 += l_区間長 >= C_短い問い合わせの最短長 ? 1 : 0;
            return (l_塩基数, l_長い区間数);
        }

        /// <summary>
        /// 種の値の順に並べる (上位のビットで振り分けてから、振り分けた先ごとに並列に並べる)
        /// </summary>
        /// <typeparam name="T">位置の型</typeparam>
        /// <param name="p_種">種の値</param>
        /// <param name="p_位置">種と組になる位置</param>
        private static void V_並べる<T>(ref uint[] p_種, ref T[] p_位置)
        {
            var l_振り分けのシフト = Math.Max(0, (2 * C_種長) - C_振り分けのビット数);
            var l_先頭 = new long[(1 << C_振り分けのビット数) + 1];
            foreach (var l_値 in p_種)
            {
                l_先頭[(l_値 >> l_振り分けのシフト) + 1]++;
            }

            for (var i = 0; i < l_先頭.Length - 1; i++)
            {
                l_先頭[i + 1] += l_先頭[i];
            }

            var l_種 = new uint[p_種.LongLength];
            var l_位置 = new T[p_位置.LongLength];
            var l_書く = (long[])l_先頭.Clone();
            for (long i = 0; i < p_種.LongLength; i++)
            {
                var l_先 = l_書く[p_種[i] >> l_振り分けのシフト]++;
                l_種[l_先] = p_種[i];
                l_位置[l_先] = p_位置[i];
            }

            _ = Parallel.For(0, l_先頭.Length - 1, b =>
            {
                var l_長さ = (int)(l_先頭[b + 1] - l_先頭[b]);
                if (l_長さ > 1)
                {
                    l_種.AsSpan((int)l_先頭[b], l_長さ).Sort(l_位置.AsSpan((int)l_先頭[b], l_長さ));
                }
            });

            p_種 = l_種;
            p_位置 = l_位置;
        }

        /// <summary>
        /// 束に入る区間の minimizer を数える、または書き込む
        /// </summary>
        /// <param name="p_語">詰めたリード</param>
        /// <param name="p_区間群">区間の開始位置と長さ</param>
        /// <param name="p_束">束の番号</param>
        /// <param name="p_窓の種数">minimizer を選ぶ窓の種数</param>
        /// <param name="p_種">書き込み先 (null なら数えるだけ)</param>
        /// <param name="p_位置32">位置の書き込み先 (32 bit)</param>
        /// <param name="p_位置64">位置の書き込み先 (64 bit)</param>
        /// <param name="p_書き始め">書き込みを始める位置</param>
        /// <returns>minimizer の数</returns>
        private static long V_集める_束(ulong[] p_語, (long A_開始, int A_長さ)[] p_区間群, int p_束, int p_窓の種数, uint[]? p_種, uint[]? p_位置32, long[]? p_位置64, long p_書き始め)
        {
            var l_書く = p_書き始め;
            var l_終わり = Math.Min(p_区間群.Length, (p_束 + 1) * C_区間の束の大きさ);
            var l_最短 = C_種長 + p_窓の種数 - 1;
            for (var i = p_束 * C_区間の束の大きさ; i < l_終わり; i++)
            {
                if (p_区間群[i].A_長さ < l_最短)
                {
                    // 窓が 1 つもできない区間は minimizer を持たない
                    continue;
                }

                V_集める_minimizer(p_語, p_区間群[i].A_開始, p_区間群[i].A_長さ, p_窓の種数, (l_値, l_場所) =>
                {
                    if (p_種 is not null)
                    {
                        p_種[l_書く] = l_値;
                        if (p_位置32 is not null)
                        {
                            p_位置32[l_書く] = (uint)l_場所;
                        }
                        else
                        {
                            p_位置64![l_書く] = l_場所;
                        }
                    }

                    l_書く++;
                });
            }

            return l_書く - p_書き始め;
        }

        /// <summary>
        /// 区間の minimizer (窓ごとにハッシュ最小で左端のもの) を重複なしで順に渡す
        /// </summary>
        /// <param name="p_語">詰めたリード</param>
        /// <param name="p_開始">区間の開始位置</param>
        /// <param name="p_長さ">区間の長さ</param>
        /// <param name="p_窓の種数">minimizer を選ぶ窓の種数</param>
        /// <param name="p_渡し先">minimizer の正準値と位置を受け取る処理</param>
        private static void V_集める_minimizer(ulong[] p_語, long p_開始, int p_長さ, int p_窓の種数, Action<uint, long> p_渡し先)
        {
            var l_種数 = p_長さ - C_種長 + 1;
            var l_正準値 = new uint[l_種数];
            var l_ハッシュ = new ulong[l_種数];
            ulong l_順 = 0;
            ulong l_逆 = 0;
            for (var i = 0; i < p_長さ; i++)
            {
                var l_文字 = Get_文字(p_語, p_開始 + i);
                l_順 = ((l_順 << 2) | (uint)l_文字) & C_種のマスク;
                l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_文字) << (2 * (C_種長 - 1)));
                if (i >= C_種長 - 1)
                {
                    var l_値 = (uint)Math.Min(l_順, l_逆);
                    l_正準値[i - C_種長 + 1] = l_値;
                    l_ハッシュ[i - C_種長 + 1] = Get_ハッシュ(l_値);
                }
            }

            var l_前 = -1;
            for (var l_窓 = 0; l_窓 + p_窓の種数 <= l_種数; l_窓++)
            {
                var l_最小位置 = l_窓;
                for (var j = l_窓 + 1; j < l_窓 + p_窓の種数; j++)
                {
                    if (l_ハッシュ[j] < l_ハッシュ[l_最小位置])
                    {
                        l_最小位置 = j;
                    }
                }

                if (l_最小位置 != l_前)
                {
                    p_渡し先(l_正準値[l_最小位置], p_開始 + l_最小位置);
                    l_前 = l_最小位置;
                }
            }
        }

        /// <summary>
        /// 位置から配列 (またはその逆相補) がそのまま並んでいるか (区間をまたがない)
        /// </summary>
        /// <param name="p_開始">照合を始める位置</param>
        /// <param name="p_配列">配列</param>
        /// <param name="p_Is逆相補">逆相補と照合するか</param>
        /// <returns></returns>
        private bool Is一致(long p_開始, ReadOnlySpan<char> p_配列, bool p_Is逆相補)
        {
            if (p_開始 < 0 || p_開始 + p_配列.Length > this.A_塩基数)
            {
                return false;
            }

            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_位置 = p_開始 + i;
                var l_期待 = p_Is逆相補 ? 3 - Get_2bit値(p_配列[p_配列.Length - 1 - i]) : Get_2bit値(p_配列[i]);
                if (Get_文字(this._語, l_位置) != l_期待)
                {
                    return false;
                }

                if (i < p_配列.Length - 1 && (this._末尾印[l_位置 / 64] & (1UL << (int)(l_位置 % 64))) != 0)
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// 並んだ値のうち、値以上になる最初の位置
        /// </summary>
        /// <param name="p_並び">昇順の値</param>
        /// <param name="p_値">探す値</param>
        /// <returns></returns>
        private static long Get_下限(uint[] p_並び, uint p_値)
        {
            long l_下 = 0;
            var l_上 = p_並び.LongLength;
            while (l_下 < l_上)
            {
                var l_中 = (l_下 + l_上) / 2;
                if (p_並び[l_中] < p_値)
                {
                    l_下 = l_中 + 1;
                }
                else
                {
                    l_上 = l_中;
                }
            }

            return l_下;
        }

        /// <summary>
        /// 種の正準値 (逆相補と小さいほう)
        /// </summary>
        /// <param name="p_種">種長の配列</param>
        /// <returns>正準値、A・C・G・T 以外を含めば -1</returns>
        private static long Get_種の値(ReadOnlySpan<char> p_種)
        {
            ulong l_順 = 0;
            ulong l_逆 = 0;
            for (var i = 0; i < p_種.Length; i++)
            {
                var l_文字 = Get_2bit値(p_種[i]);
                if (l_文字 < 0)
                {
                    return -1;
                }

                l_順 = (l_順 << 2) | (uint)l_文字;
                l_逆 |= (ulong)(3 - l_文字) << (2 * i);
            }

            return (long)Math.Min(l_順, l_逆);
        }

        /// <summary>
        /// 種の正準値を並べ替える順 (値ごとに異なる、並びの偏りを崩すための混ぜ合わせ)
        /// </summary>
        /// <param name="p_値">正準値</param>
        /// <returns></returns>
        private static ulong Get_ハッシュ(uint p_値)
        {
            var l_値 = (ulong)p_値;
            l_値 ^= l_値 >> 33;
            l_値 *= 0xFF51_AFD7_ED55_8CCDUL;
            l_値 ^= l_値 >> 33;
            l_値 *= 0xC4CE_B9FE_1A85_EC53UL;
            l_値 ^= l_値 >> 33;
            return l_値;
        }

        /// <summary>
        /// 2 bit で詰めた配列の位置の文字
        /// </summary>
        /// <param name="p_配列">詰めた配列</param>
        /// <param name="p_位置">位置</param>
        /// <returns></returns>
        private static int Get_文字(ulong[] p_配列, long p_位置)
        {
            return (int)((p_配列[p_位置 / C_語あたりの文字数] >> (62 - (int)(2 * (p_位置 % C_語あたりの文字数)))) & 3);
        }

        /// <summary>
        /// 塩基 1 文字の 2 bit の値
        /// </summary>
        /// <param name="p_塩基">塩基</param>
        /// <returns>A・C・G・T なら 0〜3、それ以外は -1</returns>
        private static int Get_2bit値(char p_塩基)
        {
            return p_塩基 switch
            {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => -1,
            };
        }

        #endregion
    }
}
