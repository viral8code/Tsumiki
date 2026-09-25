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
        /// 1 語に詰める文字数
        /// </summary>
        private const int C_語あたりの文字数 = 32;

        /// <summary>
        /// minimizer を並列に集めるときの 1 束の区間数
        /// </summary>
        private const int C_区間の束の大きさ = 4_096;

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

        #endregion

        #region プロパティ

        /// <summary>
        /// 索引に入れた塩基数
        /// </summary>
        public long A_塩基数 { get; }

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
        private ReadMinimizerIndex(ulong[] p_語, ulong[] p_末尾印, uint[] p_種, uint[]? p_位置32, long[]? p_位置64, long p_塩基数)
        {
            this._語 = p_語;
            this._末尾印 = p_末尾印;
            this._種 = p_種;
            this._位置32 = p_位置32;
            this._位置64 = p_位置64;
            this.A_塩基数 = p_塩基数;
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群から索引を作る
        /// </summary>
        /// <param name="p_配列列">配列を先頭から返すもの (2 回呼ぶ)</param>
        /// <returns>作った索引</returns>
        public static ReadMinimizerIndex V_構築(Func<IEnumerable<string>> p_配列列)
        {
            long l_長さ = 0;
            foreach (var l_配列 in p_配列列())
            {
                foreach (var (_, l_区間長) in Get_区間列(l_配列))
                {
                    l_長さ += l_区間長;
                }
            }

            var l_語 = new ulong[(l_長さ / C_語あたりの文字数) + 2];
            var l_末尾印 = new ulong[(l_長さ / 64) + 2];
            List<(long A_開始, int A_長さ)> l_区間群 = [];
            long l_位置 = 0;
            foreach (var l_配列 in p_配列列())
            {
                foreach (var (l_開始, l_区間長) in Get_区間列(l_配列))
                {
                    for (var i = 0; i < l_区間長; i++)
                    {
                        l_語[l_位置 / C_語あたりの文字数] |= (ulong)Get_2bit値(l_配列[l_開始 + i]) << (62 - (int)(2 * (l_位置 % C_語あたりの文字数)));
                        l_位置++;
                    }

                    l_末尾印[(l_位置 - 1) / 64] |= 1UL << (int)((l_位置 - 1) % 64);
                    if (l_区間長 >= C_最短の問い合わせ長)
                    {
                        l_区間群.Add((l_位置 - l_区間長, l_区間長));
                    }
                }
            }

            var l_束数 = (l_区間群.Count + C_区間の束の大きさ - 1) / C_区間の束の大きさ;
            var l_束の先頭 = new long[l_束数 + 1];
            _ = Parallel.For(0, l_束数, l_束 =>
            {
                l_束の先頭[l_束 + 1] = V_集める_束(l_語, l_区間群, l_束, null, null, null, 0);
            });
            for (var i = 0; i < l_束数; i++)
            {
                l_束の先頭[i + 1] += l_束の先頭[i];
            }

            var l_種 = new uint[l_束の先頭[l_束数]];
            var l_Is32bit = l_長さ <= uint.MaxValue;
            var l_位置32 = l_Is32bit ? new uint[l_種.LongLength] : null;
            var l_位置64 = l_Is32bit ? null : new long[l_種.LongLength];
            _ = Parallel.For(0, l_束数, l_束 =>
            {
                _ = V_集める_束(l_語, l_区間群, l_束, l_種, l_位置32, l_位置64, l_束の先頭[l_束]);
            });
            if (l_位置32 is not null)
            {
                Array.Sort(l_種, l_位置32);
            }
            else
            {
                Array.Sort(l_種, l_位置64);
            }

            return new ReadMinimizerIndex(l_語, l_末尾印, l_種, l_位置32, l_位置64, l_長さ);
        }

        /// <summary>
        /// 配列がリードかその逆相補に 1 回以上出てくるか
        /// </summary>
        /// <param name="p_配列">調べる配列 (最短の問い合わせ長以上、A・C・G・T だけ)</param>
        /// <returns>出てくれば true</returns>
        public bool Has出現(ReadOnlySpan<char> p_配列)
        {
            if (p_配列.Length < C_最短の問い合わせ長)
            {
                throw new ArgumentException($"配列は {C_最短の問い合わせ長} 塩基以上が要る");
            }

            Span<ulong> l_ハッシュ = stackalloc ulong[C_窓の種数];
            Span<uint> l_正準値 = stackalloc uint[C_窓の種数];
            for (var j = 0; j < C_窓の種数; j++)
            {
                var l_値 = Get_種の値(p_配列.Slice(j, C_種長));
                if (l_値 < 0)
                {
                    return false;
                }

                l_正準値[j] = (uint)l_値;
                l_ハッシュ[j] = Get_ハッシュ((uint)l_値);
            }

            var l_最小 = ulong.MaxValue;
            for (var j = 0; j < C_窓の種数; j++)
            {
                l_最小 = Math.Min(l_最小, l_ハッシュ[j]);
            }

            var l_種 = uint.MaxValue;
            for (var j = 0; j < C_窓の種数; j++)
            {
                if (l_ハッシュ[j] == l_最小)
                {
                    l_種 = l_正準値[j];
                    break;
                }
            }

            var l_下 = Get_下限(this._種, l_種);
            for (var l_項 = l_下; l_項 < this._種.LongLength && this._種[l_項] == l_種; l_項++)
            {
                var l_場所 = this._位置32?[l_項] ?? this._位置64![l_項];
                for (var j = 0; j < C_窓の種数; j++)
                {
                    if (l_ハッシュ[j] != l_最小)
                    {
                        continue;
                    }

                    if (this.Is一致(l_場所 - j, p_配列, false) || this.Is一致(l_場所 - (p_配列.Length - j - C_種長), p_配列, true))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 束に入る区間の minimizer を数える、または書き込む
        /// </summary>
        /// <param name="p_語">詰めたリード</param>
        /// <param name="p_区間群">区間の開始位置と長さ</param>
        /// <param name="p_束">束の番号</param>
        /// <param name="p_種">書き込み先 (null なら数えるだけ)</param>
        /// <param name="p_位置32">位置の書き込み先 (32 bit)</param>
        /// <param name="p_位置64">位置の書き込み先 (64 bit)</param>
        /// <param name="p_書き始め">書き込みを始める位置</param>
        /// <returns>minimizer の数</returns>
        private static long V_集める_束(ulong[] p_語, List<(long A_開始, int A_長さ)> p_区間群, int p_束, uint[]? p_種, uint[]? p_位置32, long[]? p_位置64, long p_書き始め)
        {
            var l_書く = p_書き始め;
            var l_終わり = Math.Min(p_区間群.Count, (p_束 + 1) * C_区間の束の大きさ);
            for (var i = p_束 * C_区間の束の大きさ; i < l_終わり; i++)
            {
                V_集める_minimizer(p_語, p_区間群[i].A_開始, p_区間群[i].A_長さ, (l_値, l_場所) =>
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
        /// <param name="p_渡し先">minimizer の正準値と位置を受け取る処理</param>
        private static void V_集める_minimizer(ulong[] p_語, long p_開始, int p_長さ, Action<uint, long> p_渡し先)
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
            for (var l_窓 = 0; l_窓 + C_窓の種数 <= l_種数; l_窓++)
            {
                var l_最小位置 = l_窓;
                for (var j = l_窓 + 1; j < l_窓 + C_窓の種数; j++)
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
        /// 配列のうち A・C・G・T が続く区間
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <returns>区間の開始位置と長さ</returns>
        private static IEnumerable<(int A_開始, int A_長さ)> Get_区間列(string p_配列)
        {
            var l_開始 = -1;
            for (var i = 0; i <= p_配列.Length; i++)
            {
                var l_Is塩基 = i < p_配列.Length && Get_2bit値(p_配列[i]) >= 0;
                if (l_Is塩基 && l_開始 < 0)
                {
                    l_開始 = i;
                }
                else if (!l_Is塩基 && l_開始 >= 0)
                {
                    yield return (l_開始, i - l_開始);
                    l_開始 = -1;
                }
            }
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
