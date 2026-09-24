using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 短い反復配列の経路付け替えに対する拒否権 (ABySS RResolver 型)
    /// </summary>
    internal sealed class RepeatRMerVerifier
    {
        #region 定数

        /// <summary>
        /// 集合を分割する数のビット数
        /// </summary>
        private const int 分割のビット数 = 6;

        /// <summary>
        /// 集合の分割数
        /// </summary>
        /// <remarks>
        /// 分割ごとに錠を取ることで、リードの走査を並列にしても登録が 1 本の錠に詰まらない
        /// </remarks>
        private const int 分割数 = 1 << 分割のビット数;

        #endregion

        #region 内部変数

        /// <summary>
        /// 128 塩基までの厳密キー集合
        /// </summary>
        private readonly HashSet<(UInt128 A_上位, UInt128 A_下位)>[]? _長集合;

        /// <summary>
        /// 128 塩基を超える厳密キー集合
        /// </summary>
        private readonly HashSet<KmerKey>[]? _大集合;

        /// <summary>
        /// 32 塩基までの厳密キー集合
        /// </summary>
        private readonly HashSet<ulong>[]? _小集合;

        /// <summary>
        /// 64 塩基までの厳密キー集合
        /// </summary>
        private readonly HashSet<UInt128>[]? _中集合;

        /// <summary>
        /// 登録を許す r-mer、null なら全部登録する
        /// </summary>
        /// <remarks>
        /// 構築が終わった後は読むだけなので、走査を並列にしても錠は要らない
        /// </remarks>
        private RepeatRMerVerifier? _候補集合;

        /// <summary>
        /// 分割ごとの錠
        /// </summary>
        private readonly Lock[] _分割錠;

        /// <summary>
        /// r-mer の長さ
        /// </summary>
        private readonly int _r長;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// コンストラクタ
        /// </summary>
        /// <param name="p_r長">r-mer の長さ</param>
        private RepeatRMerVerifier(int p_r長)
        {
            this._r長 = p_r長;
            this._分割錠 = [.. Enumerable.Range(0, 分割数).Select(_ => new Lock())];
            if (p_r長 <= 32)
            {
                this._小集合 = [.. Enumerable.Range(0, 分割数).Select(_ => new HashSet<ulong>())];
            }
            else if (p_r長 <= 64)
            {
                this._中集合 = [.. Enumerable.Range(0, 分割数).Select(_ => new HashSet<UInt128>())];
            }
            else if (p_r長 <= 128)
            {
                this._長集合 = [.. Enumerable.Range(0, 分割数).Select(_ => new HashSet<(UInt128 A_上位, UInt128 A_下位)>())];
            }
            else
            {
                this._大集合 = [.. Enumerable.Range(0, 分割数).Select(_ => new HashSet<KmerKey>())];
            }
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 生リードファイル群を 1 回走査し、出現した r-mer (正準形) の集合を作る
        /// </summary>
        /// <remarks>
        /// 空、存在しないパスは片側リードのみの実行に対応するため無視する<br/>
        /// 検査する窓はこの k のグラフ上の経路なので、中の k-mer はすべて p_kmerインデックス にある<br/>
        /// 両端の k-mer が集合に無い r-mer は照合されることがなく、登録を省いても判定は変わらない<br/>
        /// 省けるのは主にリードのエラーを含む r-mer で、ゲノムの数十倍に膨らむ集合がゲノム規模に収まる<br/>
        /// p_問い合わせ配列 を渡す場合、後の問い合わせがその配列の r-mer に限られることを呼び出し側が保証すること
        /// </remarks>
        /// <param name="p_リードパス一覧">走査するリードファイルのパス</param>
        /// <param name="p_r長">r-mer の長さ</param>
        /// <param name="p_kmerインデックス">登録する r-mer を両端の k-mer で絞り込む集合、null なら絞り込まない</param>
        /// <param name="p_k長">p_kmerインデックス の k</param>
        /// <param name="p_問い合わせ配列">後でこの検証器に問い合わせる配列、渡すとその r-mer だけを登録する</param>
        /// <returns>構築した検証器</returns>
        public static RepeatRMerVerifier V_構築(IEnumerable<string> p_リードパス一覧, int p_r長, TrustedKmerIndex? p_kmerインデックス = null, int p_k長 = 0, IEnumerable<string>? p_問い合わせ配列 = null)
        {
            using var l_計測 = new StageTimer($"repeat-index r={p_r長}");
            if (p_r長 <= 0)
            {
                throw new ArgumentException("r-mer length must be positive");
            }

            var l_パス群 = p_リードパス一覧
                .Where(中間データ置き場.Is存在)
                .ToList();

            var l_検証器 = new RepeatRMerVerifier(p_r長);
            if (p_問い合わせ配列 is not null)
            {
                var l_候補 = new RepeatRMerVerifier(p_r長);
                foreach (var l_配列 in p_問い合わせ配列)
                {
                    l_候補.V_登録_rMer(l_配列, p_r長, null, 0);
                }
                l_検証器._候補集合 = l_候補;
            }
            var l_絞り込み = p_k長 > 0 && p_k長 < p_r長 ? p_kmerインデックス : null;
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            foreach (var l_パス in l_パス群)
            {
                ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, FastqReader.Get_生リード列(l_パス), (l_リード, _) =>
                {
                    if (l_リード is not null)
                    {
                        l_検証器.V_登録_rMer(l_リード, p_r長, l_絞り込み, p_k長);
                    }
                });
            }
            return l_検証器;
        }

        /// <summary>
        /// 接合点の支持が閾値に届いているか
        /// </summary>
        /// <param name="p_head配列">接合点の手前の配列</param>
        /// <param name="p_repeat配列">反復配列</param>
        /// <param name="p_tail配列">接合点の先の配列</param>
        /// <param name="p_閾値">支持しているとみなす r-mer の本数</param>
        /// <returns>閾値に届いていれば true</returns>
        public bool Has接合点支持(string p_head配列, string p_repeat配列, string p_tail配列, int p_閾値)
        {
            var (A_全体, A_入口, A_出口) = this.Get_接合点別支持数(p_head配列, p_repeat配列, p_tail配列);
            return A_入口 > 0 && A_出口 > 0 && A_全体 >= p_閾値;
        }

        /// <summary>
        /// 配列のうち、リードで一度も観測されていない r-mer が続く範囲を集める
        /// </summary>
        /// <param name="p_配列">調べる配列</param>
        /// <param name="p_連続の下限">未観測の窓がこの数だけ続いたら範囲として拾う</param>
        /// <param name="p_範囲">見つけた範囲 (窓の開始位置、終了位置は含む) の書き留め先</param>
        /// <remarks>
        /// 低カバレッジでは未観測の窓が散発するので、連続した長さで反復由来の継ぎ目と区別する<br/>
        /// 曖昧塩基を含む窓は判定できないため観測済みとして扱い、連続を切る
        /// </remarks>
        public void V_収集_未観測の連続範囲(string p_配列, int p_連続の下限, List<(int A_開始, int A_終了)> p_範囲)
        {
            if (p_連続の下限 <= 0 || p_配列.Length < this._r長)
            {
                return;
            }

            var l_直近の曖昧位置 = -1;
            var l_連続開始 = -1;
            for (var i = 0; i + this._r長 <= p_配列.Length; i++)
            {
                var l_新規末尾 = i + this._r長 - 1;
                if (i == 0)
                {
                    for (var j = 0; j < this._r長; j++)
                    {
                        if (Util.Is曖昧塩基(p_配列[j]))
                        {
                            l_直近の曖昧位置 = j;
                        }
                    }
                }
                else if (Util.Is曖昧塩基(p_配列[l_新規末尾]))
                {
                    l_直近の曖昧位置 = l_新規末尾;
                }

                var l_Is未観測 = l_直近の曖昧位置 < i && !this.Has観測(p_配列.AsSpan(i, this._r長));
                if (l_Is未観測)
                {
                    if (l_連続開始 < 0)
                    {
                        l_連続開始 = i;
                    }
                    continue;
                }

                if (l_連続開始 >= 0 && i - l_連続開始 >= p_連続の下限)
                {
                    p_範囲.Add((l_連続開始, i - 1));
                }
                l_連続開始 = -1;
            }

            var l_窓数 = p_配列.Length - this._r長 + 1;
            if (l_連続開始 >= 0 && l_窓数 - l_連続開始 >= p_連続の下限)
            {
                p_範囲.Add((l_連続開始, l_窓数 - 1));
            }
        }

        /// <summary>
        /// 配列とその逆相補のうち、順鎖と逆鎖どちらから読んでも同一になるキーを返す
        /// </summary>
        /// <remarks>
        /// 128 塩基を超える長さは 2 語に収まらないため <see cref="KmerKey"/> を使う
        /// </remarks>
        /// <param name="p_配列">元の配列</param>
        /// <returns>正準化したキー</returns>
        public static (UInt128 A_上位, UInt128 A_下位) Get_正準値(ReadOnlySpan<char> p_配列)
        {
            ArgumentOutOfRangeException.ThrowIfGreaterThan(p_配列.Length, 128);

            var l_Is逆鎖 = false;
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_順 = Util.Get_塩基ID(p_配列[i]);
                var l_逆 = 5 - Util.Get_塩基ID(p_配列[p_配列.Length - 1 - i]);
                if (l_順 != l_逆)
                {
                    l_Is逆鎖 = l_逆 < l_順;
                    break;
                }
            }
            var l_上位 = (UInt128)0;
            var l_下位 = (UInt128)0;
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_値 = l_Is逆鎖 ? 4 - Util.Get_塩基ID(p_配列[p_配列.Length - 1 - i]) : Util.Get_塩基ID(p_配列[i]) - 1;
                l_上位 = (l_上位 << 2) | (l_下位 >> 126);
                l_下位 = (l_下位 << 2) | (uint)l_値;
            }
            return (l_上位, l_下位);
        }

        #endregion

        #region テストメソッド

        /// <summary>
        /// head → repeat → tail の経路上で、head-repeat 接合点と repeat-tail 接合点の両方を実際に跨いだ r-mer のうち、生リード由来の集合に見つかった本数を返す
        /// </summary>
        /// <param name="p_head配列">接合点の手前の配列</param>
        /// <param name="p_repeat配列">反復配列</param>
        /// <param name="p_tail配列">接合点の先の配列</param>
        /// <returns>支持している r-mer の本数</returns>
        public int Get_接合点の支持数(string p_head配列, string p_repeat配列, string p_tail配列)
        {
            return this.Get_接合点別支持数(p_head配列, p_repeat配列, p_tail配列).A_全体;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// head → repeat → tail の経路上で、head-repeat 接合点と repeat-tail 接合点の両方を実際に跨いだ r-mer のうち、生リード由来の集合に見つかった本数を返す
        /// </summary>
        /// <param name="p_head配列">接合点の手前の配列</param>
        /// <param name="p_repeat配列">反復配列</param>
        /// <param name="p_tail配列">接合点の先の配列</param>
        /// <returns>支持している r-mer の本数</returns>
        private (int A_全体, int A_入口, int A_出口) Get_接合点別支持数(string p_head配列, string p_repeat配列, string p_tail配列)
        {
            // head-repeat、repeat-tail は de Bruijn グラフの辺である以上、
            // 必ず (このアセンブリの) k-1 塩基を共有しており、そのコピーは repeat 自身の配列の両端にもそのまま現れる
            // head と tail をそのまま margin に使うと、接合点に重なり区間が二重に並ぶ
            // (head 側のコピーの直後に repeat 自身のコピー) テスト配列を作ってしまい、
            // 本物のゲノム配列 (重なりは一度しか現れない) には存在しない配列になるため、
            // head と tail からはこの重なりを除いた固有部分だけを margin に使う
            var l_重なり長 = Math.Max(0, ConfigurationManager.A_実行時引数.A_k長 - 1);
            var l_head固有 = p_head配列.Length > l_重なり長 ? p_head配列[..^l_重なり長] : string.Empty;
            var l_tail固有 = p_tail配列.Length > l_重なり長 ? p_tail配列[l_重なり長..] : string.Empty;

            var l_margin長 = this._r長 - 1;
            var l_head側 = l_head固有.Length <= l_margin長 ? l_head固有 : l_head固有[^l_margin長..];
            var l_tail側 = l_tail固有.Length <= l_margin長 ? l_tail固有 : l_tail固有[..l_margin長];

            var l_テスト配列 = l_head側 + p_repeat配列 + l_tail側;
            var l_接合点1 = l_head側.Length;
            var l_接合点2 = l_head側.Length + p_repeat配列.Length;

            var l_支持数 = 0;
            var l_入口支持数 = 0;
            var l_出口支持数 = 0;
            // 窓ごとに Is曖昧塩基 を r 回呼ぶと O (n*r) になるため、
            // V_登録_rMer と同じく「直近に見た曖昧塩基の位置」を
            // 窓のスライドに合わせて償却 O (n) で更新する
            // (接合点を跨がない窓は continue するが、曖昧判定は
            // スキップせず続けないと以降の窓の判定がずれる)
            var l_直近の曖昧位置 = -1;
            for (var i = 0; i + this._r長 <= l_テスト配列.Length; i++)
            {
                var l_窓終端 = i + this._r長; // exclusive
                var l_新規末尾 = l_窓終端 - 1;
                if (i == 0)
                {
                    for (var j = 0; j < this._r長; j++)
                    {
                        if (Util.Is曖昧塩基(l_テスト配列[j]))
                        {
                            l_直近の曖昧位置 = j;
                        }
                    }
                }
                else if (Util.Is曖昧塩基(l_テスト配列[l_新規末尾]))
                {
                    l_直近の曖昧位置 = l_新規末尾;
                }

                // repeat の先頭 k-1 塩基は head の末尾のコピーなので、
                // そこまでしか踏み込まない窓は head の部分文字列そのもので、
                // head を読んだだけのリードでも必ず真になる (tail 側も同じ)
                var l_接合点1を跨ぐ = i < l_接合点1 && l_窓終端 > l_接合点1 + l_重なり長;
                var l_接合点2を跨ぐ = l_接合点2 < l_窓終端 && i < l_接合点2 - l_重なり長;
                if (!l_接合点1を跨ぐ && !l_接合点2を跨ぐ)
                {
                    continue;
                }

                if (l_直近の曖昧位置 < i && this.Has観測(l_テスト配列.AsSpan(i, this._r長)))
                {
                    l_支持数++;
                    l_入口支持数 += l_接合点1を跨ぐ ? 1 : 0;
                    l_出口支持数 += l_接合点2を跨ぐ ? 1 : 0;
                }
            }
            return (l_支持数, l_入口支持数, l_出口支持数);
        }

        /// <summary>
        /// r-mer の値を厳密な集合へ登録する
        /// </summary>
        /// <param name="p_値">正準化した r-mer の値</param>
        private void V_登録((UInt128 A_上位, UInt128 A_下位) p_値)
        {
            var l_分割 = Get_分割番号(p_値);
            lock (this._分割錠[l_分割])
            {
                _ = this._小集合?[l_分割].Add((ulong)p_値.A_下位) ?? this._中集合?[l_分割].Add(p_値.A_下位) ?? this._長集合![l_分割].Add(p_値);
            }
        }

        /// <summary>
        /// 128 塩基を超える r-mer を厳密な集合へ登録する
        /// </summary>
        /// <param name="p_キー">正準化した r-mer のキー</param>
        /// <remarks>
        /// キーは窓の作業領域を指すので、集合へ入れるときだけ複製する<br/>
        /// 先に所属を見るのは、リードの r-mer はほとんどが登録済みで、複製すると窓ごとに配列を捨てることになるため
        /// </remarks>
        private void V_登録(KmerKey p_キー)
        {
            var l_分割 = Get_分割番号(p_キー);
            var l_集合 = this._大集合![l_分割];
            lock (this._分割錠[l_分割])
            {
                if (!l_集合.Contains(p_キー))
                {
                    _ = l_集合.Add(p_キー.Get_複製());
                }
            }
        }

        /// <summary>
        /// r-mer をリードで見たか
        /// </summary>
        /// <param name="p_値">正準化した r-mer の値</param>
        /// <remarks>
        /// 構築を終えた後は読み取りだけなので錠を取らない
        /// </remarks>
        /// <returns>完全に一致する配列を見ていれば true</returns>
        private bool Has観測((UInt128 A_上位, UInt128 A_下位) p_値)
        {
            var l_分割 = Get_分割番号(p_値);
            return this._小集合?[l_分割].Contains((ulong)p_値.A_下位) ?? this._中集合?[l_分割].Contains(p_値.A_下位) ?? this._長集合![l_分割].Contains(p_値);
        }

        /// <summary>
        /// 窓の配列をリードで見たか
        /// </summary>
        /// <param name="p_窓">調べる配列 (長さは r 長、曖昧塩基を含まないこと)</param>
        /// <returns>完全に一致する配列を見ていれば true</returns>
        private bool Has観測(ReadOnlySpan<char> p_窓)
        {
            return this._大集合 is null ? this.Has観測(Get_正準値(p_窓)) : this.Has観測(new KmerKey(p_窓).Get_正規形());
        }

        /// <summary>
        /// 128 塩基を超える r-mer をリードで見たか
        /// </summary>
        /// <param name="p_キー">正準化した r-mer のキー</param>
        /// <returns>完全に一致する配列を見ていれば true</returns>
        private bool Has観測(KmerKey p_キー)
        {
            return this._大集合![Get_分割番号(p_キー)].Contains(p_キー);
        }

        /// <summary>
        /// r-mer の値から、登録先の分割を決める
        /// </summary>
        /// <param name="p_値">正準化した r-mer の値</param>
        /// <remarks>
        /// パック値の下位ビットは末尾の数塩基そのもので偏るため、混ぜてから上位ビットを使う
        /// </remarks>
        /// <returns></returns>
        private static int Get_分割番号((UInt128 A_上位, UInt128 A_下位) p_値)
        {
            var l_混合 = (ulong)p_値.A_下位 ^ (ulong)(p_値.A_下位 >> 64) ^ (ulong)p_値.A_上位;
            return (int)((l_混合 * 0x9E37_79B9_7F4A_7C15UL) >> (64 - 分割のビット数));
        }

        /// <summary>
        /// 128 塩基を超える r-mer のキーから、登録先の分割を決める
        /// </summary>
        /// <param name="p_キー">正準化した r-mer のキー</param>
        /// <returns></returns>
        private static int Get_分割番号(KmerKey p_キー)
        {
            return (int)(((ulong)p_キー.GetHashCode() * 0x9E37_79B9_7F4A_7C15UL) >> (64 - 分割のビット数));
        }

        /// <summary>
        /// 1 本のリードに現れる r-mer をすべて登録する
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_r長">r-mer の長さ</param>
        /// <param name="p_絞り込み">両端の k-mer がこの集合にある r-mer だけを登録する、null なら絞り込まない</param>
        /// <param name="p_k長">p_絞り込み の k</param>
        /// <remarks>
        /// 窓は 1 塩基ずつ転がして更新する (窓ごとに詰め直すとリード 1 本あたり O (n*r) になり、この呼び出しは全リード分繰り返される) <br/>
        /// 曖昧塩基が入ると窓は空に戻るので、それを跨ぐ窓は登録されない
        /// </remarks>
        private void V_登録_rMer(string p_リード, int p_r長, TrustedKmerIndex? p_絞り込み, int p_k長)
        {
            var l_k窓の信頼 = p_絞り込み is null ? default : p_リード.Length <= 1_024 ? stackalloc bool[p_リード.Length] : new bool[p_リード.Length];
            if (p_絞り込み is not null)
            {
                V_収集_k窓の信頼(p_リード, p_絞り込み, p_k長, l_k窓の信頼);
            }

            if (p_r長 <= 128)
            {
                var l_窓 = new RollingKmer(p_r長);
                for (var i = 0; i < p_リード.Length; i++)
                {
                    if (l_窓.Try追加(p_リード[i], out var l_キー)
                        && Is登録対象(l_k窓の信頼, p_絞り込み, i - p_r長 + 1, p_r長, p_k長)
                        && (this._候補集合?.Has観測(l_キー) ?? true))
                    {
                        this.V_登録(l_キー);
                    }
                }
                return;
            }

            var l_広い窓 = new WideRollingKmer(p_r長);
            for (var i = 0; i < p_リード.Length; i++)
            {
                if (l_広い窓.Try追加(p_リード[i], out var l_キー)
                    && Is登録対象(l_k窓の信頼, p_絞り込み, i - p_r長 + 1, p_r長, p_k長)
                    && (this._候補集合?.Has観測(l_キー) ?? true))
                {
                    this.V_登録(l_キー);
                }
            }
        }

        /// <summary>
        /// リード上の各位置から始まる k 窓が、信頼できる k-mer 集合にあるかを調べる
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_絞り込み">信頼できる k-mer 集合</param>
        /// <param name="p_k長">p_絞り込み の k</param>
        /// <param name="p_信頼">書き留め先</param>
        private static void V_収集_k窓の信頼(string p_リード, TrustedKmerIndex p_絞り込み, int p_k長, Span<bool> p_信頼)
        {
            if (p_k長 <= TrustedKmerIndex.パック値のk上限)
            {
                var l_k窓 = new RollingKmer(p_k長);
                for (var i = 0; i < p_リード.Length; i++)
                {
                    if (l_k窓.Try追加(p_リード[i], out var l_kキー))
                    {
                        p_信頼[i - p_k長 + 1] = p_絞り込み.Haskmer_正規形(l_kキー.A_上位, l_kキー.A_下位);
                    }
                }
                return;
            }

            var l_広いk窓 = new WideRollingKmer(p_k長);
            for (var i = 0; i < p_リード.Length; i++)
            {
                if (l_広いk窓.Try追加(p_リード[i], out var l_kキー))
                {
                    p_信頼[i - p_k長 + 1] = p_絞り込み.Haskmer_正規形(l_kキー);
                }
            }
        }

        /// <summary>
        /// その位置から始まる r 窓を登録するか
        /// </summary>
        /// <param name="p_k窓の信頼">リード上の各位置から始まる k 窓が信頼できる k-mer 集合にあるか</param>
        /// <param name="p_絞り込み">絞り込みに使う集合、null なら絞り込まない</param>
        /// <param name="p_開始">r 窓の開始位置</param>
        /// <param name="p_r長">r-mer の長さ</param>
        /// <param name="p_k長">p_絞り込み の k</param>
        /// <returns>登録するなら true</returns>
        private static bool Is登録対象(ReadOnlySpan<bool> p_k窓の信頼, TrustedKmerIndex? p_絞り込み, int p_開始, int p_r長, int p_k長)
        {
            return p_絞り込み is null || (p_k窓の信頼[p_開始] && p_k窓の信頼[p_開始 + p_r長 - p_k長]);
        }

        #endregion
    }
}
