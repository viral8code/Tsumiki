using Tsumiki.Commons;
using Tsumiki.IO;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 短い反復配列の経路付け替えに対する拒否権 (ABySS RResolver 型)
    /// </summary>
    /// <remarks>
    /// UnitigGraph.V_解決_短い反復 の対応付けは集計されたペア支持の多数決で決めており、その接続を実際に跨いだリードが存在するかを直接は確かめていない<br/>
    /// 生リードから作った r-mer の集合を使い、head → repeat → tail という経路の接合点を実際に跨いだリードが存在するかを検定する<br/>
    /// 接合点を跨がない (どちらか一方の配列の内部だけに収まる) 窓は、その unitig 自身が既に実在の配列である以上、内部の k-mer は常に読まれているため、経路の正しさに関わらず必ず真になり、判定対象から除く<br/>
    /// 跨いでいても踏み込みが k-1 塩基以内の窓は同じ理由で必ず真になる (repeat の先頭 k-1 塩基は head 末尾のコピーなので、その窓は head の部分文字列そのもの) <br/>
    /// 数えてよいのは共有区間を越えた窓だけで、その本数は接合点ごとに r - k 本になる<br/>
    /// head / repeat / tail は de Bruijn グラフの辺である以上、隣接するもの同士は必ずアセンブリの k-1 塩基を共有しており、この共有区間は repeat 自身の配列の両端にもそのまま現れる<br/>
    /// r がこの k 以下だと接合点を跨ぐと判定した窓も共有区間の内側に収まってしまい、head だけ、repeat だけを独立に読んだリードでも真になって、常に素通りする拒否権になる<br/>
    /// 呼び出し側は r をこの反復解決で使う k より確実に長く取る責任を持つ (AssemblyPipeline は k + <see cref="Consts.rMer長のk超過分の既定値"/> で決めている) <br/>
    /// r-mer は存在確認だけに使うため、TrustedKmerIndex のような外部ソート付きカウンタは不要になる<br/>
    /// 128 塩基までは 2 bit パックした厳密キーを使い、それより長い配列は正準文字列で照合する<br/>
    /// 偽陽性を持つ集合を経路支持に使わないため、長い r-mer の検査には追加メモリが必要になる
    /// </remarks>
    internal sealed class RepeatRMerVerifier
    {
        #region 内部変数

        /// <summary>
        /// 長さに依存せず完全一致を判定する正準キー集合
        /// </summary>
        private readonly HashSet<(UInt128 A_上位, UInt128 A_下位, string? A_長い配列)>? _rMer集合;

        /// <summary>
        /// 32 塩基までの厳密キー集合
        /// </summary>
        private readonly HashSet<ulong>? _小集合;

        /// <summary>
        /// 64 塩基までの厳密キー集合
        /// </summary>
        private readonly HashSet<UInt128>? _中集合;

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
            if (p_r長 <= 32)
            {
                this._小集合 = [];
            }
            else if (p_r長 <= 64)
            {
                this._中集合 = [];
            }
            else
            {
                this._rMer集合 = [];
            }
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 生リードファイル群を 1 回走査し、出現した r-mer (正準形) の集合を作る
        /// </summary>
        /// <remarks>
        /// 空、存在しないパスは片側リードのみの実行に対応するため無視する
        /// </remarks>
        /// <param name="p_リードパス一覧">走査するリードファイルのパス</param>
        /// <param name="p_r長">r-mer の長さ</param>
        /// <returns>構築した検証器</returns>
        public static RepeatRMerVerifier V_構築(IEnumerable<string> p_リードパス一覧, int p_r長)
        {
            if (p_r長 <= 0)
            {
                throw new ArgumentException("r-mer length must be positive");
            }

            var l_パス群 = p_リードパス一覧
                .Where(x => !string.IsNullOrWhiteSpace(x) && File.Exists(x))
                .ToList();

            var l_検証器 = new RepeatRMerVerifier(p_r長);

            foreach (var l_パス in l_パス群)
            {
                using var l_読み込み = new FastqReader(l_パス);
                while (l_読み込み.Get_続きがあるか())
                {
                    var l_リード = l_読み込み.Get_次のリード().A_生リード;
                    if (l_リード is not null)
                    {
                        l_検証器.V_登録_rMer(l_リード, p_r長);
                    }
                }
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
        public bool Get_接合点に支持があるか(string p_head配列, string p_repeat配列, string p_tail配列, int p_閾値)
        {
            var (A_全体, A_入口, A_出口) = this.Get_接合点別支持数(p_head配列, p_repeat配列, p_tail配列);
            return A_入口 > 0 && A_出口 > 0 && A_全体 >= p_閾値;
        }

        /// <summary>
        /// 配列とその逆相補のうち、順鎖と逆鎖どちらから読んでも同一になるキーを返す
        /// </summary>
        /// <remarks>
        /// KmerKey は現在の実行時引数の k 長を前提にするため r-mer (k とは別の長さ) には使えず、KmerPacking の長さに依存しない API を使う
        /// </remarks>
        /// <param name="p_配列">元の配列</param>
        /// <returns>正準化したキー</returns>
        public static (UInt128 A_上位, UInt128 A_下位, string? A_長い配列) Get_正準値(ReadOnlySpan<char> p_配列)
        {
            var l_逆鎖か = false;
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_順 = Util.Get_塩基ID(p_配列[i]);
                var l_逆 = 5 - Util.Get_塩基ID(p_配列[p_配列.Length - 1 - i]);
                if (l_順 != l_逆)
                {
                    l_逆鎖か = l_逆 < l_順;
                    break;
                }
            }
            var l_上位 = (UInt128)0;
            var l_下位 = (UInt128)0;
            var l_長い配列 = p_配列.Length > 128 ? new char[p_配列.Length] : null;
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_値 = l_逆鎖か ? 4 - Util.Get_塩基ID(p_配列[p_配列.Length - 1 - i]) : Util.Get_塩基ID(p_配列[i]) - 1;
                if (l_長い配列 is not null)
                {
                    l_長い配列[i] = "ACGT"[l_値];
                }
                else
                {
                    l_上位 = (l_上位 << 2) | (l_下位 >> 126);
                    l_下位 = (l_下位 << 2) | (uint)l_値;
                }
            }
            return (l_上位, l_下位, l_長い配列 is null ? null : new string(l_長い配列));
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
            // 窓ごとに Get_曖昧塩基か を r 回呼ぶと O (n*r) になるため、
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
                        if (Util.Get_曖昧塩基か(l_テスト配列[j]))
                        {
                            l_直近の曖昧位置 = j;
                        }
                    }
                }
                else if (Util.Get_曖昧塩基か(l_テスト配列[l_新規末尾]))
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

                if (l_直近の曖昧位置 < i && this.Get_見たか(Get_正準値(l_テスト配列.AsSpan(i, this._r長))))
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
        private void V_登録((UInt128 A_上位, UInt128 A_下位, string? A_長い配列) p_値)
        {
            _ = this._小集合?.Add((ulong)p_値.A_下位) ?? this._中集合?.Add(p_値.A_下位) ?? this._rMer集合!.Add(p_値);
        }

        /// <summary>
        /// r-mer をリードで見たか
        /// </summary>
        /// <param name="p_値">正準化した r-mer の値</param>
        /// <returns>完全に一致する配列を見ていれば true</returns>
        private bool Get_見たか((UInt128 A_上位, UInt128 A_下位, string? A_長い配列) p_値)
        {
            return this._小集合?.Contains((ulong)p_値.A_下位) ?? this._中集合?.Contains(p_値.A_下位) ?? this._rMer集合!.Contains(p_値);
        }

        /// <summary>
        /// 1 本のリードに現れる r-mer をすべて登録する
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_r長">r-mer の長さ</param>
        /// <remarks>
        /// 窓ごとに Get_曖昧塩基か を r 回呼ぶと全体で O (n*r) になる (この呼び出しは V_構築 から全リード分繰り返されるためリード規模でそのまま効く) <br/>
        /// 窓をスライドさせる際は新しく入る 1 塩基だけを見て「直近に見た曖昧塩基の位置」を更新すれば、その位置が現在の窓の左端以降にある間は判定を使い回せる (曖昧塩基は稀なので償却 O (n) で済む)
        /// </remarks>
        private void V_登録_rMer(string p_リード, int p_r長)
        {
            var l_直近の曖昧位置 = -1;
            for (var i = 0; i + p_r長 <= p_リード.Length; i++)
            {
                var l_新規末尾 = i + p_r長 - 1;
                if (i == 0)
                {
                    for (var j = 0; j < p_r長; j++)
                    {
                        if (Util.Get_曖昧塩基か(p_リード[j]))
                        {
                            l_直近の曖昧位置 = j;
                        }
                    }
                }
                else if (Util.Get_曖昧塩基か(p_リード[l_新規末尾]))
                {
                    l_直近の曖昧位置 = l_新規末尾;
                }

                if (l_直近の曖昧位置 >= i)
                {
                    continue;
                }

                this.V_登録(Get_正準値(p_リード.AsSpan(i, p_r長)));
            }
        }

        #endregion
    }
}
