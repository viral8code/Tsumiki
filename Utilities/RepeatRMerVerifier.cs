using Tsumiki.Commons;
using Tsumiki.IO;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 短い反復配列の経路付け替えに対する拒否権 (ABySS RResolver 型)
    /// </summary>
    /// <remarks>
    /// UnitigGraph.V_解決_短い反復 の対応付けは集計されたペア支持の多数決で決めており、
    /// その接続を実際に跨いだリードが存在するかを直接は確かめていない<br/>
    /// 生リードから作った r-mer の集合を使い、head → repeat → tail という経路の接合点を
    /// 実際に跨いだリードが存在するかを検定する<br/>
    /// 接合点を跨がない (どちらか一方の配列の内部だけに収まる) 窓は、
    /// その unitig 自身が既に実在の配列である以上、内部の k-mer は常に読まれているため、
    /// 経路の正しさに関わらず必ず真になり、判定対象から除く<br/>
    /// 跨いでいても踏み込みが k-1 塩基以内の窓は同じ理由で必ず真になる
    /// (repeat の先頭 k-1 塩基は head 末尾のコピーなので、その窓は head の部分文字列そのもの)<br/>
    /// 数えてよいのは共有区間を越えた窓だけで、その本数は接合点ごとに r - k 本になる<br/>
    /// head / repeat / tail は de Bruijn グラフの辺である以上、
    /// 隣接するもの同士は必ずアセンブリの k-1 塩基を共有しており、
    /// この共有区間は repeat 自身の配列の両端にもそのまま現れる<br/>
    /// r がこの k 以下だと接合点を跨ぐと判定した窓も共有区間の内側に収まってしまい、
    /// head だけ、repeat だけを独立に読んだリードでも真になって、常に素通りする拒否権になる<br/>
    /// 呼び出し側は r をこの反復解決で使う k より確実に長く取る責任を持つ
    /// (AssemblyPipeline は k + <see cref="Consts.rMer長のk超過分の既定値"/> で決めている)<br/>
    /// r-mer は存在確認だけに使うため、TrustedKmerIndex のような外部ソート付きカウンタは不要になる<br/>
    /// 32 塩基までは 2 bit パックが ulong に収まるので厳密な HashSet を使い、
    /// それを超える長さでは正規形のハッシュをブルームフィルタへ入れる<br/>
    /// 長い r-mer は誤り由来の種類数が膨らんで厳密な集合が数 GB になる一方、
    /// 偽陽性は棄却し損ねる方向にしか働かないため安全側に倒れる
    /// (調べないという以前の扱いは、偽陽性率 100 % と同じことだった)
    /// </remarks>
    internal sealed class RepeatRMerVerifier
    {
        #region 内部変数

        /// <summary>
        /// r-mer の厳密な集合、r が 32 を超える場合は null
        /// </summary>
        private readonly HashSet<ulong>? _rMer集合;

        /// <summary>
        /// r-mer のふるい、r が 32 以下の場合は null
        /// </summary>
        private readonly BloomFilter? _rMerふるい;

        /// <summary>
        /// r-mer の長さ
        /// </summary>
        private readonly int _r長;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// コンストラクタ
        /// </summary>
        /// <param name="p_rMer集合">r-mer の厳密な集合、r が 32 を超える場合は null</param>
        /// <param name="p_rMerふるい">r-mer のふるい、r が 32 以下の場合は null</param>
        /// <param name="p_r長">r-mer の長さ</param>
        private RepeatRMerVerifier(HashSet<ulong>? p_rMer集合, BloomFilter? p_rMerふるい, int p_r長)
        {
            this._rMer集合 = p_rMer集合;
            this._rMerふるい = p_rMerふるい;
            this._r長 = p_r長;
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

            var l_検証器 = p_r長 <= 32
                ? new RepeatRMerVerifier([], null, p_r長)
                : new RepeatRMerVerifier(null, Get_ふるい(l_パス群), p_r長);

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
        /// head → repeat → tail の経路上で、head-repeat 接合点と repeat-tail 接合点の
        /// 両方を実際に跨いだ r-mer のうち、生リード由来の集合に見つかった本数を返す
        /// </summary>
        /// <param name="p_head配列">接合点の手前の配列</param>
        /// <param name="p_repeat配列">反復配列</param>
        /// <param name="p_tail配列">接合点の先の配列</param>
        /// <returns>支持している r-mer の本数</returns>
        public int Get_接合点の支持数(string p_head配列, string p_repeat配列, string p_tail配列)
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
            for (var i = 0; i + this._r長 <= l_テスト配列.Length; i++)
            {
                var l_窓終端 = i + this._r長; // exclusive

                // repeat の先頭 k-1 塩基は head の末尾のコピーなので、
                // そこまでしか踏み込まない窓は head の部分文字列そのもので、
                // head を読んだだけのリードでも必ず真になる (tail 側も同じ)
                var l_接合点1を跨ぐ = i < l_接合点1 && l_窓終端 > l_接合点1 + l_重なり長;
                var l_接合点2を跨ぐ = l_接合点2 < l_窓終端 && i < l_接合点2 - l_重なり長;
                if (!l_接合点1を跨ぐ && !l_接合点2を跨ぐ)
                {
                    continue;
                }

                var l_曖昧か = false;
                for (var j = 0; j < this._r長; j++)
                {
                    if (Util.Get_曖昧塩基か(l_テスト配列[i + j]))
                    {
                        l_曖昧か = true;
                        break;
                    }
                }
                if (!l_曖昧か && this.Get_見たか(Get_正準値(l_テスト配列.AsSpan(i, this._r長))))
                {
                    l_支持数++;
                }
            }
            return l_支持数;
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
            return this.Get_接合点の支持数(p_head配列, p_repeat配列, p_tail配列) >= p_閾値;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// r-mer の値を集合、あるいはふるいへ登録する
        /// </summary>
        /// <param name="p_値">正準化した r-mer の値</param>
        private void V_登録(ulong p_値)
        {
            _ = (this._rMer集合?.Add(p_値));
            this._rMerふるい?.V_登録(p_値);
        }

        /// <summary>
        /// r-mer をリードで見たか
        /// </summary>
        /// <param name="p_値">正準化した r-mer の値</param>
        /// <returns>見ていれば true、ふるいを使う場合は偽陽性で true になりうる</returns>
        private bool Get_見たか(ulong p_値)
        {
            return this._rMer集合?.Contains(p_値) ?? this._rMerふるい!.Get_含まれるか(p_値);
        }

        /// <summary>
        /// ふるいの大きさをリードファイルの総量から決める
        /// </summary>
        /// <remarks>
        /// 実際の種類数は数え終わるまで分からないので、塩基数を上限とみなして
        /// 1 件あたり 4 ビットを見込み、確保量に上限を設ける
        /// </remarks>
        /// <param name="p_パス群">走査するリードファイルのパス</param>
        /// <returns>確保したふるい</returns>
        private static BloomFilter Get_ふるい(IEnumerable<string> p_パス群)
        {
            // FASTQ は塩基とクオリティで 1 塩基あたり約 2 バイトになる
            var l_見込み塩基数 = p_パス群.Sum(x => new FileInfo(x).Length) / 2;
            var l_ビット数 = Math.Clamp(4 * l_見込み塩基数, 1L << 20, Consts.rMerふるいのビット数上限);
            return new BloomFilter(l_ビット数, Consts.rMerふるいのハッシュ数);
        }

        /// <summary>
        /// 1 本のリードに現れる r-mer をすべて登録する
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_r長">r-mer の長さ</param>
        private void V_登録_rMer(string p_リード, int p_r長)
        {
            for (var i = 0; i + p_r長 <= p_リード.Length; i++)
            {
                var l_曖昧か = false;
                for (var j = 0; j < p_r長; j++)
                {
                    if (Util.Get_曖昧塩基か(p_リード[i + j]))
                    {
                        l_曖昧か = true;
                        break;
                    }
                }
                if (!l_曖昧か)
                {
                    this.V_登録(Get_正準値(p_リード.AsSpan(i, p_r長)));
                }
            }
        }

        /// <summary>
        /// 配列とその逆相補のうち、順鎖と逆鎖どちらから読んでも同一になるキーを返す
        /// </summary>
        /// <remarks>
        /// KmerKey は現在の実行時引数の k 長を前提にするため r-mer (k とは別の長さ) には使えず、
        /// KmerPacking の長さに依存しない API を使う
        /// </remarks>
        /// <param name="p_配列">元の配列</param>
        /// <returns>正準化したキー</returns>
        private static ulong Get_正準値(ReadOnlySpan<char> p_配列)
        {
            Span<byte> l_塩基ID列 = stackalloc byte[p_配列.Length];
            for (var i = 0; i < p_配列.Length; i++)
            {
                l_塩基ID列[i] = Util.Get_塩基ID(p_配列[i]);
            }
            return KmerPacking.Get_正規化ハッシュ_64(l_塩基ID列);
        }

        #endregion
    }
}
