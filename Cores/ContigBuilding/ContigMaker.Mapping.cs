using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.IO;
using Tsumiki.Models.ContigBuilding;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.UnitigBuilding;
using Tsumiki.Utilities;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker のうち、unitig への k-mer 索引構築とリードマッピングを担う部分
    /// </summary>
    /// <remarks>
    /// (contig 結合そのものは ContigMaker.cs、walk 構築は ContigMaker.Walk.cs、フラグメント長標本の収集は ContigMaker.FragmentSampling.cs を参照)
    /// </remarks>
    internal partial class ContigMaker
    {
        #region 定数

        /// <summary>
        /// 曖昧 kmer の番兵
        /// </summary>
        private const int 曖昧kmerの番兵 = int.MinValue;

        #endregion

        #region 内部変数

        /// <summary>
        /// k-mer から、それが載る unitig と開始位置を引く辞書
        /// </summary>
        private readonly Dictionary<KmerKey, (int A_ユニティグID, int A_開始位置)> _kmer辞書;

        /// <summary>
        /// unitig 長
        /// </summary>
        private readonly Dictionary<int, int> _ユニティグ長;

        /// <summary>
        /// 正逆両鎖の unitig 配列
        /// </summary>
        private readonly List<string> _ユニティグ配列;

        /// <summary>
        /// リードが跨いだ unitig の組と、その本数
        /// </summary>
        private readonly Dictionary<(int, int), ulong> _リード隣接;

        /// <summary>
        /// ペアが跨いだ unitig の組と、その間に通った頂点
        /// </summary>
        private readonly Dictionary<(int, int), List<int>> _ペア経路;

        /// <summary>
        /// unitig 配置
        /// </summary>
        private readonly Dictionary<int, ユニティグ配置> _ユニティグ配置 = [];

        #endregion

        #region コンストラクタ

        /// <summary>
        /// コンストラクタ
        /// </summary>
        /// <param name="p_ユニティグファイルパス"></param>
        public ContigMaker(string p_ユニティグファイルパス)
        {
            this._kmer辞書 = [];
            this._ユニティグ長 = [];
            this._ユニティグ配列 = [string.Empty, string.Empty];
            this._リード隣接 = [];
            this._ペア経路 = [];
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            using FastaReader l_読み込み = new(p_ユニティグファイルパス);
            var l_ID = 1;
            var l_短すぎるユニティグ数 = 0;
            var l_曖昧数 = 0;
            while (l_読み込み.Has続き())
            {
                var l_ユニティグ = l_読み込み.Get_次の配列();
                this._ユニティグ長[l_ID] = l_ユニティグ.A_配列.Length;
                this._ユニティグ配列.Add(l_ユニティグ.A_配列);
                this._ユニティグ配列.Add(Util.V_逆相補(l_ユニティグ.A_配列));

                if (l_ユニティグ.A_配列.Length < l_k長)
                {
                    // k 未満の unitig は k-mer を持てずマッピング対象から漏れる
                    // 黙って漏れないよう数だけ可視化しておく
                    l_短すぎるユニティグ数++;
                    l_ID++;
                    continue;
                }
                for (var i = l_k長; i <= l_ユニティグ.A_配列.Length; i++)
                {
                    var l_開始位置 = i - l_k長;
                    var l_キー = new KmerKey(l_ユニティグ.A_配列.AsSpan(l_開始位置, l_k長));
                    var l_逆鎖キー = l_キー.Get_逆相補();

                    // 逆鎖キーは unitig 全体を逆相補した (=逆鎖の向きで読んだ) 場合の
                    // 配列に対応する
                    // 区間 [開始位置, 開始位置+k 長) を
                    // 長さ L の配列の逆側に写すと [L-i, L-開始位置) になるため、
                    // 逆鎖側での開始位置は L-i
                    var l_逆鎖開始位置 = l_ユニティグ.A_配列.Length - i;
                    l_曖昧数 += V_登録_kmer(this._kmer辞書, l_キー, l_ID, l_開始位置);
                    l_曖昧数 += V_登録_kmer(this._kmer辞書, l_逆鎖キー, -l_ID, l_逆鎖開始位置);
                }
                l_ID++;
            }

            if (l_短すぎるユニティグ数 > 0)
            {
                Logger.V_出力(メッセージID.短すぎるユニティグの除外, l_短すぎるユニティグ数);
            }

            if (l_曖昧数 > 0)
            {
                Logger.V_出力(メッセージID.曖昧なkmer登録, l_曖昧数);
            }
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// unitig 間の隣接を de Bruijn グラフから厳密に構築する
        /// </summary>
        /// <remarks>
        /// <see cref="V_結合_コンティグ"/> を呼ぶ前 (コピー数推定の接続伝播など) でも独立に呼べるよう公開している<br/>
        /// 呼ぶたびに FASTA を読み直して新しいグラフを作る (unitig 数の規模では軽量なので使い捨てで構わない)
        /// </remarks>
        /// <returns></returns>
        public UnitigGraph Get_グラフ()
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            return UnitigGraph.Get_グラフ(this._ユニティグ配列, this._kmer辞書, l_k長, 曖昧kmerの番兵);
        }

        /// <summary>
        /// リードを unitig へ貼り付け、隣接とペアの支持を集める
        /// </summary>
        /// <param name="p_リードパス">貼り付けるリードのパス</param>
        public void V_マッピング_リード(string p_リードパス)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            // k-mer 辞書 は構築後に変更されない読み取り専用データなので、
            // 複数スレッドから安全に参照できる
            // 隣接への書き込みはスレッドごとにローカルな辞書に集計し、
            // 最後にマージすることでロックを避ける
            var l_ローカル隣接 = new Dictionary<(int, int), ulong>[l_スレッド数];
            for (var i = 0; i < l_スレッド数; i++)
            {
                l_ローカル隣接[i] = [];
            }

            ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, Get_生リード列(p_リードパス), (l_リード, l_ワーカー番号) => this.V_マッピング_1リード(l_リード, l_ローカル隣接[l_ワーカー番号]));

            foreach (var l_ローカル in l_ローカル隣接)
            {
                foreach (var (l_キー, l_値) in l_ローカル)
                {
                    this._リード隣接[l_キー] = this._リード隣接.TryGetValue(l_キー, out var l_既存) ? l_既存 + l_値 : l_値;
                }
            }
        }

        /// <summary>
        /// ペアエンドから unitig 間の隣接を検出する
        /// </summary>
        /// <param name="p_リード1のパス"></param>
        /// <param name="p_リード2のパス"></param>
        /// <remarks>
        /// 単一リードでは unitig 境界を跨げない場合でも、フラグメント長ぶん離れた 2 つの unitig の隣接なら検出できる
        /// </remarks>
        public void V_マッピング_ペアリード(string p_リード1のパス, string p_リード2のパス)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            var l_ローカル隣接 = new Dictionary<(int, int), ulong>[l_スレッド数];

            // ローカルペア経路: (始点,終点) -> このワーカーで観測した各ペアの
            // 「既に見えている長さ」のリスト
            var l_ローカルペア経路 = new Dictionary<(int, int), List<int>>[l_スレッド数];

            // ライブラリの向き (FR/RF/FF/RR) は決め打ちできないため、符号が
            // 一致するヒットと不一致のヒットを別々に集計し、多数派を採用する
            var l_ローカル同一向き標本 = new List<int>[l_スレッド数];
            var l_ローカル逆向き標本 = new List<int>[l_スレッド数];
            for (var i = 0; i < l_スレッド数; i++)
            {
                l_ローカル隣接[i] = [];
                l_ローカルペア経路[i] = [];
                l_ローカル同一向き標本[i] = [];
                l_ローカル逆向き標本[i] = [];
            }

            ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, Get_ペアリード列(p_リード1のパス, p_リード2のパス), (l_ペア, l_ワーカー番号) => this.V_処理_1ペア(l_ペア.A_リード1, l_ペア.A_リード2, l_ローカル隣接[l_ワーカー番号], l_ローカルペア経路[l_ワーカー番号], l_ローカル同一向き標本[l_ワーカー番号], l_ローカル逆向き標本[l_ワーカー番号]));

            foreach (var l_ローカル in l_ローカル隣接)
            {
                foreach (var (l_キー, l_値) in l_ローカル)
                {
                    this._リード隣接[l_キー] = this._リード隣接.TryGetValue(l_キー, out var l_既存) ? l_既存 + l_値 : l_値;
                }
            }

            foreach (var l_ローカルペア in l_ローカルペア経路)
            {
                foreach (var (l_キー, l_値群) in l_ローカルペア)
                {
                    if (this._ペア経路.TryGetValue(l_キー, out var l_一覧))
                    {
                        l_一覧.AddRange(l_値群);
                    }
                    else
                    {
                        this._ペア経路[l_キー] = [.. l_値群];
                    }
                }
            }

            // 「符号一致」「符号不一致」それぞれの総標本数を集計し、
            // 多数派の側だけを実際のライブラリ配置として採用する
            // 少数派側は測定ノイズ・誤マッピング・稀な異常配置とみなして捨てる
            var l_同一向き合計 = l_ローカル同一向き標本.Sum(x => x.Count);
            var l_逆向き合計 = l_ローカル逆向き標本.Sum(x => x.Count);

            IEnumerable<List<int>> l_採用する標本群;
            string l_採用ラベル;

            if (l_同一向き合計 == 0 && l_逆向き合計 == 0)
            {
                l_採用する標本群 = [];
                l_採用ラベル = "none";
            }
            else if (l_同一向き合計 >= l_逆向き合計)
            {
                l_採用する標本群 = l_ローカル同一向き標本;
                l_採用ラベル = "same-orientation";
            }
            else
            {
                l_採用する標本群 = l_ローカル逆向き標本;
                l_採用ラベル = "opposite-orientation";
            }

            var l_同一ユニティグ標本 = new List<int>();
            foreach (var l_標本 in l_採用する標本群)
            {
                l_同一ユニティグ標本.AddRange(l_標本);
            }
            this.A_インサートサイズ標本.AddRange(l_同一ユニティグ標本);
            this.A_同一ユニティグ標本.AddRange(l_同一ユニティグ標本);

            var l_ペア支持数 = this._ペア経路.Values.Sum(x => x.Count);
            Logger.V_出力(メッセージID.ペア隣接候補数, this._ペア経路.Count, l_ペア支持数);
            Logger.V_出力(メッセージID.同一ユニティグのペア向き集計, l_同一向き合計, l_逆向き合計, l_採用ラベル, l_同一ユニティグ標本.Count);
            if (l_同一ユニティグ標本.Count > 0)
            {
                // 同一 unitig 内標本は、unitig 自体がフラグメント長より短い場合
                // 両端が同じ unitig 内に収まるペアしか観測できず、より短い
                // フラグメントに偏った標本になりやすい (unitig が短いほど顕著)
                Logger.V_出力(メッセージID.同一ユニティグの断片長分布, Get_分布要約(l_同一ユニティグ標本));
                Logger.V_出力(メッセージID.同一ユニティグの断片長中央値, StatsUtil.Get_中央値(l_同一ユニティグ標本), l_同一ユニティグ標本.Count);
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k-mer 辞書へ 1 件登録する
        /// </summary>
        /// <param name="p_辞書"></param>
        /// <param name="p_キー"></param>
        /// <param name="p_ID"></param>
        /// <param name="p_位置"></param>
        /// <remarks>
        /// 衝突した k-mer は後勝ちで上書きすると別 unitig 由来のリードが同じ ID に見え、偽の隣接を作る<br/>
        /// そのため曖昧としてマークし、マッピング時のヒットから除く<br/>
        /// 戻り値は新たに曖昧マークを付けた件数 (0 か 1)
        /// </remarks>
        /// <returns></returns>
        private static int V_登録_kmer(Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存))
            {
                if (l_既存.Item1 == 曖昧kmerの番兵 || l_既存.Item1 == p_ID)
                {
                    return 0;
                }
                p_辞書[p_キー] = (曖昧kmerの番兵, 0);
                return 1;
            }
            p_辞書[p_キー] = (p_ID, p_位置);
            return 0;
        }

        /// <summary>
        /// 1 本のリードが代表としてどの unitig にマップされるかを判定する
        /// </summary>
        /// <remarks>
        /// 最多得票の unitig ID と、ギャップ長推定に使う最終ヒット位置を返す
        /// </remarks>
        /// <param name="p_リード"></param>
        /// <returns></returns>
        internal 代表ユニティグヒット Get_代表ユニティグ(string p_リード)
        {
            if (string.IsNullOrEmpty(p_リード))
            {
                return 代表ユニティグヒット.A_ヒットなし;
            }

            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            if (p_リード.Length < l_k長)
            {
                return 代表ユニティグヒット.A_ヒットなし;
            }

            var l_得票 = new Dictionary<int, int>();

            // 記録するのは unitig 内での終端位置であって read 内での位置ではない
            // 両者は unitig が read より十分長いと大きく食い違う
            var l_最終終端位置 = new Dictionary<int, int>();
            var l_曖昧塩基数 = 0;
            for (var i = 0; i < l_k長; i++)
            {
                if (Util.Is曖昧塩基(p_リード[i]))
                {
                    l_曖昧塩基数++;
                }
            }
            for (var i = l_k長; i <= p_リード.Length; i++)
            {
                if (Util.Is曖昧塩基(p_リード[i - l_k長]))
                {
                    l_曖昧塩基数--;
                }

                if (l_曖昧塩基数 == 0)
                {
                    var l_キー = new KmerKey(p_リード.AsSpan(i - l_k長, l_k長));
                    if (this._kmer辞書.TryGetValue(l_キー, out var l_項目) && l_項目.A_ユニティグID != 曖昧kmerの番兵)
                    {
                        var l_ID = l_項目.A_ユニティグID;
                        l_得票[l_ID] = l_得票.GetValueOrDefault(l_ID) + 1;
                        l_最終終端位置[l_ID] = l_項目.A_開始位置 + l_k長;
                    }
                }
            }

            if (l_得票.Count == 0)
            {
                return 代表ユニティグヒット.A_ヒットなし;
            }

            var l_最良 = 0;
            var l_最多得票 = 0;
            foreach (var (l_ID, l_票数) in l_得票)
            {
                if (l_票数 > l_最多得票)
                {
                    l_最良 = l_ID;
                    l_最多得票 = l_票数;
                }
            }

            var l_ユニティグ長 = this._ユニティグ長.GetValueOrDefault(Math.Abs(l_最良), 0);
            return new 代表ユニティグヒット(l_最良, l_最多得票, l_最終終端位置[l_最良], l_ユニティグ長);
        }

        /// <summary>
        /// read1/read2 を同時に読み進めて対応するペアを返す
        /// </summary>
        /// <param name="p_リード1のパス"></param>
        /// <param name="p_リード2のパス"></param>
        /// <remarks>
        /// ID の対応が取れないものと片側だけ残ったものは A_リード2 を空文字にし、単一リード内の隣接検出だけは通常どおり行えるようにする
        /// </remarks>
        /// <returns></returns>
        private static IEnumerable<(string A_リード1, string A_リード2)> Get_ペアリード列(string p_リード1のパス, string p_リード2のパス)
        {
            using var l_読み込み1 = new FastqReader(p_リード1のパス);
            using var l_読み込み2 = new FastqReader(p_リード2のパス);

            var l_Is不一致警告済み = false;
            while (l_読み込み1.Has続き() && l_読み込み2.Has続き())
            {
                var l_データ1 = l_読み込み1.Get_次のリード();
                var l_データ2 = l_読み込み2.Get_次のリード();

                if (Util.Get_ペア共通ID(l_データ1.A_ID) != Util.Get_ペア共通ID(l_データ2.A_ID))
                {
                    if (!l_Is不一致警告済み)
                    {
                        Logger.V_出力(メッセージID.ペアリードIDの不一致, l_データ1.A_ID, l_データ2.A_ID);
                        l_Is不一致警告済み = true;
                    }

                    // お互いを誤ってペアとして扱わないよう、別々に流す
                    yield return (l_データ1.A_生リード!, string.Empty);
                    yield return (l_データ2.A_生リード!, string.Empty);
                    continue;
                }

                yield return (l_データ1.A_生リード!, l_データ2.A_生リード!);
            }

            // 片方のファイルだけ残っている場合は単一リードとして処理する
            while (l_読み込み1.Has続き())
            {
                yield return (l_読み込み1.Get_次のリード().A_生リード!, string.Empty);
            }
            while (l_読み込み2.Has続き())
            {
                yield return (l_読み込み2.Get_次のリード().A_生リード!, string.Empty);
            }
        }

        /// <summary>
        /// ペア 1 組を処理する
        /// </summary>
        /// <param name="p_リード1"></param>
        /// <param name="p_リード2"></param>
        /// <param name="p_ローカル隣接"></param>
        /// <param name="p_ローカルペア経路"></param>
        /// <param name="p_同一向き標本"></param>
        /// <param name="p_逆向き標本"></param>
        /// <remarks>
        /// ペアエンド由来の隣接は直接のオーバーラップを保証しない弱い証拠なので、リード隣接とは分けて集計する
        /// </remarks>
        private void V_処理_1ペア(string p_リード1, string p_リード2, Dictionary<(int, int), ulong> p_ローカル隣接, Dictionary<(int, int), List<int>> p_ローカルペア経路, List<int> p_同一向き標本, List<int> p_逆向き標本)
        {
            this.V_マッピング_1リード(p_リード1, p_ローカル隣接);
            this.V_マッピング_1リード(p_リード2, p_ローカル隣接);

            var l_ヒット1 = this.Get_代表ユニティグ(p_リード1);
            var l_ヒット2 = this.Get_代表ユニティグ(p_リード2);

            if (l_ヒット1.A_ユニティグID == 0 || l_ヒット2.A_ユニティグID == 0)
            {
                return;
            }

            if (Math.Abs(l_ヒット1.A_ユニティグID) == Math.Abs(l_ヒット2.A_ユニティグID))
            {
                V_収集_同一ユニティグ標本(l_ヒット1, l_ヒット2, p_リード1, p_リード2, p_同一向き標本, p_逆向き標本);
            }
            else
            {
                V_収集_ペア経路(l_ヒット1, l_ヒット2, p_リード1, p_リード2, p_ローカルペア経路);
            }
        }

        /// <summary>
        /// FASTQ を順に読み進めて生リード文字列だけを返す
        /// </summary>
        /// <param name="p_リードパス"></param>
        /// <returns></returns>
        private static IEnumerable<string> Get_生リード列(string p_リードパス)
        {
            using var l_読み込み = new FastqReader(p_リードパス);
            while (l_読み込み.Has続き())
            {
                yield return l_読み込み.Get_次のリード().A_生リード!;
            }
        }

        /// <summary>
        /// 1 本のリードを貼り付ける
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_ローカル隣接">このワーカーが集めた隣接</param>
        private void V_マッピング_1リード(string p_リード, Dictionary<(int, int), ulong> p_ローカル隣接)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            // k 未満のリードからは k-mer を取れない
            // この判定が無いと
            // 下の初期化ループが p_リード[i] を i = k-1 まで舐めて範囲外になる
            if (p_リード.Length < l_k長)
            {
                return;
            }

            var l_直前 = 0;
            var l_曖昧塩基数 = 0;

            // 索引には両鎖があり、辺重みは後段で逆鎖対称にするため入力方向だけを走査する
            for (var i = 0; i < l_k長; i++)
            {
                if (Util.Is曖昧塩基(p_リード[i]))
                {
                    l_曖昧塩基数++;
                }
            }
            for (var i = l_k長; i <= p_リード.Length; i++)
            {
                if (Util.Is曖昧塩基(p_リード[i - l_k長]))
                {
                    l_曖昧塩基数--;
                }

                if (l_曖昧塩基数 == 0)
                {
                    var l_キー = new KmerKey(p_リード.AsSpan(i - l_k長, l_k長));
                    if (this._kmer辞書.TryGetValue(l_キー, out var l_項目) && l_項目.A_ユニティグID != 曖昧kmerの番兵)
                    {
                        var l_ID = l_項目.A_ユニティグID;
                        if (l_直前 == 0)
                        {
                            l_直前 = l_ID;
                        }
                        else if (l_直前 != l_ID)
                        {
                            var l_経路キー = (l_直前, l_ID);
                            p_ローカル隣接[l_経路キー] = p_ローカル隣接.TryGetValue(l_経路キー, out var l_件数) ? l_件数 + 1UL : 1UL;

                            // 直前にヒットした unitig を更新する
                            // これを怠ると、
                            // リード内で 3 つ以上の unitig にまたがった場合でも
                            // 常に「最初にヒットした unitig」との組しか記録されず、
                            // 実際の隣接関係 (直前→直後) を反映できない
                            l_直前 = l_ID;
                        }
                    }
                }

            }
        }

        #endregion
    }
}
