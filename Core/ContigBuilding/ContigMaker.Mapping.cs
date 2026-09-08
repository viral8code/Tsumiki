using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker のうち、unitig への k-mer 索引構築とリードマッピングを担う部分。
    /// (contig 結合そのものは ContigMaker.cs、walk 構築は ContigMaker.Walk.cs、
    /// フラグメント長標本の収集は ContigMaker.FragmentSampling.cs を参照)
    /// </summary>
    internal partial class ContigMaker
    {
        // 同一 k-mer が複数の unitig にまたがって出現した(=反復配列等に
        // 由来する曖昧な k-mer である)ことを示す番兵値。
        // unitig ID は 1 始まりの正数、逆鎖側はその負数を使うため int.MinValue と衝突しない。
        private const int 曖昧kmerの番兵 = int.MinValue;

        // 値は (符号付きunitig ID, そのunitig内でのk-mer開始位置(0始まり、
        // 符号が示す向きの座標系))。位置情報は代表ユニティグの判定が
        // 「read内での最後のヒット位置」ではなく「unitig内での最後のヒット
        // 位置」を正しく求めるために必要(ギャップ長・インサートサイズ推定に使う)。
        private readonly Dictionary<KmerKey, (int A_ユニティグID, int A_開始位置)> _kmer辞書;

        // unitig ID(1始まり) -> unitig の塩基長。ギャップ長推定で
        // 「unitig の末尾からリードのヒット位置までの残り長」を求めるのに使う。
        private readonly Dictionary<int, int> _ユニティグ長;

        private readonly string _ユニティグファイルパス;

        // 単一リード内で直接検出された隣接(=k-1塩基のオーバーラップで
        // 実際に結合できる可能性が高い辺)。
        private readonly Dictionary<(int, int), ulong> _リード隣接;

        // ペアエンド情報(read1/read2 がそれぞれ別 unitig にマップされたこと)由来の
        // 隣接候補。キーは リード隣接 と同じ (始点, 終点) 形式(符号がunitigの向きを表す)。
        // 値は「観測されたペアの一覧」で、各観測ごとの既知長を保持し、
        // Scaffolder 側で代表値(中央値)を計算できるようにする。
        private readonly Dictionary<(int, int), List<int>> _ペア経路;

        // unitig ID(1始まり) -> その unitig が最終的にどの contig の
        // どの位置に配置されたか。contig 結合の実行後、Scaffolder から参照される。
        private readonly Dictionary<int, ユニティグ配置> _ユニティグ配置 = [];

        public ContigMaker(string p_ユニティグファイルパス)
        {
            this._ユニティグファイルパス = p_ユニティグファイルパス;
            this._kmer辞書 = [];
            this._ユニティグ長 = [];
            this._リード隣接 = [];
            this._ペア経路 = [];
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            using FastaReader l_読み込み = new(p_ユニティグファイルパス);
            var l_ID = 1;
            var l_短すぎるユニティグ数 = 0;
            var l_曖昧数 = 0;
            while (l_読み込み.Get_続きがあるか())
            {
                var l_ユニティグ = l_読み込み.Get_次の配列();
                this._ユニティグ長[l_ID] = l_ユニティグ.A_配列.Length;
                if (l_ユニティグ.A_配列.Length < l_k長)
                {
                    // k 未満の unitig は k-mer を持てずマッピング対象から漏れる。
                    // 黙って漏れないよう数だけ可視化しておく。
                    l_短すぎるユニティグ数++;
                    l_ID++;
                    continue;
                }
                for (var i = l_k長; i <= l_ユニティグ.A_配列.Length; i++)
                {
                    var l_開始位置 = i - l_k長;
                    var l_キー = new KmerKey(l_ユニティグ.A_配列.AsSpan(l_開始位置, l_k長));
                    var l_逆鎖キー = l_キー.Get_逆相補();
                    // 逆鎖キーは unitig 全体を逆相補した(=逆鎖の向きで読んだ)場合の
                    // 配列に対応する。区間 [開始位置, 開始位置+k長) を
                    // 長さ L の配列の逆側に写すと [L-i, L-開始位置) になるため、
                    // 逆鎖側での開始位置は L-i。
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

        /// <summary>
        /// k-mer辞書へ1件登録する。衝突した k-mer は後勝ちで上書きすると
        /// 別 unitig 由来のリードが同じ ID に見え、偽の隣接を作る。
        /// そのため曖昧としてマークし、マッピング時のヒットから除く。
        /// 戻り値は新たに曖昧マークを付けた件数(0 か 1)。
        /// </summary>
        private static int V_登録_kmer(
            Dictionary<KmerKey, (int, int)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
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
        /// unitig 間の隣接を de Bruijn グラフから厳密に構築する。
        /// V_結合_コンティグ を呼ぶ前(コピー数推定の接続伝播など)でも独立に
        /// 呼べるよう公開している。呼ぶたびに FASTA を読み直して新しい
        /// グラフを作る(unitig 数の規模では軽量なので使い捨てで構わない)。
        /// </summary>
        public UnitigGraph Get_グラフ()
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            List<string> l_ユニティグ配列 = [string.Empty, string.Empty];
            using (FastaReader l_読み込み = new(this._ユニティグファイルパス))
            {
                while (l_読み込み.Get_続きがあるか())
                {
                    var l_ユニティグ = l_読み込み.Get_次の配列().A_配列;
                    l_ユニティグ配列.Add(l_ユニティグ);
                    l_ユニティグ配列.Add(Util.V_逆相補(l_ユニティグ));
                }
            }

            return UnitigGraph.Get_グラフ(l_ユニティグ配列, this._kmer辞書, l_k長, 曖昧kmerの番兵);
        }

        public void V_マッピング_リード(string p_リードパス)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            // k-mer辞書 は構築後に変更されない読み取り専用データなので、
            // 複数スレッドから安全に参照できる。
            // 隣接への書き込みはスレッドごとにローカルな辞書に集計し、
            // 最後にマージすることでロックを避ける。
            var l_ローカル隣接 = new Dictionary<(int, int), ulong>[l_スレッド数];
            for (var i = 0; i < l_スレッド数; i++)
            {
                l_ローカル隣接[i] = [];
            }

            ReadPipeline.V_実行(
                l_スレッド数,
                l_スレッド数 * 256,
                Get_生リード列(p_リードパス),
                (l_リード, l_ワーカー番号) => this.V_マッピング_1リード(l_リード, l_ローカル隣接[l_ワーカー番号]));

            foreach (var l_ローカル in l_ローカル隣接)
            {
                foreach (var (l_キー, l_値) in l_ローカル)
                {
                    this._リード隣接[l_キー] =
                        this._リード隣接.TryGetValue(l_キー, out var l_既存) ? l_既存 + l_値 : l_値;
                }
            }
        }

        /// <summary>
        /// ペアエンドから unitig 間の隣接を検出する。
        /// 単一リードでは unitig 境界を跨げない場合でも、フラグメント長ぶん
        /// 離れた2つの unitig の隣接なら検出できる。
        /// </summary>
        public void V_マッピング_ペアリード(string p_リード1のパス, string p_リード2のパス)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            var l_ローカル隣接 = new Dictionary<(int, int), ulong>[l_スレッド数];
            // ローカルペア経路: (始点,終点) -> このワーカーで観測した各ペアの
            // 「既に見えている長さ」のリスト。
            var l_ローカルペア経路 = new Dictionary<(int, int), List<int>>[l_スレッド数];
            // ライブラリの向き(FR/RF/FF/RR)は決め打ちできないため、符号が
            // 一致するヒットと不一致のヒットを別々に集計し、多数派を採用する。
            var l_ローカル同一向き標本 = new List<int>[l_スレッド数];
            var l_ローカル逆向き標本 = new List<int>[l_スレッド数];
            for (var i = 0; i < l_スレッド数; i++)
            {
                l_ローカル隣接[i] = [];
                l_ローカルペア経路[i] = [];
                l_ローカル同一向き標本[i] = [];
                l_ローカル逆向き標本[i] = [];
            }

            ReadPipeline.V_実行(
                l_スレッド数,
                l_スレッド数 * 256,
                Get_ペアリード列(p_リード1のパス, p_リード2のパス),
                (l_ペア, l_ワーカー番号) => this.V_処理_1ペア(
                    l_ペア.A_リード1,
                    l_ペア.A_リード2,
                    l_ローカル隣接[l_ワーカー番号],
                    l_ローカルペア経路[l_ワーカー番号],
                    l_ローカル同一向き標本[l_ワーカー番号],
                    l_ローカル逆向き標本[l_ワーカー番号]));

            foreach (var l_ローカル in l_ローカル隣接)
            {
                foreach (var (l_キー, l_値) in l_ローカル)
                {
                    this._リード隣接[l_キー] =
                        this._リード隣接.TryGetValue(l_キー, out var l_既存) ? l_既存 + l_値 : l_値;
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
            // 多数派の側だけを実際のライブラリ配置として採用する。
            // 少数派側は測定ノイズ・誤マッピング・稀な異常配置とみなして捨てる。
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
                // 同一unitig内標本は、unitig自体がフラグメント長より短い場合
                // 両端が同じunitig内に収まるペアしか観測できず、より短い
                // フラグメントに偏った標本になりやすい(unitigが短いほど顕著)。
                Logger.V_出力(メッセージID.同一ユニティグの断片長分布, Get_分布要約(l_同一ユニティグ標本));
                Logger.V_出力(メッセージID.同一ユニティグの断片長中央値, StatsUtil.Get_中央値(l_同一ユニティグ標本), l_同一ユニティグ標本.Count);
            }
        }

        /// <summary>
        /// 1本のリードが代表としてどの unitig にマップされるかを判定する。
        /// 最多得票の unitig ID と、ギャップ長推定に使う最終ヒット位置を返す。
        /// </summary>
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
            // 記録するのは unitig 内での終端位置であって read 内での位置ではない。
            // 両者は unitig が read より十分長いと大きく食い違う。
            var l_最終終端位置 = new Dictionary<int, int>();
            var l_曖昧塩基数 = 0;
            for (var i = 0; i < l_k長; i++)
            {
                if (Util.Get_曖昧塩基か(p_リード[i]))
                {
                    l_曖昧塩基数++;
                }
            }
            for (var i = l_k長; i <= p_リード.Length; i++)
            {
                if (Util.Get_曖昧塩基か(p_リード[i - l_k長]))
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
        /// read1/read2 を同時に読み進めて対応するペアを返す。
        /// ID の対応が取れないものと片側だけ残ったものは A_リード2 を空文字にし、
        /// 単一リード内の隣接検出だけは通常どおり行えるようにする。
        /// </summary>
        private static IEnumerable<(string A_リード1, string A_リード2)> Get_ペアリード列(
            string p_リード1のパス, string p_リード2のパス)
        {
            using var l_読み込み1 = new FastqReader(p_リード1のパス);
            using var l_読み込み2 = new FastqReader(p_リード2のパス);

            var l_不一致を警告済みか = false;
            while (l_読み込み1.Get_続きがあるか() && l_読み込み2.Get_続きがあるか())
            {
                var l_データ1 = l_読み込み1.Get_次のリード();
                var l_データ2 = l_読み込み2.Get_次のリード();

                if (Util.Get_ペア共通ID(l_データ1.A_ID) != Util.Get_ペア共通ID(l_データ2.A_ID))
                {
                    if (!l_不一致を警告済みか)
                    {
                        Logger.V_出力(メッセージID.ペアリードIDの不一致, l_データ1.A_ID, l_データ2.A_ID);
                        l_不一致を警告済みか = true;
                    }
                    // お互いを誤ってペアとして扱わないよう、別々に流す。
                    yield return (l_データ1.A_生リード!, string.Empty);
                    yield return (l_データ2.A_生リード!, string.Empty);
                    continue;
                }

                yield return (l_データ1.A_生リード!, l_データ2.A_生リード!);
            }

            // 片方のファイルだけ残っている場合は単一リードとして処理する。
            while (l_読み込み1.Get_続きがあるか())
            {
                yield return (l_読み込み1.Get_次のリード().A_生リード!, string.Empty);
            }
            while (l_読み込み2.Get_続きがあるか())
            {
                yield return (l_読み込み2.Get_次のリード().A_生リード!, string.Empty);
            }
        }

        /// <summary>
        /// ペア1組を処理する。ペアエンド由来の隣接は直接のオーバーラップを
        /// 保証しない弱い証拠なので、リード隣接とは分けて集計する。
        /// </summary>
        private void V_処理_1ペア(
            string p_リード1,
            string p_リード2,
            Dictionary<(int, int), ulong> p_ローカル隣接,
            Dictionary<(int, int), List<int>> p_ローカルペア経路,
            List<int> p_同一向き標本,
            List<int> p_逆向き標本)
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
                V_収集_同一ユニティグ標本(
                    l_ヒット1, l_ヒット2, p_リード1, p_リード2, p_同一向き標本, p_逆向き標本);
            }
            else
            {
                V_収集_ペア経路(l_ヒット1, l_ヒット2, p_リード1, p_リード2, p_ローカルペア経路);
            }
        }

        /// <summary>FASTQ を順に読み進めて生リード文字列だけを返す。</summary>
        private static IEnumerable<string> Get_生リード列(string p_リードパス)
        {
            using var l_読み込み = new FastqReader(p_リードパス);
            while (l_読み込み.Get_続きがあるか())
            {
                yield return l_読み込み.Get_次のリード().A_生リード!;
            }
        }

        private void V_マッピング_1リード(string p_リード, Dictionary<(int, int), ulong> p_ローカル隣接)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;

            // k 未満のリードからは k-mer を取れない。この判定が無いと
            // 下の初期化ループが p_リード[i] を i = k-1 まで舐めて範囲外になる。
            if (p_リード.Length < l_k長)
            {
                return;
            }

            // FASTQ の生リードには N 等の曖昧塩基が混入しうるため、
            // A/C/G/T のみを前提とする厳密版ではなく曖昧塩基を許容する版を使う。
            // 曖昧塩基を含む区間の k-mer は後段のカウントによるスキップで除外される。
            var l_逆鎖リード = Util.V_逆相補_曖昧塩基あり(p_リード);
            var l_直前 = 0;
            var l_逆鎖の直前 = 0;
            var l_曖昧塩基数 = 0;
            var l_逆鎖の曖昧塩基数 = 0;
            for (var i = 0; i < l_k長; i++)
            {
                if (Util.Get_曖昧塩基か(p_リード[i]))
                {
                    l_曖昧塩基数++;
                }
                if (Util.Get_曖昧塩基か(l_逆鎖リード[i]))
                {
                    l_逆鎖の曖昧塩基数++;
                }
            }
            for (var i = l_k長; i <= p_リード.Length; i++)
            {
                if (Util.Get_曖昧塩基か(p_リード[i - l_k長]))
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
                            p_ローカル隣接[l_経路キー] =
                                p_ローカル隣接.TryGetValue(l_経路キー, out var l_件数) ? l_件数 + 1 : 1;
                            // 直前にヒットした unitig を更新する。これを怠ると、
                            // リード内で3つ以上の unitig にまたがった場合でも
                            // 常に「最初にヒットした unitig」との組しか記録されず、
                            // 実際の隣接関係(直前→直後)を反映できない。
                            l_直前 = l_ID;
                        }
                    }
                }
                if (Util.Get_曖昧塩基か(l_逆鎖リード[i - l_k長]))
                {
                    l_逆鎖の曖昧塩基数--;
                }
                if (l_逆鎖の曖昧塩基数 == 0)
                {
                    var l_逆鎖キー = new KmerKey(l_逆鎖リード.AsSpan(i - l_k長, l_k長));
                    if (this._kmer辞書.TryGetValue(l_逆鎖キー, out var l_逆鎖項目) && l_逆鎖項目.A_ユニティグID != 曖昧kmerの番兵)
                    {
                        var l_逆鎖ID = l_逆鎖項目.A_ユニティグID;
                        if (l_逆鎖の直前 == 0)
                        {
                            l_逆鎖の直前 = l_逆鎖ID;
                        }
                        else if (l_逆鎖の直前 != l_逆鎖ID)
                        {
                            var l_経路キー = (l_逆鎖の直前, l_逆鎖ID);
                            p_ローカル隣接[l_経路キー] =
                                p_ローカル隣接.TryGetValue(l_経路キー, out var l_件数) ? l_件数 + 1 : 1;
                            l_逆鎖の直前 = l_逆鎖ID;
                        }
                    }
                }
            }
        }
    }
}
