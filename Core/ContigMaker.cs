using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    internal class ContigMaker
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
                    // k-mer 長より短い unitig は本来 k-mer を1つも持てないため、
                    // 辞書に登録できずリードマッピングの対象から漏れる。
                    // 完全な解決にはより短い k-mer での再マッピング等が必要だが、
                    // ここでは少なくとも「登録が1件もされないまま黙って ID が進む」
                    // 状態を可視化するためカウントしておく。
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
                Console.WriteLine($"[Warning] {l_短すぎるユニティグ数} unitig(s) shorter than k-mer length were skipped in mapping.");
            }
            if (l_曖昧数 > 0)
            {
                Console.WriteLine($"[Warning] {l_曖昧数} k-mer registration(s) were ambiguous (shared by multiple unitigs) and will be ignored during mapping.");
            }
        }

        /// <summary>
        /// k-mer辞書へ1件登録する。既に別の unitig ID が登録されていた場合、
        /// そのままでは後勝ちで上書きされてしまい、実際には異なる unitig 由来の
        /// リードが同じ ID にマップされたかのように誤って隣接関係を作ってしまう。
        /// これを防ぐため、衝突した k-mer は曖昧としてマークし、
        /// マッピング時にはヒットとして扱わないようにする。
        /// 戻り値: 新たに曖昧マークを付けた場合は 1、そうでなければ 0。
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
        /// ペアエンドリードの情報を使って unitig 間の隣接関係を検出する。
        /// 単一リード内で複数 unitig をまたぐ場合に加え、
        /// read1 と read2 がそれぞれ別々の unitig に(単独で)マップされた場合も、
        /// 「インサートサイズ程度の距離で隣接している」という情報として記録する。
        /// これにより、unitig 長がリード長よりずっと長く単一リードでは境界を
        /// またげないケースでも隣接関係を検出できる。
        ///
        /// read1/read2 が本当にペアであることを、リード ID の対応(/1,/2 や
        /// Casava 1.8+ の記法)で検証する。対応が取れない場合は警告を出し、
        /// そのペアはペアエンド由来の隣接検出をスキップする
        /// (単一リード内の隣接検出は通常どおり行う)。
        /// </summary>
        public void V_マッピング_ペアリード(string p_リード1のパス, string p_リード2のパス)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            var l_ローカル隣接 = new Dictionary<(int, int), ulong>[l_スレッド数];
            // ローカルペア経路: (始点,終点) -> このワーカーで観測した各ペアの
            // 「既に見えている長さ」のリスト。
            var l_ローカルペア経路 = new Dictionary<(int, int), List<int>>[l_スレッド数];
            // インサートサイズ自動推定用: 両リードが同一unitigに単独マップされた
            // ペアについて、そのunitig内での距離を標本抽出する。
            // ライブラリの向きの組み合わせ(FR/RF/FF/RR)を決め打ちできないため、
            // 「符号が一致するヒット」と「符号が不一致のヒット」を別々に集計し、
            // 全体マージ後に標本数が多い方を採用する。
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
            Console.WriteLine($"[Info] Paired-end adjacency candidates detected: {this._ペア経路.Count} edges ({l_ペア支持数} supporting pairs total).");
            Console.WriteLine($"[Info] Same-unitig pair orientation counts: same-orientation={l_同一向き合計}, opposite-orientation={l_逆向き合計}. Using '{l_採用ラベル}' as the library's observed orientation for InsertSize estimation ({l_同一ユニティグ標本.Count} samples).");
            if (l_同一ユニティグ標本.Count > 0)
            {
                // 同一unitig内標本は、unitig自体がフラグメント長より短い場合
                // 両端が同じunitig内に収まるペアしか観測できず、より短い
                // フラグメントに偏った標本になりやすい(unitigが短いほど顕著)。
                Console.WriteLine($"[Info] Same-unitig fragment-length distribution: {Get_分布要約(l_同一ユニティグ標本)}.");
                Console.WriteLine($"[Info] Same-unitig fragment-length median: {Get_中央値(l_同一ユニティグ標本)} (from {l_同一ユニティグ標本.Count} samples; read lengths added back to the inner distance, so this is a true fragment length. May still be biased short if unitigs are shorter than the true insert size).");
            }
        }

        /// <summary>
        /// 両リードが同一 unitig にマップされたペアから、フラグメント長の標本を取る。
        ///
        /// 順鎖座標への変換が返すのは「そのヒットの向きで見た既知長」を順鎖座標へ
        /// 写した値、すなわち順鎖から見たリードの「内側の端」の座標になる。
        /// したがって2つの差はフラグメント長ではなく、2リードに挟まれた内側の
        /// 未読区間(inner distance)の長さである。
        ///
        /// 実データ(150bpリード・IS350ライブラリ)でこの差の中央値が58になり、
        /// リード長150bpより短いという物理的にありえない推定値になっていた。
        /// 内側距離 d と真のフラグメント長 F の関係は FR配置で
        /// F = d + len(read1) + len(read2) であり、58 + 150 + 150 = 358 で
        /// ライブラリ名(IS350)と一致する。ここで両リード長を足し戻し、
        /// 以降の推定値が一貫して「フラグメント長」の単位になるようにする。
        /// </summary>
        private static void V_収集_同一ユニティグ標本(
            代表ユニティグヒット p_ヒット1, 代表ユニティグヒット p_ヒット2,
            string p_リード1, string p_リード2,
            List<int> p_同一向き標本, List<int> p_逆向き標本)
        {
            if ((p_ヒット1.A_ユニティグID > 0) == (p_ヒット2.A_ユニティグID > 0))
            {
                // 同じ向き同士(FF/RR相当)。両リードの内側の端はどちらも
                // 同じ側を向いているため、差は「開始位置の差」に相当する。
                // 下流側リード1本分を足すとフラグメント長になる。
                var l_内側距離 = Math.Abs(Get_順鎖座標(p_ヒット1) - Get_順鎖座標(p_ヒット2));
                var l_フラグメント長 = l_内側距離 + Math.Max(p_リード1.Length, p_リード2.Length);
                if (l_フラグメント長 > 0)
                {
                    p_同一向き標本.Add(l_フラグメント長);
                }
            }
            else
            {
                // 互いに逆向き(FR相当、Illuminaペアエンドの通常配置)。
                // 順鎖側ヒットのリードがフラグメントの左端、
                // 逆鎖側ヒットのリードが右端を占める。
                var l_ヒット1が順鎖か = p_ヒット1.A_ユニティグID > 0;
                var l_順鎖側の端 = Get_順鎖座標(l_ヒット1が順鎖か ? p_ヒット1 : p_ヒット2);
                var l_逆鎖側の端 = Get_順鎖座標(l_ヒット1が順鎖か ? p_ヒット2 : p_ヒット1);
                var l_順鎖側リード長 = l_ヒット1が順鎖か ? p_リード1.Length : p_リード2.Length;
                var l_逆鎖側リード長 = l_ヒット1が順鎖か ? p_リード2.Length : p_リード1.Length;

                // フラグメントの左端 = 順鎖リードの開始位置、
                // 右端 = 逆鎖リードの終了位置。
                var l_フラグメント長 = (l_逆鎖側の端 + l_逆鎖側リード長) - (l_順鎖側の端 - l_順鎖側リード長);
                if (l_フラグメント長 > 0)
                {
                    p_逆向き標本.Add(l_フラグメント長);
                }
            }
        }

        /// <summary>
        /// 両リードが別々の unitig にマップされたペアを、隣接候補として記録する。
        ///
        /// read2 はリード分子の逆鎖側から読まれるため、read1 の向きに揃えるには
        /// read2 の逆相補鎖が実際に「read1 の下流」に来る、という関係になる。
        /// read2 側でヒットした unitig ID の符号を反転させ、read1 の向きに
        /// 揃えたうえでペアを記録する。
        ///
        /// 記録するのは「フラグメントのうち、2つのunitigの内側に既に見えている分の長さ」:
        ///   read1の長さ + unitig1末端までの残り + unitig2先頭からの残り + read2の長さ
        /// 未知区間(ギャップ)長を G とすると
        ///   フラグメント長 = この値 + G
        /// という関係が常に成り立つ(直接k-1で結合された場合は G = -(k-1))。
        ///
        /// 以前は read1/read2 の長さを含めない残り長だけを記録していたため、
        /// ここから逆算されるインサートサイズもギャップ長も
        /// 両リード長ぶん(実データで300bp)ずれていた。
        /// </summary>
        private static void V_収集_ペア経路(
            代表ユニティグヒット p_ヒット1, 代表ユニティグヒット p_ヒット2,
            string p_リード1, string p_リード2,
            Dictionary<(int, int), List<int>> p_ローカルペア経路)
        {
            var l_キー = (p_ヒット1.A_ユニティグID, -p_ヒット2.A_ユニティグID);

            var l_残り1 = p_ヒット1.A_末尾までの残り長;
            var l_残り2 = Get_反転後の残り長(p_ヒット2);
            var l_既知長 = l_残り1 + l_残り2 + p_リード1.Length + p_リード2.Length;

            if (p_ローカルペア経路.TryGetValue(l_キー, out var l_一覧))
            {
                l_一覧.Add(l_既知長);
            }
            else
            {
                p_ローカルペア経路[l_キー] = [l_既知長];
            }
        }

        /// <summary>
        /// ヒット(ある unitig ID の符号が示す向きで計算された最終一致終端位置)を、
        /// その unitig を「逆向き」に見た座標系での残り長に変換する。
        /// 元の向きでの残り長(終端までの残り)が、逆向きで見たときの
        /// 「先頭からの既知長」に相当するため、逆向きでの残り長は
        /// 元の向きでの最終一致終端位置(先頭からの既知長)がそのまま使える。
        /// </summary>
        private static int Get_反転後の残り長(代表ユニティグヒット p_ヒット)
        {
            return Math.Max(0, p_ヒット.A_最終一致終端位置);
        }

        /// <summary>
        /// ヒットの最終一致終端位置(ヒット自身の符号が示す向きの座標系での値)を、
        /// 常に unitig の「順鎖」座標系での位置に変換する。順鎖ヒットはそのまま、
        /// 逆鎖ヒットは ユニティグ長 - 最終一致終端位置 に変換する。
        /// 同一unitig上の2ヒット間の距離を求める際、座標系を揃えるために使う。
        /// </summary>
        private static int Get_順鎖座標(代表ユニティグヒット p_ヒット)
        {
            return p_ヒット.A_ユニティグID > 0
                ? p_ヒット.A_最終一致終端位置
                : p_ヒット.A_ユニティグ長 - p_ヒット.A_最終一致終端位置;
        }

        /// <summary>
        /// インサートサイズ自動推定用に標本抽出された距離の一覧
        /// (同一ユニティグ標本と確定辺標本の結合)。
        /// 後方互換のため残しているが、Scaffolder は標本の出所による
        /// バイアスの違いを考慮するため個別の一覧を優先的に参照する。
        /// </summary>
        public List<int> A_インサートサイズ標本 { get; } = [];

        /// <summary>
        /// 単一unitig内で両リードがヒットしたペアからの標本。
        /// unitig自体がフラグメント長より短い場合、両端が収まるペアしか
        /// 観測できないため、より短いフラグメントに偏りやすい。
        /// </summary>
        public List<int> A_同一ユニティグ標本 { get; } = [];

        /// <summary>
        /// unitig同士がk-1オーバーラップで直接結合されたペアからの標本。
        /// 同一ユニティグ標本のような長さバイアスを受けない。
        /// </summary>
        public List<int> A_確定辺標本 { get; } = [];

        /// <summary>
        /// ペアエンド由来の隣接候補。キーは (始点, 終点) の unitig ID(符号は向き)、
        /// 値は各観測ペアの既知長の一覧。Scaffolder から参照される。
        /// </summary>
        public IReadOnlyDictionary<(int, int), List<int>> A_ペア経路 => this._ペア経路;

        /// <summary>
        /// unitig ID(1始まり、符号なし)からその塩基長を引く。Scaffolder が
        /// contig 側の末端 unitig の長さを参照する際に使う。
        /// </summary>
        public IReadOnlyDictionary<int, int> A_ユニティグ長 => this._ユニティグ長;

        public IReadOnlyDictionary<int, ユニティグ配置> A_ユニティグ配置 => this._ユニティグ配置;

        /// <summary>
        /// 1本のリードが「代表として」どの unitig にマップされるかを判定する。
        /// リード中で最も安定して(連続して)ヒットし続けた unitig ID に加え、
        /// スキャフォールディングのギャップ長推定に使うための最終ヒット位置も返す。
        /// どの unitig にもヒットしなかった場合はヒットなしを返す。
        /// ペアエンドの隣接検出でのみ使用する軽量な単方向スキャン。
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
            // 各候補 ID ごとに、その ID として最後にヒットしたk-merの
            // 「unitig内での」終端位置を記録する。k-mer辞書 が
            // (unitigID, unitig内開始位置) を保持するようになったため、
            // read内での相対位置ではなく辞書から得た本物のunitig内位置を使う
            // (以前はread内終端位置をそのままunitig内終端位置として誤用しており、
            //  unitigがread長より十分短い場合はたまたま近い値になり問題が
            //  表面化しにくかったが、unitigが長くなると全く違う値になっていた)。
            var l_最終終端位置 = new Dictionary<int, int>();
            var l_曖昧塩基数 = 0;
            for (var i = 0; i < l_k長; i++)
            {
                if (Util.Get_塩基ID候補(p_リード[i]).Count > 1)
                {
                    l_曖昧塩基数++;
                }
            }
            for (var i = l_k長; i <= p_リード.Length; i++)
            {
                if (Util.Get_塩基ID候補(p_リード[i - l_k長]).Count > 1)
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
        /// 1リード分の k-mer マッピングを行い、隣接関係を(スレッドローカルな)
        /// 辞書に集計する。
        /// </summary>
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
                        Console.WriteLine($"[Warning] Paired read IDs do not match at this position (\"{l_データ1.A_ID}\" vs \"{l_データ2.A_ID}\"). " +
                            "Paired-end adjacency detection may be unreliable for reads after this point; " +
                            "single-read adjacency detection is unaffected.");
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
                if (Util.Get_塩基ID候補(p_リード[i]).Count > 1)
                {
                    l_曖昧塩基数++;
                }
                if (Util.Get_塩基ID候補(l_逆鎖リード[i]).Count > 1)
                {
                    l_逆鎖の曖昧塩基数++;
                }
            }
            for (var i = l_k長; i <= p_リード.Length; i++)
            {
                if (Util.Get_塩基ID候補(p_リード[i - l_k長]).Count > 1)
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
                if (Util.Get_塩基ID候補(l_逆鎖リード[i - l_k長]).Count > 1)
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

        /// <param name="p_コピー数">
        /// unitig ID -> 推定コピー数。先読み探索で「この unitig を何回まで通ってよいか」の
        /// 予算に使う。渡さない場合はすべて1コピーとして扱い、先読み探索も控えめになる。
        /// </param>
        public void V_結合_コンティグ(
            string p_コンティグパス,
            decimal p_優勢閾値,
            ulong p_最小証拠数,
            IReadOnlyDictionary<int, int>? p_コピー数 = null)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            var l_重なり長 = l_k長 - 1;

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

            // 隣接は de Bruijn グラフから厳密に導く(UnitigGraph の説明を参照)。
            // リードマッピング由来の隣接情報は「辺を作る」ためではなく、
            // 分岐点でどの辺を選ぶかの「重み」としてのみ使う。
            var l_グラフ = UnitigGraph.Get_グラフ(l_ユニティグ配列, this._kmer辞書, l_k長, 曖昧kmerの番兵);

            var l_辺数 = 0;
            var l_分岐頂点数 = 0;
            for (var v = 2; v < l_グラフ.A_出辺.Count; v++)
            {
                l_辺数 += l_グラフ.A_出辺[v].Count;
                if (l_グラフ.A_出辺[v].Count > 1)
                {
                    l_分岐頂点数++;
                }
            }
            Console.WriteLine($"[Debug] Exact de Bruijn unitig graph: {l_辺数} directed edge(s), {l_分岐頂点数} branching vertex(es) out of {l_グラフ.A_出辺.Count - 2}.");

            // リード由来の支持数を逆鎖対称に集計する。辺 v→w と w^1→v^1 は
            // 同一の物理的な隣接を表すため、重みも同一でなければ順鎖側と
            // 逆鎖側で異なる経路が選ばれ、同じ領域が 2 通りに組み立てられてしまう。
            Dictionary<(int, int), ulong> l_支持 = [];
            foreach (var ((l_始点, l_終点), l_件数) in this._リード隣接)
            {
                if (l_始点 == l_終点)
                {
                    continue;
                }
                var v = Get_頂点番号(l_始点);
                var w = Get_頂点番号(l_終点);
                l_支持[(v, w)] = l_支持.GetValueOrDefault((v, w)) + l_件数;
                l_支持[(w ^ 1, v ^ 1)] = l_支持.GetValueOrDefault((w ^ 1, v ^ 1)) + l_件数;
            }

            // ペアエンド由来の支持も同じ重みに足し込む。1本のリード(150bp)では
            // それより長い反復配列を跨げず、分岐でどちらへ進むべきか決められないが、
            // フラグメント(実データで約350bp)なら跨げる場合がある。
            //
            // ここで参照されるのは、あくまで de Bruijn グラフ上に実在する辺に
            // ついての値だけである。フラグメント長ぶん離れているだけで隣接して
            // いない unitig 対もペア経路には入るが、それらは辺が無いため
            // 選択に影響しない。
            Dictionary<(int, int), ulong> l_ペア連結 = [];
            var l_ペア支持を足した数 = 0;
            foreach (var ((l_始点, l_終点), l_標本) in this._ペア経路)
            {
                if (l_始点 == l_終点)
                {
                    continue;
                }
                var v = Get_頂点番号(l_始点);
                var w = Get_頂点番号(l_終点);
                var l_件数 = (ulong)l_標本.Count;
                l_支持[(v, w)] = l_支持.GetValueOrDefault((v, w)) + l_件数;
                l_支持[(w ^ 1, v ^ 1)] = l_支持.GetValueOrDefault((w ^ 1, v ^ 1)) + l_件数;
                l_ペア連結[(v, w)] = l_ペア連結.GetValueOrDefault((v, w)) + l_件数;
                l_ペア連結[(w ^ 1, v ^ 1)] = l_ペア連結.GetValueOrDefault((w ^ 1, v ^ 1)) + l_件数;
                l_ペア支持を足した数++;
            }
            Console.WriteLine($"[Debug] Branch-selection weights: {this._リード隣接.Count} single-read adjacency pair(s) + {l_ペア支持を足した数} paired-end pair(s).");

            // 単純バブルを潰してから辺を選ぶ。相互一意性を課す以上、
            // 再合流点の入次数が2以上のまま残っているとその経路全体が
            // 結合されなくなるため、先に枝を1本に絞っておく必要がある。
            var l_除去バブル数 = l_グラフ.V_除去_単純バブル(l_ユニティグ配列, l_支持);
            if (l_除去バブル数 > 0)
            {
                Console.WriteLine($"[Debug] Popped {l_除去バブル数} simple bubble branch(es) (kept as standalone contigs; only their graph edges were removed).");
            }

            // 跨げる見込みのある長さの上限。フラグメント長の実測中央値を使う
            // (これより長い反復は、そもそも両端を別々の unitig に載せた
            //  ペアが存在しえない)。標本が無い場合は控えめな既定値。
            var l_反復長の上限 = this.A_同一ユニティグ標本.Count > 0
                ? Get_中央値(this.A_同一ユニティグ標本)
                : l_k長 * 4;
            var l_解決した反復数 = l_グラフ.V_解決_短い反復(
                l_ユニティグ配列, l_支持, l_ペア連結, l_反復長の上限, p_優勢閾値, p_最小証拠数);
            Console.WriteLine(
                $"[Debug] Repeat resolution: {l_解決した反復数} short repeat(s) (<= {l_反復長の上限}bp) were duplicated " +
                "and untangled using read pairs that span them.");

            // 各頂点について「出て行く先」を高々 1 つに絞る。
            var l_選択 = new int[l_グラフ.A_出辺.Count];
            Array.Fill(l_選択, -1);
            var l_一意な頂点数 = 0;
            var l_支持で解決した数 = 0;
            var l_反復由来で未解決の数 = 0;
            for (var v = 2; v < l_グラフ.A_出辺.Count; v++)
            {
                var l_出辺 = l_グラフ.A_出辺[v];
                if (l_出辺.Count == 0)
                {
                    continue;
                }
                if (l_出辺.Count == 1)
                {
                    l_選択[v] = l_出辺[0];
                    l_一意な頂点数++;
                    continue;
                }

                // 分岐元が多コピー(反復配列)の場合、そこから出るリード支持は
                // どのコピー由来か区別できないため、行き先を選ぶ根拠にならない。
                // 反復の各コピーの続きが全部この1頂点に集まっているので、
                // 支持はすべての行き先に付いてしまう。
                //
                // これを見落として支持で選んでいたため、A-R-B-R-C という構造の
                // 合成ゲノム(R は2コピーの反復)で、R から出る分岐を誤って選び
                // A-R-C という「中間を飛ばして端同士を繋いだ」contig を出力していた
                // (真値照合で発覚)。正しい続きは、反復の外側にある単一コピー領域
                // からのペアエンドでしか決められない。
                if (p_コピー数 is not null && p_コピー数.GetValueOrDefault(v >> 1, 1) > 1)
                {
                    l_反復由来で未解決の数++;
                    continue;
                }

                var l_合計 = 0UL;
                var l_最良 = -1;
                var l_最良の支持 = 0UL;
                foreach (var w in l_出辺)
                {
                    var l_件数 = l_支持.GetValueOrDefault((v, w));
                    l_合計 += l_件数;
                    if (l_件数 > l_最良の支持)
                    {
                        l_最良の支持 = l_件数;
                        l_最良 = w;
                    }
                }
                if (l_最良 >= 0 && l_最良の支持 >= p_最小証拠数 && l_合計 > 0
                    && (decimal)l_最良の支持 / l_合計 >= p_優勢閾値)
                {
                    l_選択[v] = l_最良;
                    l_支持で解決した数++;
                }
            }

            // 相互一意な辺だけを実際の結合として採用する。
            // v→w を結合してよいのは「v の唯一の行き先が w」であり、かつ
            // 「w の唯一の来訪元が v」であるときに限る。後者は逆鎖対称性より
            // 選択[w^1] == v^1 と同値。この条件を欠くと、複数の異なる
            // unitig が同じ次の unitig を指し(実データで 1550 頂点)、
            // 先着 1 本だけが結合されて残りが千切れる形になっていた。
            var l_結合 = new int[l_グラフ.A_出辺.Count];
            Array.Fill(l_結合, -1);
            var l_結合数 = 0;
            var l_反復通り抜けで棄却した数 = 0;
            for (var v = 2; v < l_グラフ.A_出辺.Count; v++)
            {
                var w = l_選択[v];
                if (w < 0 || l_選択[w ^ 1] != (v ^ 1))
                {
                    continue;
                }
                // 結合は逆鎖側と対で成立する。片側だけ許すと結合の対称性が
                // 崩れ、walk の始点判定が壊れるため、
                // どちらかが通り抜け不可なら対ごと採用しない。
                if (!Get_通り抜けてよいか(l_グラフ, p_コピー数, v) || !Get_通り抜けてよいか(l_グラフ, p_コピー数, w ^ 1))
                {
                    l_反復通り抜けで棄却した数++;
                    continue;
                }
                l_結合[v] = w;
                l_結合数++;
            }
            if (l_反復通り抜けで棄却した数 > 0)
            {
                Console.WriteLine(
                    $"[Debug] {l_反復通り抜けで棄却した数 / 2} join(s) were refused because they would chain through a multi-copy " +
                    "repeat that has not been untangled (doing so skips whatever lies between the repeat's copies).");
            }
            Console.WriteLine($"[Debug] Edge selection: {l_一意な頂点数} vertex(es) had a single out-edge, {l_支持で解決した数} branch(es) resolved by read support, " +
                $"{l_反復由来で未解決の数} branch(es) left unresolved because they leave a multi-copy repeat (reads inside a repeat cannot tell the copies apart); {l_結合数} directed merge(s) survived the mutual-uniqueness check ({l_結合数 / 2} undirected join(s)).");

            // 1歩だけを見る相互一意性の判定では決めきれなかった分岐を、
            // 数kb先まで複数経路を並行して伸ばして(ビームサーチ)解けるだけ解く。
            // 分岐の直後だけを見ると五分五分でも、少し先まで進めると片方だけが
            // ペアエンドの証拠と整合する、という状況を拾える。
            var l_先読みで解決した数 = BeamSearchExtender.V_延長_先読み(
                l_グラフ,
                l_ユニティグ配列,
                l_結合,
                l_ペア連結,
                p_コピー数 ?? new Dictionary<int, int>(),
                p_インサートサイズ: l_反復長の上限,
                p_優勢閾値: p_優勢閾値,
                p_最小証拠数: p_最小証拠数);
            if (l_先読みで解決した数 > 0)
            {
                Console.WriteLine(
                    $"[Debug] Beam-search lookahead resolved {l_先読みで解決した数 / 2} further junction(s) that the " +
                    "single-step mutual-uniqueness rule could not decide.");
            }

            this.V_収集_確定辺標本(l_結合);

            // 双子(v と v^1)は同一 unitig の裏表なので、unitig 単位で訪問済みを
            // 管理する。これを頂点単位でやっていたため、順鎖側の walk と逆鎖側の
            // walk が同じ unitig を別々に出力し、contig 総長が unitig 総長の
            // ちょうど 2 倍に膨れていた。
            var l_ユニティグ数 = (l_ユニティグ配列.Count - 2) / 2;
            var l_訪問済み = new bool[l_ユニティグ数 + 1];

            List<string> l_コンティグ群 = [];
            List<List<int>> l_walk順群 = [];
            List<bool> l_環状フラグ群 = [];

            // 結合グラフ上で「入ってくる結合を持たない」頂点が経路の始点。
            // v への結合が存在することは、逆鎖対称性より 結合[v^1] != -1 と同値。
            for (var v = 2; v < l_グラフ.A_出辺.Count; v++)
            {
                if (l_結合[v ^ 1] != -1 || l_訪問済み[v >> 1])
                {
                    continue;
                }
                V_実行_walk(l_ユニティグ配列, l_結合, l_訪問済み, l_重なり長, v, l_コンティグ群, l_walk順群, l_環状フラグ群);
            }

            // 始点を持たない=循環している経路を拾う(環状ゲノム/プラスミド等)。
            for (var v = 2; v < l_グラフ.A_出辺.Count; v += 2)
            {
                if (l_訪問済み[v >> 1])
                {
                    continue;
                }
                V_実行_walk(l_ユニティグ配列, l_結合, l_訪問済み, l_重なり長, v, l_コンティグ群, l_walk順群, l_環状フラグ群);
            }

            using var l_書き込み = new FastaWriter(p_コンティグパス);
            var l_ID = 1;
            var l_総延長 = 0L;
            for (var c = 0; c < l_コンティグ群.Count; c++)
            {
                var l_コンティグ = l_コンティグ群[c];
                var l_walk順 = l_walk順群[c];
                var l_逆相補 = Util.V_逆相補(l_コンティグ);
                var l_逆相補を採用するか = string.CompareOrdinal(l_コンティグ, l_逆相補) > 0;
                // 環状に閉じた contig は、その複製単位(染色体・プラスミド)を
                // 完全に組み上げられたことを意味するため、名前に明示する。
                var l_名前 = l_環状フラグ群[c] ? $"NODE{l_ID}_circular" : $"NODE{l_ID}";
                l_書き込み.V_書き込み(l_名前, l_逆相補を採用するか ? l_逆相補 : l_コンティグ);

                // walk順に含まれる各頂点(unitig の向き付き番号)を配置として記録する。
                // walk順は「実際に配列へ連結された順」なので、そのままこの contig 内での
                // 並び順になる。逆相補を採用した場合、contigs.fasta 上の配列は
                // walk 順と逆向きになっているため、位置(先頭/末尾)の解釈は
                // Scaffolder 側で反転させる。
                for (var w = 0; w < l_walk順.Count; w++)
                {
                    var l_頂点番号 = l_walk順[w];
                    this._ユニティグ配置[l_頂点番号 >> 1] = new ユニティグ配置(
                        p_コンティグID: l_ID,
                        p_コンティグが逆相補か: l_逆相補を採用するか,
                        p_walk順の位置: w,
                        p_walk順の総数: l_walk順.Count,
                        p_walk中で逆鎖か: (l_頂点番号 & 1) == 1);
                }

                l_ID++;
                l_総延長 += l_コンティグ.Length;
            }
            Console.WriteLine("Total Length of contigs : " + l_総延長);

            var l_環状コンティグ = Enumerable.Range(0, l_コンティグ群.Count).Where(x => l_環状フラグ群[x]).ToList();
            if (l_環状コンティグ.Count > 0)
            {
                var l_長さ一覧 = string.Join(", ", l_環状コンティグ.Select(x => $"{l_コンティグ群[x].Length}bp"));
                Console.WriteLine(
                    $"[Info] {l_環状コンティグ.Count} contig(s) closed into a circle ({l_長さ一覧}). " +
                    "A closed circle means that replicon (chromosome or plasmid) was assembled end to end.");
            }
            else
            {
                Console.WriteLine("[Info] No contig closed into a circle; every replicon is still fragmented.");
            }
        }

        /// <summary>
        /// その頂点を「通り抜けて」よいか。
        ///
        /// 反復配列 R がゲノム中に2回現れ A-R-B-R-C という並びだった場合、
        /// A→R と R→C はどちらも個別には本物の隣接である(前者は1つ目の
        /// コピー、後者は2つ目のコピー)。しかし walk は各 unitig を1回しか
        /// 使えないため、この2つを R 経由で連鎖させると A-R-C という
        /// 「中間の B を飛ばして端同士を繋いだ」配列ができてしまう。
        /// 合成ゲノムの真値照合で実際にこれが起きていた。
        ///
        /// 通り抜けてよいのは、反復解決によって解きほぐされ、入次数・出次数が
        /// どちらも1になった場合だけ。その状態なら「どのコピーにいるか」が確定している。
        /// </summary>
        private static bool Get_通り抜けてよいか(
            UnitigGraph p_グラフ, IReadOnlyDictionary<int, int>? p_コピー数, int p_頂点)
        {
            if ((p_コピー数?.GetValueOrDefault(p_頂点 >> 1, 1) ?? 1) <= 1)
            {
                return true;
            }
            return p_グラフ.A_出辺[p_頂点].Count == 1 && p_グラフ.Get_入次数(p_頂点) == 1;
        }

        /// <summary>
        /// 始点から結合を辿って1本の contig を組み立て、結果を各一覧へ追加する。
        /// </summary>
        private static void V_実行_walk(
            List<string> p_ユニティグ配列, int[] p_結合, bool[] p_訪問済み, int p_重なり長, int p_始点,
            List<string> p_コンティグ群, List<List<int>> p_walk順群, List<bool> p_環状フラグ群)
        {
            List<int> l_walk順 = [];
            var (l_配列, l_環状か) = Get_walk結果(p_ユニティグ配列, p_結合, p_訪問済み, p_重なり長, p_始点, l_walk順);
            p_コンティグ群.Add(l_配列);
            p_walk順群.Add(l_walk順);
            p_環状フラグ群.Add(l_環状か);
        }

        /// <summary>
        /// 始点から結合を辿って配列を組み立てる。経路が始点へ戻ってきた場合は
        /// 環状として報告する。
        /// </summary>
        private static (string A_配列, bool A_環状か) Get_walk結果(
            List<string> p_ユニティグ配列, int[] p_結合, bool[] p_訪問済み, int p_重なり長, int p_始点,
            List<int> p_walk順)
        {
            var l_出力 = new StringBuilder(p_ユニティグ配列[p_始点]);
            p_walk順.Add(p_始点);
            p_訪問済み[p_始点 >> 1] = true;
            var l_現在 = p_始点;
            var l_環状か = false;
            while (true)
            {
                var l_次 = p_結合[l_現在];
                if (l_次 < 0)
                {
                    break;
                }
                if (p_訪問済み[l_次 >> 1])
                {
                    // 始点へ戻ってきた = 経路が閉じている。細菌の染色体と
                    // プラスミドは環状なので、これは「その複製単位を
                    // 完全に1周組み上げられた」ことを意味する。
                    l_環状か = l_次 == p_始点;
                    break;
                }
                var l_配列 = p_ユニティグ配列[l_次];
                if (l_配列.Length < p_重なり長 || l_出力.Length < p_重なり長)
                {
                    break;
                }
                // 構築方法より k-1 のオーバーラップは保証されているが、
                // 万一崩れていた場合に誤った配列を作らないよう検証する。
                if (!Get_重なりが一致するか(l_出力, l_配列, p_重なり長))
                {
                    break;
                }
                _ = l_出力.Append(l_配列[p_重なり長..]);
                p_訪問済み[l_次 >> 1] = true;
                p_walk順.Add(l_次);
                l_現在 = l_次;
            }

            // 環状の場合、末尾 unitig は「始点 unitig と重なる k-1 塩基」を
            // 自分の末尾に含んでいる。その k-1 塩基は配列の先頭にも現れて
            // いるので、そのまま出すと円周が k-1 塩基ぶん長くなってしまう。
            // 線状 contig の連結では次の unitig 側から重なりを取り除いて
            // いるが、環状の場合は「次」が既に出力済みの始点なので、
            // ここで末尾から取り除く。
            if (l_環状か && l_出力.Length > p_重なり長)
            {
                _ = l_出力.Remove(l_出力.Length - p_重なり長, p_重なり長);
            }

            return (l_出力.ToString(), l_環状か);
        }

        /// <summary>
        /// 相互一意性の検査を通って実際に結合が確定した辺を走査し、
        /// その関係をペア経路のキー形式(符号付き unitig ID のペア)に変換する。
        /// 変換後、該当するペア経路エントリの「既知長」から
        /// フラグメント長 = 既知長 - (k-1) を計算して標本に積む。
        /// </summary>
        private void V_収集_確定辺標本(int[] p_結合)
        {
            var l_重なり長 = ConfigurationManager.A_実行時引数.A_k長 - 1;
            List<int> l_確定辺標本 = [];

            for (var v = 2; v < p_結合.Length; v++)
            {
                var l_次 = p_結合[v];
                if (l_次 < 0)
                {
                    continue;
                }

                // 頂点番号 -> 符号付き unitig ID。
                var l_始点ユニティグ = (v >> 1) * ((v & 1) == 0 ? 1 : -1);
                var l_終点ユニティグ = (l_次 >> 1) * ((l_次 & 1) == 0 ? 1 : -1);

                if (!this._ペア経路.TryGetValue((l_始点ユニティグ, l_終点ユニティグ), out var l_既知長標本))
                {
                    continue;
                }

                foreach (var l_既知長 in l_既知長標本)
                {
                    // 直接結合された辺では2つのunitigがk-1塩基重なるので、
                    // 未知区間の長さは G = -(k-1)。よって
                    // フラグメント長 = 既知長 - (k-1)。
                    var l_フラグメント長 = l_既知長 - l_重なり長;
                    if (l_フラグメント長 > 0)
                    {
                        l_確定辺標本.Add(l_フラグメント長);
                    }
                }
            }

            this.A_インサートサイズ標本.AddRange(l_確定辺標本);
            this.A_確定辺標本.AddRange(l_確定辺標本);

            Console.WriteLine($"[Info] InsertSize samples derived from resolved (actually-joined) unitig adjacency: {l_確定辺標本.Count}.");
            if (l_確定辺標本.Count > 0)
            {
                // このプールは「unitig同士がk-1オーバーラップで直接結合された」
                // ペアのみを対象とするため、同一unitig標本のような
                // 「フラグメントが1つのunitigに収まる必要がある」制約が
                // なく、短いunitigによる短フラグメントへの偏りを受けにくい。
                Console.WriteLine($"[Info] Resolved-edge sample median: {Get_中央値(l_確定辺標本)} (from {l_確定辺標本.Count} samples; not subject to the same-unitig length bias).");
            }
        }

        /// <summary>
        /// 符号付き unitig ID(正=順鎖、負=逆鎖)をグラフの頂点番号に変換する。
        /// </summary>
        internal static int Get_頂点番号(int p_符号付きユニティグID)
        {
            return (Math.Abs(p_符号付きユニティグID) << 1) | (p_符号付きユニティグID > 0 ? 0 : 1);
        }

        /// <summary>
        /// 出力の末尾 p_重なり長 文字と unitig の先頭 p_重なり長 文字が
        /// 一致するかどうかを判定する。
        /// </summary>
        private static bool Get_重なりが一致するか(StringBuilder p_出力, string p_ユニティグ, int p_重なり長)
        {
            var l_開始位置 = p_出力.Length - p_重なり長;
            for (var j = 0; j < p_重なり長; j++)
            {
                if (p_出力[l_開始位置 + j] != p_ユニティグ[j])
                {
                    return false;
                }
            }
            return true;
        }

        private static int Get_中央値(List<int> p_値一覧)
        {
            var l_整列済み = p_値一覧.OrderBy(x => x).ToList();
            var l_中央 = l_整列済み.Count / 2;
            return l_整列済み.Count % 2 == 0 ? (l_整列済み[l_中央 - 1] + l_整列済み[l_中央]) / 2 : l_整列済み[l_中央];
        }

        /// <summary>整列済みの一覧から分位点の値を取り出す。</summary>
        private static int Get_分位点(List<int> p_整列済み, double p_分位)
        {
            return p_整列済み[Math.Clamp((int)(p_分位 * (p_整列済み.Count - 1)), 0, p_整列済み.Count - 1)];
        }

        /// <summary>
        /// フラグメント長分布の分位点を要約する。中央値だけでは
        /// 「このライブラリがどれだけの長さのギャップを跨げるか」が分からない。
        /// リード長の2倍を超える分だけがスキャフォールディングで橋渡しできる
        /// 未知区間の長さなので、分布の裾(特に上側)が実際の橋渡し能力を決める。
        /// </summary>
        private static string Get_分布要約(List<int> p_値一覧)
        {
            var l_整列済み = p_値一覧.OrderBy(x => x).ToList();
            return $"p1={Get_分位点(l_整列済み, 0.01)}, p10={Get_分位点(l_整列済み, 0.10)}, " +
                $"p25={Get_分位点(l_整列済み, 0.25)}, p50={Get_分位点(l_整列済み, 0.50)}, " +
                $"p75={Get_分位点(l_整列済み, 0.75)}, p90={Get_分位点(l_整列済み, 0.90)}, " +
                $"p99={Get_分位点(l_整列済み, 0.99)}, max={l_整列済み[^1]}";
        }
    }
}
