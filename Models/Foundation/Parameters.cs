using Tsumiki.Commons;

namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// 実行時引数
    /// </summary>
    internal class Parameters
    {
        #region 定数

        /// <summary>
        /// インサートサイズ未指定表示
        /// </summary>
        private const string インサートサイズ未指定表示 = "unspecified";

        /// <summary>
        /// k-mer カットオフの既定値
        /// </summary>
        private const int kmerカットオフの既定値 = 2;

        /// <summary>
        /// k 長の既定値
        /// </summary>
        private const int k長の既定値 = 31;

        #endregion

        #region 内部変数

        /// <summary>
        /// リード 1 のパス
        /// </summary>
        private List<string> _リード1のパス群 = [];

        /// <summary>
        /// リード 2 のパス
        /// </summary>
        private List<string> _リード2のパス群 = [];

        /// <summary>
        /// シングルエンドのライブラリのパス
        /// </summary>
        private List<string> _シングルのパス群 = [];

        /// <summary>
        /// 前処理などで差し替えたライブラリ、未差し替えなら null
        /// </summary>
        private List<(string A_リード1, string A_リード2)>? _差し替えたライブラリ群;

        /// <summary>
        /// k 長
        /// </summary>
        private int _k長 = k長の既定値;

        /// <summary>
        /// -k に指定された k の一覧
        /// </summary>
        private List<int> _k長一覧 = [];

        /// <summary>
        /// k-mer カットオフ
        /// </summary>
        private ulong _kmerカットオフ = kmerカットオフの既定値;

        /// <summary>
        /// Phred オフセット
        /// </summary>
        private int _Phredオフセット = Consts.Phredオフセットの既定値;

        /// <summary>
        /// ライブラリごとの Phred オフセット
        /// </summary>
        /// <remarks>
        /// 符号化がライブラリで違うと、1 つの値で読んだ側の品質判定が静かに狂う
        /// </remarks>
        private List<int> _Phredオフセット群 = [];

        /// <summary>
        /// ライブラリごとの代表リード長
        /// </summary>
        private List<int> _ライブラリのリード長 = [];

        /// <summary>
        /// 並列に使うスレッド数
        /// </summary>
        private int _スレッド数 = Environment.ProcessorCount;

        /// <summary>
        /// ペアの支持で結合を確定させる優勢の閾値
        /// </summary>
        private decimal _ペア結合閾値 = Consts.ペア結合閾値の既定値;

        /// <summary>
        /// 結合を確定させるために要求するペアの支持数
        /// </summary>
        private ulong _ペア支持数閾値 = Consts.ペア支持数閾値の既定値;

        #endregion

        #region プロパティ

        /// <summary>
        /// リード 1 のパス (複数ライブラリならカンマ区切り)
        /// </summary>
        /// <remarks>
        /// 1 本だけを要る処理は <see cref="A_ライブラリ群"/> を回すこと<br/>
        /// この値をそのままファイルパスとして使うと、複数ライブラリのときに存在しないパスとして失敗する
        /// </remarks>
        public string A_リード1のパス
        {
            get => string.Join(",", this.A_ライブラリ群.Select(x => x.A_リード1));
            set
            {
                var l_群 = Get_パス群(value);
                foreach (var l_パス in l_群)
                {
                    if (!File.Exists(l_パス))
                    {
                        throw new ArgumentException($"Read1's path {l_パス} is not found");
                    }
                }
                this._リード1のパス群 = l_群;
                this._差し替えたライブラリ群 = null;
            }
        }

        /// <summary>
        /// リード 2 のパス (複数ライブラリならカンマ区切り)
        /// </summary>
        public string A_リード2のパス
        {
            get => string.Join(",", this.A_ライブラリ群.Select(x => string.IsNullOrWhiteSpace(x.A_リード2) ? "-" : x.A_リード2));
            set
            {
                var l_群 = Get_パス群(value);
                foreach (var l_パス in l_群)
                {
                    if (!File.Exists(l_パス))
                    {
                        throw new ArgumentException($"Read2's path {l_パス} is not found");
                    }
                }
                this._リード2のパス群 = l_群;
                this._差し替えたライブラリ群 = null;
            }
        }

        /// <summary>
        /// シングルエンドのライブラリのパス (カンマ区切り)
        /// </summary>
        public string A_シングルのパス
        {
            get => string.Join(",", this._シングルのパス群);
            set
            {
                var l_群 = Get_パス群(value);
                foreach (var l_パス in l_群)
                {
                    if (!File.Exists(l_パス))
                    {
                        throw new ArgumentException($"Single-end read's path {l_パス} is not found");
                    }
                }
                this._シングルのパス群 = l_群;
                this._差し替えたライブラリ群 = null;
            }
        }

        /// <summary>
        /// ライブラリごとのリードの組
        /// </summary>
        /// <remarks>
        /// ペアのライブラリが先、シングルエンドのライブラリが後<br/>
        /// シングルエンドのライブラリはリード 2 が空文字列になる<br/>
        /// -2 を伴わない -1 は、従来どおりシングルエンドとして扱う
        /// </remarks>
        public IReadOnlyList<(string A_リード1, string A_リード2)> A_ライブラリ群
            => this._差し替えたライブラリ群 ?? Get_入力からのライブラリ群();

        /// <summary>
        /// ライブラリの数
        /// </summary>
        public int A_ライブラリ数 => this._リード1のパス群.Count;

        /// <summary>
        /// ペアを持つライブラリがあるか
        /// </summary>
        public bool Hasペア => this.A_ライブラリ群.Any(x => !string.IsNullOrWhiteSpace(x.A_リード2));

        /// <summary>
        /// リード 1 と 2 でライブラリの数が揃っているか
        /// </summary>
        public bool Isライブラリ数が一致 => this._リード2のパス群.Count is 0 || this._リード2のパス群.Count == this._リード1のパス群.Count;

        /// <summary>
        /// ライブラリごとの代表リード長 (観測値)
        /// </summary>
        /// <remarks>
        /// 期待ペア数のモデルはリード長と断片長の差で位置数を数えるので、長さの違うライブラリに 1 つの値を当てると期待が 0 になって証拠が消える
        /// </remarks>
        public IReadOnlyList<int> A_ライブラリのリード長 => this._ライブラリのリード長;

        /// <summary>
        /// ライブラリごとの代表リード長を置く
        /// </summary>
        /// <param name="p_長さ群">ライブラリ順のリード長</param>
        public void Set_ライブラリのリード長(IEnumerable<int> p_長さ群)
        {
            this._ライブラリのリード長 = [.. p_長さ群];
        }

        /// <summary>
        /// -k が明示的に指定されたかどうか
        /// </summary>
        /// <remarks>
        /// 指定されていない場合に限り、実際のリード長から求めた k を自動採用する
        /// </remarks>
        public bool A_Isk長明示指定 { get; private set; }

        /// <summary>
        /// k 長
        /// </summary>
        public int A_k長
        {
            get => this._k長;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of kmer a positive integer");
                }
                this._k長 = value;
                this.A_Isk長明示指定 = true;
            }
        }

        /// <summary>
        /// -k にカンマ区切りで指定された k の一覧 (昇順・重複なし)
        /// </summary>
        /// <remarks>
        /// 未指定なら空<br/>
        /// 2 個以上あれば multi-k として扱う
        /// </remarks>
        public IReadOnlyList<int> A_k長一覧 => this._k長一覧;

        /// <summary>
        /// -kc が明示的に指定されたかどうか
        /// </summary>
        /// <remarks>
        /// 指定されていない場合に限り、k-mer スペクトルの谷から求めたカットオフを自動採用する
        /// </remarks>
        public bool A_Iskmerカットオフ明示指定 { get; private set; }

        /// <summary>
        /// k-mer カットオフ
        /// </summary>
        public ulong A_kmerカットオフ
        {
            get => this._kmerカットオフ;
            set
            {
                if (value <= 0UL)
                {
                    throw new ArgumentException("Please make the value of kmer cut off a positive integer");
                }
                this._kmerカットオフ = value;
                this.A_Iskmerカットオフ明示指定 = true;
            }
        }

        /// <summary>
        /// -p が明示的に指定されたかどうか
        /// </summary>
        /// <remarks>
        /// 指定されていない場合に限り、FASTQ のクオリティ文字列から推定したオフセットを自動採用する<br/>
        /// 明示指定はユーザーの判断なので、推定結果で上書きはしない
        /// </remarks>
        public bool A_IsPhred明示指定 { get; private set; }

        /// <summary>
        /// Phred オフセット
        /// </summary>
        public int A_Phredオフセット
        {
            get => this._Phredオフセット;
            set
            {
                if (!Consts.許容Phredオフセット.Contains(value))
                {
                    throw new ArgumentException($"Phred value is must {string.Join(" or ", Consts.許容Phredオフセット)}");
                }
                this._Phredオフセット = value;
                this._Phredオフセット群 = [];
                this.A_IsPhred明示指定 = true;
            }
        }

        /// <summary>
        /// そのライブラリの Phred オフセット
        /// </summary>
        /// <param name="p_ライブラリ番号">0 起点のライブラリ番号</param>
        /// <remarks>
        /// 明示指定されていれば全ライブラリでその値を使う
        /// </remarks>
        /// <returns></returns>
        public int Get_Phredオフセット(int p_ライブラリ番号)
        {
            return this.A_IsPhred明示指定 || p_ライブラリ番号 >= this._Phredオフセット群.Count
                ? this._Phredオフセット
                : this._Phredオフセット群[p_ライブラリ番号];
        }

        /// <summary>
        /// そのライブラリの Phred オフセットを推定値として置く
        /// </summary>
        /// <param name="p_ライブラリ番号">0 起点のライブラリ番号</param>
        /// <param name="p_オフセット">置くオフセット</param>
        public void Set_推定Phredオフセット(int p_ライブラリ番号, int p_オフセット)
        {
            while (this._Phredオフセット群.Count <= p_ライブラリ番号)
            {
                this._Phredオフセット群.Add(this._Phredオフセット);
            }
            this._Phredオフセット群[p_ライブラリ番号] = p_オフセット;

            // 先頭ライブラリの値は、ライブラリを指定しない経路 (表示や既定) の代表として置く
            if (p_ライブラリ番号 == 0)
            {
                this.Set_推定Phredオフセット(p_オフセット);
            }
        }

        /// <summary>
        /// クオリティカットオフ
        /// </summary>
        public int A_クオリティカットオフ { get; set; } = Consts.クオリティカットオフの既定値;

        /// <summary>
        /// k-mer カウント時にメモリ上へ保持するカウントの総量 (バイト)
        /// </summary>
        /// <remarks>
        /// メモリとディスク I/O のトレードオフを環境に合わせて調整するためのもの<br/>
        /// 増やすとフラッシュ回数が減って I/O が軽くなり、減らすとメモリが軽くなる
        /// </remarks>
        public long A_メモリ予算バイト数 { get; private set; } = Consts.メモリ予算の既定値;

        /// <summary>
        /// メモリ量の指定
        /// </summary>
        /// <remarks>
        /// "2G" / "512M" / "1024" (接尾辞なしは MB) を受け付ける
        /// </remarks>
        public string A_メモリ予算
        {
            get => Util.Get_表示用メモリサイズ(this.A_メモリ予算バイト数);
            set => this.A_メモリ予算バイト数 = Util.V_変換_メモリサイズ(value);
        }

        /// <summary>
        /// 期待インサートサイズ
        /// </summary>
        /// <remarks>
        /// CLI で明示指定されなかった場合は null のままとし、scaffolding 実行時にマップ済みペアから標本推定を試みる<br/>
        /// (自動推定できた値はこのプロパティには反映せず、Scaffolder 側で別途保持する<br/>
        /// CLI 指定値と自動推定値を区別するため)
        /// </remarks>
        public int? A_インサートサイズ { get; set; } = null;

        /// <summary>
        /// ヘルプモードか
        /// </summary>
        public bool A_Isヘルプモード { get; set; } = false;

        /// <summary>
        /// バージョンだけ表示して終わるか
        /// </summary>
        public bool A_Isバージョンモード { get; set; } = false;

        /// <summary>
        /// 進行状況メッセージの言語
        /// </summary>
        public 言語 A_言語 { get; set; } = 言語.日本語;

        /// <summary>
        /// 画面へ出す量
        /// </summary>
        /// <remarks>
        /// ファイルへの記録はこれに関わらず全量を残すので、静かにしても後から原因を追う手掛かりは失われない
        /// </remarks>
        public ログ水準 A_ログ水準 { get; set; } = ログ水準.標準;

        /// <summary>
        /// 曖昧塩基を許容するか
        /// </summary>
        public bool A_Is曖昧塩基許容 { get; set; } = false;

        /// <summary>
        /// エラー訂正するか
        /// </summary>
        public bool A_Isエラー訂正 { get; set; } = false;

        /// <summary>
        /// ペアエンドのオーバーラップ解析 (アダプタ除去 + 相互訂正) をエラー訂正・アセンブリの前に行うか
        /// </summary>
        public bool A_Is前処理 { get; set; } = false;

        /// <summary>
        /// 複数の k でアセンブリし、リファレンス無しの評価で最良のものを選ぶか
        /// </summary>
        /// <remarks>
        /// 最適な k はゲノムの反復構造で決まり、リードからは事前に分からないため、精度を求めるなら試すしかない<br/>
        /// 実行時間と引き換えになるので既定は false
        /// </remarks>
        public bool A_Isマルチk { get; set; } = false;

        /// <summary>
        /// multi-k で、前段の k の配列を次の k へ引き継ぐか
        /// </summary>
        /// <remarks>
        /// 引き継ぐのは配列だけで、繋ぐという決定は引き継がない
        /// </remarks>
        public bool A_Is引き継ぎ { get; set; } = true;

        /// <summary>
        /// multi-k で、各 k の信頼できる k-mer 集合の中でペアを橋渡しして合成リード (SuperRead) を作り、次の k への引き継ぎに加えるか
        /// </summary>
        public bool A_IsSuperRead作成 { get; set; } = false;

        /// <summary>
        /// 短い反復の解決 (V_解決_短い反復) で、対応付けを確定させる前に r-mer (アセンブリの k とは独立の短い長さ) による接合点の検証を課すか
        /// </summary>
        /// <remarks>
        /// 生リードの追加走査が 1 回 k 毎に要る (既定は false)
        /// </remarks>
        public bool A_Is反復rMer検証 { get; set; } = false;

        /// <summary>
        /// GapFiller が埋められなかった scaffold のギャップを、その両端に実際にマップされた局所リードだけを使う局所アセンブリ (LocalAssembler) で埋めるか
        /// </summary>
        /// <remarks>
        /// AssemblyMerger (-mg) の安全な代替
        /// </remarks>
        public bool A_Is局所アセンブリ { get; set; } = false;

        /// <summary>
        /// バブル除去・反復解決後の unitig グラフを GFA1 形式でも出力するか
        /// </summary>
        /// <remarks>
        /// 決められない分岐がなぜそこで打ち切られたかを、Bandage 等のビューアで直接確認できるようにする
        /// </remarks>
        public bool A_IsGFA出力 { get; set; } = false;

        /// <summary>
        /// multi-k の結果を統合するか
        /// </summary>
        /// <remarks>
        /// 既定は false<br/>
        /// 同じリードから作ったアセンブリは同じ反復配列で同じ誤りをするため、統合しても新しい情報がほとんど入らず、誤アセンブリだけが持ち込まれる
        /// </remarks>
        public bool A_Isマージ { get; set; } = false;

        /// <summary>
        /// 最終成果物に元リードを貼り直し、多数決で置換を直すか
        /// </summary>
        /// <remarks>
        /// リードを 1 回余分に走査するぶん時間がかかるため既定は false
        /// </remarks>
        public bool A_Isポリッシュ { get; set; } = false;

        /// <summary>
        /// 環状に閉じたと判定した配列について、その閉じ目を跨ぐリードが実在するかを確かめるか
        /// </summary>
        /// <remarks>
        /// 完全長を名乗るには必須の検査だが、リードの追加走査が要るため既定は false
        /// </remarks>
        public bool A_Is環状閉鎖検証 { get; set; } = false;

        /// <summary>
        /// カットオフで落ちた k-mer のうち、リードの中で信頼できる k-mer に挟まれているものを救い上げるか
        /// </summary>
        public bool A_Is救済kmer使用 { get; set; } = false;

        /// <summary>
        /// コピー数推定に使う単一コピー深度基準
        /// </summary>
        /// <remarks>
        /// cutoff の選択とは独立に保持し、既定値は従来どおり spectrum にする
        /// </remarks>
        public Tsumiki.Models.UnitigBuilding.コピー数基準の出所 A_コピー数基準の出所 { get; set; } = Tsumiki.Models.UnitigBuilding.コピー数基準の出所.Spectrum;

        /// <summary>
        /// 低カバレッジ unitig 端をトリミングするか
        /// </summary>
        /// <remarks>
        /// E2 用の切替で、false は tip 除去と cutoff を変えず端トリミングだけを止める
        /// </remarks>
        public bool A_Is低カバレッジ端トリミング { get; set; } = true;

        /// <summary>
        /// 一時ディレクトリに残っている前回の成果を再利用して途中から続けるか
        /// </summary>
        /// <remarks>
        /// 同じ入力・同じオプションで作り終えた k だけを飛ばす
        /// </remarks>
        public bool A_Is再開 { get; set; } = false;

        /// <summary>
        /// 一時ディレクトリ
        /// </summary>
        public string A_一時ディレクトリ { get; set; } = Consts.一時ディレクトリの既定値;

        /// <summary>
        /// 実行後に一時ディレクトリを消すか
        /// </summary>
        /// <remarks>
        /// k ごとの成果物が入っており後から見比べたくなるため、既定では残す
        /// </remarks>
        public bool A_Is一時ディレクトリ削除 { get; set; } = false;

        /// <summary>
        /// 並列に使うスレッド数
        /// </summary>
        public int A_スレッド数
        {
            get => this._スレッド数;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of thread count a positive integer");
                }
                this._スレッド数 = value;
            }
        }

        /// <summary>
        /// ペアの支持で結合を確定させる優勢の閾値
        /// </summary>
        public decimal A_ペア結合閾値
        {
            get => this._ペア結合閾値;
            set
            {
                if (value is <= 0M or > 1M)
                {
                    throw new ArgumentException("Please make the value of pair unite threshold a ratio between 0 (exclusive) and 1");
                }
                this._ペア結合閾値 = value;
            }
        }

        /// <summary>
        /// 結合を確定させるために要求するペアの支持数
        /// </summary>
        public ulong A_ペア支持数閾値
        {
            get => this._ペア支持数閾値;
            set
            {
                if (value <= 0UL)
                {
                    throw new ArgumentException("Please make the value of pair count threshold a positive integer");
                }
                this._ペア支持数閾値 = value;
            }
        }

        /// <summary>
        /// v0.2 仕上げ経路の設定
        /// </summary>
        /// <remarks>
        /// T01 では内部からだけ設定し、旧 CLI の引数なし動作を維持する
        /// </remarks>
        public 仕上げ設定 A_仕上げ設定 { get; set; } = new();

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 明示指定の状態を保って作業用設定を複製する
        /// </summary>
        /// <returns>元の設定と k 長一覧を共有しない複製</returns>
        /// <summary>
        /// 前処理などで作り直したパスへ差し替える
        /// </summary>
        /// <param name="p_群">ライブラリごとのリードの組</param>
        public void Set_ライブラリ群(IEnumerable<(string A_リード1, string A_リード2)> p_群)
        {
            // ペアの有無ごと保つ。リード 2 だけを詰めて持つと、シングルが混ざったときに対応がずれる
            this._差し替えたライブラリ群 = [.. p_群];
        }

        /// <summary>
        /// 実行時引数からライブラリの一覧を組み立てる
        /// </summary>
        /// <returns></returns>
        private List<(string A_リード1, string A_リード2)> Get_入力からのライブラリ群()
        {
            List<(string A_リード1, string A_リード2)> l_群 = this._リード2のパス群.Count > 0
                ? [.. this._リード1のパス群.Select((x, i) => (x, i < this._リード2のパス群.Count ? this._リード2のパス群[i] : string.Empty))]
                : [.. this._リード1のパス群.Select(x => (x, string.Empty))];
            l_群.AddRange(this._シングルのパス群.Select(x => (x, string.Empty)));
            return l_群;
        }

        /// <summary>
        /// カンマ区切りのパス指定を 1 本ずつに分ける
        /// </summary>
        /// <param name="p_指定">カンマ区切りのパス</param>
        /// <returns></returns>
        private static List<string> Get_パス群(string? p_指定)
        {
            return string.IsNullOrWhiteSpace(p_指定)
                ? []
                : [.. p_指定.Split(',').Select(x => x.Trim()).Where(x => x.Length > 0)];
        }

        public Parameters Get_複製()
        {
            var l_複製 = (Parameters)this.MemberwiseClone();
            l_複製._k長一覧 = [.. this._k長一覧];
            l_複製._リード1のパス群 = [.. this._リード1のパス群];
            l_複製._Phredオフセット群 = [.. this._Phredオフセット群];
            l_複製._ライブラリのリード長 = [.. this._ライブラリのリード長];
            l_複製._リード2のパス群 = [.. this._リード2のパス群];
            l_複製._シングルのパス群 = [.. this._シングルのパス群];
            l_複製._差し替えたライブラリ群 = this._差し替えたライブラリ群 is null ? null : [.. this._差し替えたライブラリ群];
            l_複製.A_仕上げ設定 = new()
            {
                A_Is有効 = this.A_仕上げ設定.A_Is有効,
                A_schemaバージョン = this.A_仕上げ設定.A_schemaバージョン,
                A_設定revision = this.A_仕上げ設定.A_設定revision,
            };
            return l_複製;
        }

        /// <summary>
        /// k の一覧を設定する
        /// </summary>
        /// <param name="p_k長一覧">設定する k の一覧</param>
        public void Set_k長一覧(IEnumerable<int> p_k長一覧)
        {
            var l_一覧 = p_k長一覧.Distinct().OrderBy(x => x).ToList();
            if (l_一覧.Count == 0)
            {
                throw new ArgumentException("Please set at least one kmer length");
            }
            foreach (var l_k長 in l_一覧)
            {
                if (l_k長 <= 0)
                {
                    throw new ArgumentException("Please make the value of kmer a positive integer");
                }
            }
            this._k長一覧 = l_一覧;
            this.A_k長 = l_一覧[^1];
        }

        /// <summary>
        /// CLI の実行機能をプロファイルで設定する
        /// </summary>
        /// <param name="p_名前"></param>
        public void V_適用_実行プロファイル(string p_名前)
        {
            var l_Is標準 = p_名前 switch
            {
                "standard" => true,
                "legacy" => false,
                _ => throw new ArgumentException($"Unknown profile \"{p_名前}\": expected standard or legacy"),
            };

            this.A_Is前処理 = l_Is標準;
            this.A_Isエラー訂正 = l_Is標準;
            this.A_Isマルチk = l_Is標準;
            this.A_Is引き継ぎ = true;
            this.A_IsSuperRead作成 = l_Is標準;
            this.A_Is反復rMer検証 = l_Is標準;
            this.A_Is局所アセンブリ = l_Is標準;
            this.A_IsGFA出力 = l_Is標準;
            this.A_Isポリッシュ = l_Is標準;
            this.A_Is環状閉鎖検証 = l_Is標準;
            this.A_Is救済kmer使用 = l_Is標準;
            this.A_Isマージ = false;
            this.A_コピー数基準の出所 = l_Is標準 ? Tsumiki.Models.UnitigBuilding.コピー数基準の出所.Weighted : Tsumiki.Models.UnitigBuilding.コピー数基準の出所.Spectrum;
            this.A_Is低カバレッジ端トリミング = true;
        }

        /// <summary>
        /// 推定結果から k 長を設定する
        /// </summary>
        /// <param name="p_k長">推定して得られた k 長</param>
        /// <remarks>
        /// A_Isk長明示指定 は立てないため、「ユーザーが明示指定した」扱いにはならない
        /// </remarks>
        public void Set_推定k長(int p_k長)
        {
            var l_Is明示指定済み = this.A_Isk長明示指定;
            this.A_k長 = p_k長;
            this.A_Isk長明示指定 = l_Is明示指定済み;
        }

        /// <summary>
        /// 推定結果から k-mer カットオフを設定する
        /// </summary>
        /// <param name="p_カットオフ">推定して得られたカットオフ</param>
        /// <remarks>
        /// A_Iskmerカットオフ明示指定 は立てない
        /// </remarks>
        public void Set_推定kmerカットオフ(ulong p_カットオフ)
        {
            var l_Is明示指定済み = this.A_Iskmerカットオフ明示指定;
            this.A_kmerカットオフ = p_カットオフ;
            this.A_Iskmerカットオフ明示指定 = l_Is明示指定済み;
        }

        /// <summary>
        /// 推定結果から Phred オフセットを設定する
        /// </summary>
        /// <param name="p_オフセット">推定して得られた Phred オフセット</param>
        /// <remarks>
        /// A_IsPhred明示指定 は立てないため、「ユーザーが明示指定した」扱いにはならない
        /// </remarks>
        public void Set_推定Phredオフセット(int p_オフセット)
        {
            var l_Is明示指定済み = this.A_IsPhred明示指定;
            this.A_Phredオフセット = p_オフセット;
            this.A_IsPhred明示指定 = l_Is明示指定済み;
        }

        /// <summary>
        /// (オーバーライド) 実行時引数を表す文字列を返す
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return $"""

                ============= Parameters =============

                libraries: {string.Join(" | ", this.A_ライブラリ群.Select(x => string.IsNullOrWhiteSpace(x.A_リード2) ? x.A_リード1 + " (single-end)" : x.A_リード1 + " + " + x.A_リード2))}
                kmer: {(this.A_k長一覧.Count > 1 ? string.Join(", ", this.A_k長一覧) : this.A_k長.ToString())}
                kmer cutoff: {this.A_kmerカットオフ}
                phred: {(this._Phredオフセット群.Count > 1 ? string.Join(", ", this._Phredオフセット群) : this.A_Phredオフセット.ToString())}
                quality cutoff: {this.A_クオリティカットオフ}
                counting memory budget: {this.A_メモリ予算}
                insert size: {this.A_インサートサイズ?.ToString() ?? インサートサイズ未指定表示}
                allow ambiguous bases : {this.A_Is曖昧塩基許容}
                error correction : {this.A_Isエラー訂正}
                preprocess (adapter trim + pair correction) : {this.A_Is前処理}
                multi-k : {this.A_Isマルチk}
                carry sequence between k : {this.A_Is引き継ぎ}
                build SuperReads : {this.A_IsSuperRead作成}
                verify repeat resolution with r-mers : {this.A_Is反復rMer検証}
                local assembly for remaining gaps : {this.A_Is局所アセンブリ}
                write GFA of the unitig graph : {this.A_IsGFA出力}
                merge multi-k results : {this.A_Isマージ}
                rescue mercy k-mers : {this.A_Is救済kmer使用}
                copy-number baseline requested : {this.A_コピー数基準の出所}
                trim low-coverage graph ends : {this.A_Is低カバレッジ端トリミング}
                polish final assembly with reads : {this.A_Isポリッシュ}
                verify circular closure with reads : {this.A_Is環状閉鎖検証}
                resume from temp directory : {this.A_Is再開}
                temp directory : {this.A_一時ディレクトリ}
                delete temp directory when finished : {this.A_Is一時ディレクトリ削除}
                thread count : {this.A_スレッド数}
                message language : {Get_言語名(this.A_言語)}
                console log level : {Get_ログ水準名(this.A_ログ水準)}
                pair unite threshold : {this.A_ペア結合閾値}
                pair count threshold : {this.A_ペア支持数閾値}

                ======================================

                """;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// -log に書く綴り
        /// </summary>
        /// <param name="p_水準">綴りへ変換するログ水準</param>
        /// <returns>-log に書く綴り</returns>
        /// <remarks>
        /// 表示は CLI で指定する形に合わせる
        /// </remarks>
        private static string Get_ログ水準名(ログ水準 p_水準)
        {
            return p_水準 switch
            {
                ログ水準.最小 => Consts.ログ水準名.最小,
                ログ水準.詳細 => Consts.ログ水準名.詳細,
                _ => Consts.ログ水準名.標準,
            };
        }

        /// <summary>
        /// -lang に書く綴り
        /// </summary>
        /// <param name="p_言語">綴りへ変換する言語</param>
        /// <returns>-lang に書く綴り</returns>
        /// <remarks>
        /// 表示は CLI で指定する形に合わせる
        /// </remarks>
        private static string Get_言語名(言語 p_言語)
        {
            return p_言語 switch
            {
                言語.日本語 => Consts.言語名.日本語,
                言語.中国語 => Consts.言語名.中国語,
                _ => Consts.言語名.英語,
            };
        }

        #endregion
    }
}
