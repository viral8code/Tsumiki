using Tsumiki.Common;

namespace Tsumiki.Model
{
    internal class Parameters
    {
        private string _リード1のパス = string.Empty;
        public string A_リード1のパス
        {
            get => this._リード1のパス;
            set
            {
                if (!Path.Exists(value))
                {
                    throw new ArgumentException($"Read1's path {value} is not found");
                }
                this._リード1のパス = value;
            }
        }

        private string _リード2のパス = string.Empty;
        public string A_リード2のパス
        {
            get => this._リード2のパス;
            set
            {
                if (!Path.Exists(value))
                {
                    throw new ArgumentException($"Read2's path {value} is not found");
                }
                this._リード2のパス = value;
            }
        }

        /// <summary>
        /// -k が明示的に指定されたかどうか。指定されていない場合に限り、
        /// 実際のリード長から求めた k を自動採用する。
        /// </summary>
        public bool A_k長が明示指定されたか { get; private set; }

        private int _k長 = Consts.k長の既定値;
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
                this.A_k長が明示指定されたか = true;
            }
        }

        private List<int> _k長一覧 = [];

        /// <summary>
        /// -k にカンマ区切りで指定された k の一覧(昇順・重複なし)。未指定なら空。
        /// 2個以上あれば multi-k として扱う。
        /// </summary>
        public IReadOnlyList<int> A_k長一覧 => this._k長一覧;

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
        /// 推定結果から k 長を設定する。A_k長が明示指定されたか は立てないため、
        /// 「ユーザーが明示指定した」扱いにはならない。
        /// </summary>
        public void Set_推定k長(int p_k長)
        {
            var l_明示指定済みか = this.A_k長が明示指定されたか;
            this.A_k長 = p_k長;
            this.A_k長が明示指定されたか = l_明示指定済みか;
        }

        /// <summary>
        /// -kc が明示的に指定されたかどうか。指定されていない場合に限り、
        /// k-mer スペクトルの谷から求めたカットオフを自動採用する。
        /// </summary>
        public bool A_kmerカットオフが明示指定されたか { get; private set; }

        private ulong _kmerカットオフ = Consts.kmerカットオフの既定値;
        public ulong A_kmerカットオフ
        {
            get => this._kmerカットオフ;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of kmer cut off a positive integer");
                }
                this._kmerカットオフ = value;
                this.A_kmerカットオフが明示指定されたか = true;
            }
        }

        /// <summary>
        /// 推定結果から k-mer カットオフを設定する。
        /// A_kmerカットオフが明示指定されたか は立てない。
        /// </summary>
        public void Set_推定kmerカットオフ(ulong p_カットオフ)
        {
            var l_明示指定済みか = this.A_kmerカットオフが明示指定されたか;
            this.A_kmerカットオフ = p_カットオフ;
            this.A_kmerカットオフが明示指定されたか = l_明示指定済みか;
        }

        /// <summary>
        /// -p が明示的に指定されたかどうか。指定されていない場合に限り、
        /// FASTQ のクオリティ文字列から推定したオフセットを自動採用する。
        /// 明示指定はユーザーの判断なので、推定結果で上書きはしない。
        /// </summary>
        public bool A_Phredが明示指定されたか { get; private set; }

        private int _Phredオフセット = Consts.Phredオフセットの既定値;
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
                this.A_Phredが明示指定されたか = true;
            }
        }

        /// <summary>
        /// 推定結果から Phred オフセットを設定する。A_Phredが明示指定されたか は
        /// 立てないため、「ユーザーが明示指定した」扱いにはならない。
        /// </summary>
        public void Set_推定Phredオフセット(int p_オフセット)
        {
            var l_明示指定済みか = this.A_Phredが明示指定されたか;
            this.A_Phredオフセット = p_オフセット;
            this.A_Phredが明示指定されたか = l_明示指定済みか;
        }

        public int A_クオリティカットオフ { get; set; } = Consts.クオリティカットオフの既定値;

        /// <summary>
        /// k-mer カウント時にメモリ上へ保持するカウントの総量(バイト)。
        /// メモリとディスク I/O のトレードオフを環境に合わせて調整するためのもの。
        /// 増やすとフラッシュ回数が減って I/O が軽くなり、減らすとメモリが軽くなる。
        /// </summary>
        public long A_メモリ予算バイト数 { get; private set; } = Consts.メモリ予算の既定値;

        /// <summary>
        /// メモリ量の指定。"2G" / "512M" / "1024"(接尾辞なしは MB)を受け付ける。
        /// </summary>
        public string A_メモリ予算
        {
            get => Util.Get_表示用メモリサイズ(this.A_メモリ予算バイト数);
            set => this.A_メモリ予算バイト数 = Util.V_変換_メモリサイズ(value);
        }

        /// <summary>
        /// 期待インサートサイズ。CLI で明示指定されなかった場合は null のままとし、
        /// スキャフォールディング実行時にマップ済みペアから標本推定を試みる。
        /// (自動推定できた値はこのプロパティには反映せず、Scaffolder 側で
        ///  別途保持する。CLI 指定値と自動推定値を区別するため。)
        /// </summary>
        public int? A_インサートサイズ { get; set; } = null;

        public bool A_ヘルプモードか { get; set; } = false;

        /// <summary>バージョンだけ表示して終わるか。</summary>
        public bool A_バージョンモードか { get; set; } = false;

        /// <summary>進行状況メッセージの言語。</summary>
        public 言語 A_言語 { get; set; } = 言語.日本語;

        /// <summary>
        /// 画面へ出す量。ファイルへの記録はこれに関わらず全量を残すので、
        /// 静かにしても後から原因を追う手掛かりは失われない。
        /// </summary>
        public ログ水準 A_ログ水準 { get; set; } = ログ水準.標準;

        /// <summary>-log に書く綴り。表示は CLI で指定する形に合わせる。</summary>
        private static string Get_ログ水準名(ログ水準 p_水準)
        {
            return p_水準 switch
            {
                ログ水準.最小 => Consts.ログ水準名.最小,
                ログ水準.詳細 => Consts.ログ水準名.詳細,
                _ => Consts.ログ水準名.標準,
            };
        }

        /// <summary>-lang に書く綴り。表示は CLI で指定する形に合わせる。</summary>
        private static string Get_言語名(言語 p_言語)
        {
            return p_言語 switch
            {
                言語.日本語 => Consts.言語名.日本語,
                言語.中国語 => Consts.言語名.中国語,
                _ => Consts.言語名.英語,
            };
        }

        public bool A_曖昧塩基を許容するか { get; set; } = false;

        public bool A_エラー訂正するか { get; set; } = false;

        /// <summary>
        /// ペアエンドのオーバーラップ解析(アダプタ除去 + 相互訂正)を
        /// エラー訂正・アセンブリの前に行うか。
        /// </summary>
        public bool A_前処理するか { get; set; } = false;

        /// <summary>
        /// 複数の k でアセンブリし、リファレンス無しの評価で最良のものを選ぶか。
        /// 最適な k はゲノムの反復構造で決まり、リードからは事前に分からないため、
        /// 精度を求めるなら試すしかない。実行時間と引き換えになるので既定は false。
        /// </summary>
        public bool A_マルチkか { get; set; } = false;

        /// <summary>
        /// multi-k で、前段の k の配列を次の k へ引き継ぐか。
        /// 引き継ぐのは配列だけで、繋ぐという決定は引き継がない。
        /// </summary>
        public bool A_引き継ぐか { get; set; } = true;

        /// <summary>
        /// multi-k で、各 k の信頼できる k-mer 集合の中でペアを橋渡しして
        /// 合成リード(SuperRead)を作り、次の k への引き継ぎに加えるか。
        /// </summary>
        public bool A_SuperReadを作るか { get; set; } = false;

        /// <summary>
        /// 短い反復の解決(V_解決_短い反復)で、対応付けを確定させる前に
        /// r-mer(アセンブリの k とは独立の短い長さ)による接合点の検証を
        /// 課すか。生リードの追加走査が1回k毎に要る(既定は false)。
        /// </summary>
        public bool A_反復をrMerで検証するか { get; set; } = false;

        /// <summary>
        /// GapFiller が埋められなかったスキャフォールドのギャップを、
        /// その両端に実際にマップされた局所リードだけを使う局所アセンブリ
        /// (LocalAssembler)で埋めるか。AssemblyMerger(-mg)の安全な代替。
        /// </summary>
        public bool A_局所アセンブリするか { get; set; } = false;

        /// <summary>
        /// バブル除去・反復解決後の unitig グラフを GFA1 形式でも出力するか。
        /// 決められない分岐がなぜそこで打ち切られたかを、Bandage 等の
        /// ビューアで直接確認できるようにする。
        /// </summary>
        public bool A_GFAを出力するか { get; set; } = false;

        /// <summary>
        /// multi-k の結果を統合するか。既定は false。
        /// 同じリードから作ったアセンブリは同じ反復配列で同じ誤りをするため、
        /// 統合しても新しい情報がほとんど入らず、誤アセンブリだけが持ち込まれる。
        /// </summary>
        public bool A_マージするか { get; set; } = false;

        /// <summary>
        /// 最終成果物に元リードを貼り直し、多数決で置換を直すか。
        /// リードを1回余分に走査するぶん時間がかかるため既定は false。
        /// </summary>
        public bool A_ポリッシュするか { get; set; } = false;

        /// <summary>
        /// 環状に閉じたと判定した配列について、その閉じ目を跨ぐリードが
        /// 実在するかを確かめるか。完全長を名乗るには必須の検査だが、
        /// リードの追加走査が要るため既定は false。
        /// </summary>
        public bool A_環状閉鎖を検証するか { get; set; } = false;

        /// <summary>
        /// カットオフで落ちた k-mer のうち、リードの中で信頼できる k-mer に
        /// 挟まれているものを救い上げるか。
        /// </summary>
        public bool A_救済kmerを使うか { get; set; } = false;

        /// <summary>
        /// 一時ディレクトリに残っている前回の成果を再利用して途中から続けるか。
        /// 同じ入力・同じオプションで作り終えた k だけを飛ばす。
        /// </summary>
        public bool A_再開するか { get; set; } = false;

        public string A_一時ディレクトリ { get; set; } = Consts.一時ディレクトリの既定値;

        /// <summary>
        /// 実行後に一時ディレクトリを消すか。k ごとの成果物が入っており
        /// 後から見比べたくなるため、既定では残す。
        /// </summary>
        public bool A_一時ディレクトリを削除するか { get; set; } = false;

        private int _スレッド数 = Environment.ProcessorCount;
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

        private decimal _ペア結合閾値 = Consts.ペア結合閾値の既定値;
        public decimal A_ペア結合閾値
        {
            get => this._ペア結合閾値;
            set
            {
                if (value is <= 0 or > 1)
                {
                    throw new ArgumentException("Please make the value of pair unite threshold a ratio between 0 (exclusive) and 1");
                }
                this._ペア結合閾値 = value;
            }
        }

        private ulong _ペア支持数閾値 = Consts.ペア支持数閾値の既定値;
        public ulong A_ペア支持数閾値
        {
            get => this._ペア支持数閾値;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of pair count threshold a positive integer");
                }
                this._ペア支持数閾値 = value;
            }
        }

        public override string ToString()
        {
            return $"""

                ============= Parameters =============

                read1: {this.A_リード1のパス}
                read2: {this.A_リード2のパス}
                kmer: {(this.A_k長一覧.Count > 1 ? string.Join(", ", this.A_k長一覧) : this.A_k長.ToString())}
                kmer cutoff: {this.A_kmerカットオフ}
                phred: {this.A_Phredオフセット}
                quality cutoff: {this.A_クオリティカットオフ}
                counting memory budget: {this.A_メモリ予算}
                insert size: {this.A_インサートサイズ?.ToString() ?? Consts.インサートサイズ未指定表示}
                allow ambiguous bases : {this.A_曖昧塩基を許容するか}
                error correction : {this.A_エラー訂正するか}
                preprocess (adapter trim + pair correction) : {this.A_前処理するか}
                multi-k : {this.A_マルチkか}
                carry sequence between k : {this.A_引き継ぐか}
                build SuperReads : {this.A_SuperReadを作るか}
                verify repeat resolution with r-mers : {this.A_反復をrMerで検証するか}
                local assembly for remaining gaps : {this.A_局所アセンブリするか}
                write GFA of the unitig graph : {this.A_GFAを出力するか}
                merge multi-k results : {this.A_マージするか}
                rescue mercy k-mers : {this.A_救済kmerを使うか}
                polish final assembly with reads : {this.A_ポリッシュするか}
                verify circular closure with reads : {this.A_環状閉鎖を検証するか}
                resume from temp directory : {this.A_再開するか}
                temp directory : {this.A_一時ディレクトリ}
                delete temp directory when finished : {this.A_一時ディレクトリを削除するか}
                thread count : {this.A_スレッド数}
                message language : {Get_言語名(this.A_言語)}
                console log level : {Get_ログ水準名(this.A_ログ水準)}
                pair unite threshold : {this.A_ペア結合閾値}
                pair count threshold : {this.A_ペア支持数閾値}

                ======================================

                """;
        }
    }
}
