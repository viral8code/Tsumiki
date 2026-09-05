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

        public bool A_曖昧塩基を許容するか { get; set; } = false;

        public bool A_エラー訂正するか { get; set; } = false;

        /// <summary>
        /// 複数の k でアセンブリし、リファレンス無しの評価で最良のものを選ぶか。
        /// 最適な k はゲノムの反復構造で決まり、リードからは事前に分からないため、
        /// 精度を求めるなら試すしかない。実行時間と引き換えになるので既定は false。
        /// </summary>
        public bool A_マルチkか { get; set; } = false;

        public string A_一時ディレクトリ { get; set; } = Consts.一時ディレクトリの既定値;

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
                kmer: {this.A_k長}
                kmer cutoff: {this.A_kmerカットオフ}
                phred: {this.A_Phredオフセット}
                quality cutoff: {this.A_クオリティカットオフ}
                counting memory budget: {this.A_メモリ予算}
                insert size: {this.A_インサートサイズ?.ToString() ?? Consts.インサートサイズ未指定表示}
                allow ambiguous bases : {this.A_曖昧塩基を許容するか}
                error correction : {this.A_エラー訂正するか}
                multi-k : {this.A_マルチkか}
                temp directory : {this.A_一時ディレクトリ}
                thread count : {this.A_スレッド数}
                pair unite threshold : {this.A_ペア結合閾値}
                pair count threshold : {this.A_ペア支持数閾値}

                ======================================

                """;
        }
    }
}
