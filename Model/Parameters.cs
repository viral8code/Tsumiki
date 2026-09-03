using Tsumiki.Common;

namespace Tsumiki.Model
{
    internal class Parameters
    {
        private string _path1 = string.Empty;
        public string ReadPath1
        {
            get => this._path1;
            set
            {
                if (!Path.Exists(value))
                {
                    throw new ArgumentException($"Read1's path {value} is not found");
                }
                this._path1 = value;
            }
        }

        private string _path2 = string.Empty;
        public string ReadPath2
        {
            get => this._path2;
            set
            {
                if (!Path.Exists(value))
                {
                    throw new ArgumentException($"Read2's path {value} is not found");
                }
                this._path2 = value;
            }
        }

        private int _kmer = Consts.DefaultKmerValue;
        public int Kmer
        {
            get => this._kmer;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of kmer a positive integer");
                }
                this._kmer = value;
            }
        }

        private ulong _kmerCutoff = Consts.DefaultKmerCutoffValue;
        public ulong KmerCutoff
        {
            get => this._kmerCutoff;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of kmer cut off a positive integer");
                }
                this._kmerCutoff = value;
            }
        }

        /// <summary>
        /// -p が明示的に指定されたかどうか。指定されていない場合に限り、
        /// FASTQ のクオリティ文字列から推定したオフセットを自動採用する。
        /// 明示指定はユーザーの判断なので、推定結果で上書きはしない。
        /// </summary>
        public bool IsPhredExplicitlySet { get; private set; }

        private int _phred = Consts.DefaultPhredValue;
        public int Phred
        {
            get => this._phred;
            set
            {
                if (!Consts.AllowedPhredValue.Contains(value))
                {
                    throw new ArgumentException($"Phred value is must {string.Join(" or ", Consts.AllowedPhredValue)}");
                }
                this._phred = value;
                this.IsPhredExplicitlySet = true;
            }
        }

        /// <summary>
        /// 推定結果から Phred オフセットを設定する。IsPhredExplicitlySet は立てないため、
        /// 「ユーザーが明示指定した」扱いにはならない。
        /// </summary>
        public void SetInferredPhred(int value)
        {
            var wasExplicit = this.IsPhredExplicitlySet;
            this.Phred = value;
            this.IsPhredExplicitlySet = wasExplicit;
        }

        public int QualityCutoff { get; set; } = Consts.DefaultQualityCutoffValue;

        public ulong RowBitSize = int.MaxValue;
        public string BitSize
        {
            get
            {
                double aboutSize = this.RowBitSize;
                var unit = "";
                if (aboutSize >= 8e12)
                {
                    aboutSize /= 8e12;
                    unit = "T";
                }
                else if (aboutSize >= 8e9)
                {
                    aboutSize /= 8e9;
                    unit = "G";
                }
                else if (aboutSize >= 8e6)
                {
                    aboutSize /= 8e6;
                    unit = "M";
                }
                else if (aboutSize >= 8e3)
                {
                    aboutSize /= 8e3;
                    unit = "K";
                }
                return $"{aboutSize:0.#} {unit}";
            }

            set => this.RowBitSize =
                      value[^1] is 'K' or 'k'
                    ? (ulong)(double.Parse(value[..^1]) * 8e3)
                    : value[^1] is 'M' or 'm'
                    ? (ulong)(double.Parse(value[..^1]) * 8e6)
                    : value[^1] is 'G' or 'g'
                    ? (ulong)(double.Parse(value[..^1]) * 8e9)
                    : value[^1] is 'T' or 't' ? (ulong)(double.Parse(value[..^1]) * 8e12) : (ulong)double.Parse(value);
        }

        /// <summary>
        /// 期待挿入サイズ。CLI で明示指定されなかった場合は null のままとし、
        /// スキャフォールディング実行時にマップ済みペアからサンプリング推定を試みる。
        /// (自動推定できた値はこのプロパティには反映せず、Scaffolder 側で
        ///  別途保持する。CLI 指定値と自動推定値を区別するため。)
        /// </summary>
        public int? InsertSize { get; set; } = null;

        public bool IsHelpMode { get; set; } = false;

        public bool AllowAmbiguousBases { get; set; } = false;

        public bool EnableErrorCorrection { get; set; } = false;

        public string TempDirectory { get; set; } = Consts.DefaultTempFolder;

        private int _threadCount = Environment.ProcessorCount;
        public int ThreadCount
        {
            get => this._threadCount;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of thread count a positive integer");
                }
                this._threadCount = value;
            }
        }

        private decimal _pairUniteThreshold = Consts.DefaultPairUniteThreshold;
        public decimal PairUniteThreshold
        {
            get => this._pairUniteThreshold;
            set
            {
                if (value is <= 0 or > 1)
                {
                    throw new ArgumentException("Please make the value of pair unite threshold a ratio between 0 (exclusive) and 1");
                }
                this._pairUniteThreshold = value;
            }
        }

        private ulong _pairCountThreshold = Consts.DefaultPairCountThreshold;
        public ulong PairCountThreshold
        {
            get => this._pairCountThreshold;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of pair count threshold a positive integer");
                }
                this._pairCountThreshold = value;
            }
        }

        public override string ToString()
        {
            return $"""

                ============= Parameters =============

                read1: {this.ReadPath1}
                read2: {this.ReadPath2}
                kmer: {this.Kmer}
                kmer cutoff: {this.KmerCutoff}
                phred: {this.Phred}
                quality cutoff: {this.QualityCutoff}
                bit size: {this.BitSize}
                insert size: {this.InsertSize?.ToString() ?? Consts.NullInsertSizeText}
                allow ambiguous bases : {this.AllowAmbiguousBases}
                error correction : {this.EnableErrorCorrection}
                temp directory : {this.TempDirectory}
                thread count : {this.ThreadCount}
                pair unite threshold : {this.PairUniteThreshold}
                pair count threshold : {this.PairCountThreshold}

                ======================================

                """;
        }
    }
}