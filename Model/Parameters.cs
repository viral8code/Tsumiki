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

        private int _kmer = 31;
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

        private ulong _kmerCutoff = 2;
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

        private int _phred = 33;
        public int Phred
        {
            get => this._phred;
            set
            {
                if (value is not 33 and not 64)
                {
                    throw new ArgumentException("Phred value is must 33 or 64");
                }
                this._phred = value;
            }
        }

        public int QualityCutoff { get; set; } = 0;

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

        public int? InsertSize { get; set; } = null;

        public bool IsHelpMode { get; set; } = false;

        public bool AllowAmbiguousBases { get; set; } = false;

        public override string ToString()
        {
            return $"""
                === Parameters ===

                read1: {this.ReadPath1}
                read2: {this.ReadPath2}
                kmer: {this.Kmer}
                kmer cutoff: {this.KmerCutoff}
                phred: {this.Phred}
                quality cutoff: {this.QualityCutoff}
                bit size: {this.BitSize}
                insert size: {this.InsertSize?.ToString() ?? "unspecified"}
                allow ambiguous bases : {this.AllowAmbiguousBases}

                ==================
                """;
        }
    }
}
