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

        public int Phred { get; set; } = 33;

        public int QualityCutoff { get; set; } = 0;

        public ulong RowBitSize = int.MaxValue;
        public string BitSize
        {
            get
            {
                double aboutSize = this.RowBitSize;
                var unit = "";
                if (aboutSize >= 1e12)
                {
                    aboutSize /= 1e12;
                    unit = "T";
                }
                else if (aboutSize >= 1e9)
                {
                    aboutSize /= 1e9;
                    unit = "G";
                }
                else if (aboutSize >= 1e6)
                {
                    aboutSize /= 1e6;
                    unit = "M";
                }
                else if (aboutSize >= 1e3)
                {
                    aboutSize /= 1e3;
                    unit = "K";
                }
                return $"{aboutSize:0.#} {unit}";
            }
            set
            {
                if (value[^1] is 'K' or 'k')
                {
                    this.RowBitSize = (ulong)(double.Parse(value[..^1]) * 1e3);
                }
                else if (value[^1] is 'M' or 'm')
                {
                    this.RowBitSize = (ulong)(double.Parse(value[..^1]) * 1e6);
                }
                else
                {
                    this.RowBitSize = value[^1] is 'G' or 'g'
                        ? (ulong)(double.Parse(value[..^1]) * 1e9)
                        : value[^1] is 'T' or 't' ? (ulong)(double.Parse(value[..^1]) * 1e12) : (ulong)double.Parse(value);
                }
            }
        }

        public int InsertSize { get; set; } = 350;

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

                ==================
                """;
        }
    }
}
