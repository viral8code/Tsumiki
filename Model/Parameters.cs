using System.IO;

namespace Tsumiki.Model
{
    internal class Parameters
    {
        private string _path1 = string.Empty;
        public string ReadPath1
        {
            get
            {
                return _path1;
            }
            set
            {
                if (!Path.Exists(value))
                {
                    throw new ArgumentException($"Read1's path {value} is not found");
                }
                _path1 = value;
            }
        }

        private string _path2 = string.Empty;
        public string ReadPath2
        {
            get
            {
                return _path2;
            }
            set
            {
                if (!Path.Exists(value))
                {
                    throw new ArgumentException($"Read2's path {value} is not found");
                }
                _path2 = value;
            }
        }

        private int _kmer = 31;
        public int Kmer
        {
            get
            {
                return _kmer;
            }
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of kmer a positive integer");
                }
                _kmer = value;
            }
        }

        private int _kmerCutoff = 2;
        public int KmerCutoff
        {
            get
            {
                return _kmerCutoff;
            }
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Please make the value of kmer cut off a positive integer");
                }
                _kmerCutoff = value;
            }
        }

        public int Phred { get; set; } = 33;

        public int QualityCutoff { get; set; } = 0;

        private ulong _bitSize = int.MaxValue;
        public string BitSize
        {
            get
            {
                double aboutSize = _bitSize;
                string unit = "";
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
                    _bitSize = (ulong)(double.Parse(value[..^1]) * 1e3);
                }
                else if (value[^1] is 'M' or 'm')
                {
                    _bitSize = (ulong)(double.Parse(value[..^1]) * 1e6);
                }
                else if (value[^1] is 'G' or 'g')
                {
                    _bitSize = (ulong)(double.Parse(value[..^1]) * 1e9);
                }
                else if (value[^1] is 'T' or 't')
                {
                    _bitSize = (ulong)(double.Parse(value[..^1]) * 1e12);
                }
                else
                {
                    _bitSize = (ulong)double.Parse(value);
                }
            }
        }

        public override string ToString()
        {
            return $"""
                === Parameters ===

                read1: {ReadPath1}
                read2: {ReadPath2}
                kmer: {Kmer}
                kmer cutoff: {KmerCutoff}
                phred: {Phred}
                quality cutoff: {QualityCutoff}
                bit size: {BitSize}

                ==================
                """;
        }
    }
}
