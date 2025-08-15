using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Tsumiki.Model
{
    internal class Parameters
    {
        public string ReadPath1 { get; set; } = string.Empty;

        public string ReadPath2 { get; set; } = string.Empty;

        public int Kmer { get; set; } = 31;

        public int KmerCutoff { get; set; } = 3;

        public int Phred { get; set; } = 33;

        public int QualityCutoff { get; set; } = 0;
    }
}
