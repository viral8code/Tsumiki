using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Tsumiki.IO
{
    internal class FastaWriter(string fileName) : IDisposable
    {
        private readonly StreamWriter _writer = new(fileName);

        public void Write(string seqID, string sequence)
        {
            _writer.Write(">");
            _writer.WriteLine(seqID);
            _writer.WriteLine(sequence);
        }

        public void Dispose()
        {
            _writer?.Dispose();
        }
    }
}
