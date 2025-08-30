namespace Tsumiki.IO
{
    internal class FastaWriter(string fileName) : IDisposable
    {
        private readonly StreamWriter _writer = new(fileName);

        public void Write(object seqID, string sequence)
        {
            this._writer.Write(">");
            this._writer.WriteLine(seqID);
            this._writer.WriteLine(sequence);
        }

        public void Dispose()
        {
            this._writer?.Dispose();
        }
    }
}
