namespace Tsumiki.IO
{
    internal class FastqWriter(string fileName) : IDisposable
    {
        private readonly StreamWriter _writer = new(fileName);

        /// <summary>
        /// id は先頭の "@" を含む形で渡すこと(FastqReader.NextRead().ID がそのまま使える)。
        /// </summary>
        public void Write(string id, string sequence, string quality)
        {
            this._writer.WriteLine(id);
            this._writer.WriteLine(sequence);
            this._writer.WriteLine("+");
            this._writer.WriteLine(quality);
        }

        public void Dispose()
        {
            this._writer?.Dispose();
        }
    }
}
