using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Utility
{
    internal class FastqReader(string path) : IDisposable
    {
        public string FilePath { get; private set; } = path;

        private readonly StreamReader reader = new(path);

        public bool HasNext()
        {
            return !this.reader.EndOfStream;
        }

        public string NextData()
        {
            var dataLine = this.reader.ReadLine();
            while (string.IsNullOrWhiteSpace(dataLine))
            {
                dataLine = this.reader.ReadLine();
            }
            return dataLine;
        }

        public ReadData NextRead()
        {
            try
            {
                var id = this.NextData();
                var read = this.NextData();
                _ = this.NextData();
                var quality = this.NextData();

                return new ReadData()
                {
                    ID = id,
                    Read = read,
                    Quality = quality,
                };
            }
            catch (Exception ex)
            {
                Logger.PrintWarning(Logger.GetMethodName(), ex);
                throw;
            }
        }

        public void Dispose()
        {
            this.reader.Dispose();
        }
    }
}
