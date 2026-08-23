using System.IO.Compression;
using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.IO
{
    internal class FastqReader : IDisposable
    {
        public string FilePath { get; private set; }

        private readonly StreamReader reader;

        private const int BufferedSize = 1 << 25;

        public FastqReader(string path)
        {
            this.FilePath = path;
            var inputFileStream = new FileStream(path, FileMode.Open, FileAccess.Read);
            if (Path.GetExtension(path)?.ToLower() == ".gz")
            {
                var decompressionStream = new GZipStream(inputFileStream, CompressionMode.Decompress);
                this.reader = new(decompressionStream, bufferSize: BufferedSize);
            }
            else
            {
                this.reader = new(inputFileStream, bufferSize: BufferedSize);
            }
        }

        public bool HasNext()
        {
            return !this.reader.EndOfStream;
        }

        private string NextData()
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
                    Read = Util.ToByteList(read),
                    RowRead = read,
                    Quality = quality,
                };
            }
            catch (Exception ex)
            {
                Logger.PrintWarning(Logger.GetMethodName(), ex);
                throw;
            }
        }

        /// <summary>
        /// 曖昧塩基を無視する経路向けの軽量版。Read(List&lt;byte[]&gt;)の代わりに
        /// SimpleRead(byte[])のみを構築する。LoadReadFileToBloomFilterIgnoreAmbiguity から使用する。
        /// </summary>
        public ReadData NextReadSimple()
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
                    SimpleRead = Util.ToSimpleByteArray(read),
                    RowRead = read,
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