using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Util
{
    internal class FastqReader : IDisposable
    {
        private readonly string filePath;

        private readonly StreamReader reader;

        public FastqReader(string path)
        {
            filePath = path;
            if (File.Exists(path))
            {
                reader = new StreamReader(path);
            }
            else
            {
                throw new ArgumentException(path + " is not found");
            }
        }

        public string GetFilePath()
        {
            return filePath;
        }

        public string nextData()
        {
            var dataLine = reader.ReadLine();
            while (string.IsNullOrWhiteSpace(dataLine))
            {
                dataLine = reader.ReadLine();
            }
            return dataLine;
        }

        public ReadData nextRead()
        {
            try
            {
                var id = nextData();
                var read = nextData();
                _ = nextData();
                var quality = nextData();

                return new ReadData()
                {
                    ID = id,
                    Read = read,
                    Quality = quality,
                };
            }
            catch (Exception ex)
            {
                Logger.PrintInnerException(Logger.GetMethodName(), ex);
                throw;
            }
        }

        public void Dispose()
        {
            reader.Dispose();
        }
    }
}
