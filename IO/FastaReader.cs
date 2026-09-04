using System.IO.Compression;
using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.IO
{
    internal class FastaReader : IDisposable
    {
        public string A_ファイルパス { get; private set; }

        private readonly StreamReader _読み込み;

        private const int バッファサイズ = 1 << 25;

        public FastaReader(string p_パス)
        {
            this.A_ファイルパス = p_パス;
            var l_入力ストリーム = new FileStream(p_パス, FileMode.Open, FileAccess.Read);
            if (Path.GetExtension(p_パス)?.ToLower() == ".gz")
            {
                var l_展開ストリーム = new GZipStream(l_入力ストリーム, CompressionMode.Decompress);
                this._読み込み = new(l_展開ストリーム, bufferSize: バッファサイズ);
            }
            else
            {
                this._読み込み = new(l_入力ストリーム, bufferSize: バッファサイズ);
            }
        }

        public bool Get_続きがあるか()
        {
            return !this._読み込み.EndOfStream;
        }

        private string Get_次の行()
        {
            var l_行 = this._読み込み.ReadLine();
            while (string.IsNullOrWhiteSpace(l_行))
            {
                l_行 = this._読み込み.ReadLine();
            }
            return l_行;
        }

        public 配列エントリ Get_次の配列()
        {
            try
            {
                var l_ID = this.Get_次の行();
                var l_配列 = this.Get_次の行();

                return new 配列エントリ(l_ID, l_配列);
            }
            catch (Exception ex)
            {
                Logger.V_出力_警告(Logger.Get_メソッド名(), ex);
                throw;
            }
        }

        public void Dispose()
        {
            this._読み込み.Dispose();
        }
    }
}
