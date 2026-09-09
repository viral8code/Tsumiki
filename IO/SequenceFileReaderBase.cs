using System.IO.Compression;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTA/FASTQ 共通の下回り<br/>
    /// ファイルを開く (.gz なら透過的に展開する) 処理と
    /// 空行を読み飛ばす行読み込みは、レコードの形式によらず同じで済む
    /// </summary>
    internal abstract class SequenceFileReaderBase : IDisposable
    {
        /// <summary>
        /// ファイルパス
        /// </summary>
        public string A_ファイルパス { get; }

        /// <summary>
        /// 読み込み
        /// </summary>
        private readonly StreamReader _読み込み;

        /// <summary>
        /// バッファサイズ
        /// </summary>
        private const int バッファサイズ = 1 << 25;

        protected SequenceFileReaderBase(string p_パス)
        {
            this.A_ファイルパス = p_パス;
            var l_入力ストリーム = new FileStream(p_パス, FileMode.Open, FileAccess.Read);
            this._読み込み = Path.GetExtension(p_パス)?.ToLower() == ".gz"
                ? new StreamReader(new GZipStream(l_入力ストリーム, CompressionMode.Decompress), bufferSize: バッファサイズ)
                : new StreamReader(l_入力ストリーム, bufferSize: バッファサイズ);
        }

        public bool Get_続きがあるか()
        {
            return !this._読み込み.EndOfStream;
        }

        protected virtual string Get_次の行()
        {
            var l_行 = this._読み込み.ReadLine();
            while (string.IsNullOrWhiteSpace(l_行))
            {
                l_行 = this._読み込み.ReadLine();
            }
            return l_行;
        }

        /// <summary>
        /// 空行の読み飛ばしをしない生の 1 行読み込み<br/>
        /// EOF 検知が要る派生クラス向け
        /// </summary>
        protected string? Get_次の行_生()
        {
            return this._読み込み.ReadLine();
        }

        public void Dispose()
        {
            this._読み込み.Dispose();
        }
    }
}
