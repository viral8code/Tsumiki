using System.IO.Compression;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTA/FASTQ 共通の下回り
    /// </summary>
    /// <remarks>
    /// ファイルを開く (.gz なら透過的に展開する) 処理と
    /// 空行を読み飛ばす行読み込みは、レコードの形式によらず同じで済む
    /// </remarks>
    internal abstract class SequenceFileReaderBase : IDisposable
    {
        #region 定数

        /// <summary>
        /// バッファサイズ
        /// </summary>
        private const int バッファサイズ = 1 << 25;

        #endregion

        #region 内部変数

        /// <summary>
        /// 読み込み
        /// </summary>
        private readonly StreamReader _読み込み;

        #endregion

        #region プロパティ

        /// <summary>
        /// ファイルパス
        /// </summary>
        public string A_ファイルパス { get; }

        #endregion

        #region コンストラクタ

        /// <summary>
        /// ファイルを開く
        /// </summary>
        /// <param name="p_パス"></param>
        /// <remarks>
        /// 拡張子が .gz なら透過的に展開する
        /// </remarks>
        protected SequenceFileReaderBase(string p_パス)
        {
            this.A_ファイルパス = p_パス;
            var l_入力ストリーム = new FileStream(p_パス, FileMode.Open, FileAccess.Read);
            this._読み込み = Path.GetExtension(p_パス)?.ToLower() == ".gz"
                ? new StreamReader(new GZipStream(l_入力ストリーム, CompressionMode.Decompress), bufferSize: バッファサイズ)
                : new StreamReader(l_入力ストリーム, bufferSize: バッファサイズ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// まだ読める行があるか
        /// </summary>
        /// <returns>続きがあれば true</returns>
        public bool Get_続きがあるか()
        {
            return !this._読み込み.EndOfStream;
        }

        /// <summary>
        /// 次の 1 行を読み込んで返す
        /// </summary>
        /// <returns>読み込んだ行</returns>
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
        /// 空行の読み飛ばしをしない生の 1 行読み込み
        /// </summary>
        /// <remarks>
        /// EOF 検知が要る派生クラス向け
        /// </remarks>
        protected string? Get_次の行_生()
        {
            return this._読み込み.ReadLine();
        }

        /// <summary>
        /// 保持している資源を解放する
        /// </summary>
        public void Dispose()
        {
            this._読み込み.Dispose();
        }

        #endregion
    }
}
