using System.IO.Compression;
using Tsumiki.Utilities;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTA/FASTQ 共通の下回り
    /// </summary>
    internal abstract class SequenceFileReaderBase : IDisposable
    {
        #region 定数

        /// <summary>
        /// バッファサイズ
        /// </summary>
        private const int C_バッファサイズ = 1 << 25;

        /// <summary>
        /// 圧縮された配列ファイルの拡張子
        /// </summary>
        private const string C_圧縮ファイル拡張子 = ".gz";

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
        protected SequenceFileReaderBase(string p_パス)
        {
            this.A_ファイルパス = p_パス;
            var l_入力ストリーム = 中間データ置き場.Get_読込ストリーム(p_パス);
            this._読み込み = Path.GetExtension(p_パス)?.ToLower() == C_圧縮ファイル拡張子
                ? new StreamReader(new GZipStream(l_入力ストリーム, CompressionMode.Decompress), bufferSize: C_バッファサイズ)
                : new StreamReader(l_入力ストリーム, bufferSize: C_バッファサイズ);
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// まだ読める行があるか
        /// </summary>
        /// <returns>続きがあれば true</returns>
        public bool Has続き()
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
        /// <returns></returns>
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
