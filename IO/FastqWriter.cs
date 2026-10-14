using System.IO.Compression;
using Tsumiki.Utilities;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTQ を 1 リードずつ書き出す
    /// </summary>
    /// <param name="p_ファイル名"></param>
    /// <param name="p_Is圧縮">真なら gzip (最速の圧縮) で書く</param>
    internal class FastqWriter(string p_ファイル名, bool p_Is圧縮 = false) : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 書き込み (閉じると、下に挟んだ gzip も閉じる)
        /// </summary>
        private readonly StreamWriter _書き込み = new(Get_書込ストリーム(p_ファイル名, p_Is圧縮));

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 1 リードを書き出す
        /// </summary>
        /// <param name="p_ID"></param>
        /// <param name="p_配列"></param>
        /// <param name="p_クオリティ"></param>
        public void V_書き込み(string p_ID, string p_配列, string p_クオリティ)
        {
            this._書き込み.WriteLine(p_ID);
            this._書き込み.WriteLine(p_配列);
            this._書き込み.WriteLine('+');
            this._書き込み.WriteLine(p_クオリティ);
        }

        /// <summary>
        /// 保持している資源を解放する
        /// </summary>
        public void Dispose()
        {
            this._書き込み?.Dispose();
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 書き込むストリームを開く (圧縮するときは gzip で包む)
        /// </summary>
        /// <param name="p_ファイル名">書き込み先</param>
        /// <param name="p_Is圧縮">真なら gzip で包む</param>
        /// <returns>書き込みストリーム</returns>
        private static Stream Get_書込ストリーム(string p_ファイル名, bool p_Is圧縮)
        {
            var l_書き込み = 中間データ置き場.Get_書込ストリーム(p_ファイル名);
            return p_Is圧縮 ? new GZipStream(l_書き込み, CompressionLevel.Fastest) : l_書き込み;
        }

        #endregion
    }
}
