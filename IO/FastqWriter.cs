namespace Tsumiki.IO
{
    /// <summary>
    /// FASTQ を 1 リードずつ書き出す
    /// </summary>
    /// <param name="p_ファイル名"></param>
    internal class FastqWriter(string p_ファイル名) : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 書き込み
        /// </summary>
        private readonly StreamWriter _書き込み = new(p_ファイル名);

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 1 リードを書き出す
        /// </summary>
        /// <param name="p_ID"></param>
        /// <param name="p_配列"></param>
        /// <param name="p_クオリティ"></param>
        /// <remarks>
        /// p_ID は先頭の "@" を含む形で渡すこと<br/>
        /// (FastqReader.Get_次のリード () .A_ID がそのまま使える)
        /// </remarks>
        public void V_書き込み(string p_ID, string p_配列, string p_クオリティ)
        {
            this._書き込み.WriteLine(p_ID);
            this._書き込み.WriteLine(p_配列);
            this._書き込み.WriteLine("+");
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
    }
}
