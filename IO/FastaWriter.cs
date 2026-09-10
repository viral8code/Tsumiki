namespace Tsumiki.IO
{
    /// <summary>
    /// FASTA を 1 配列ずつ書き出す
    /// </summary>
    internal class FastaWriter(string p_ファイル名) : IDisposable
    {
        #region 内部変数

        /// <summary>
        /// 書き込み
        /// </summary>
        private readonly StreamWriter _書き込み = new(p_ファイル名);

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 1 配列を書き出す
        /// </summary>
        /// <param name="p_配列ID">配列 ID</param>
        /// <param name="p_配列">配列</param>
        public void V_書き込み(object p_配列ID, string p_配列)
        {
            this._書き込み.Write(">");
            this._書き込み.WriteLine(p_配列ID);
            this._書き込み.WriteLine(p_配列);
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
