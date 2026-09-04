namespace Tsumiki.IO
{
    internal class FastqWriter(string p_ファイル名) : IDisposable
    {
        private readonly StreamWriter _書き込み = new(p_ファイル名);

        /// <summary>
        /// p_ID は先頭の "@" を含む形で渡すこと
        /// (FastqReader.Get_次のリード().A_ID がそのまま使える)。
        /// </summary>
        public void V_書き込み(string p_ID, string p_配列, string p_クオリティ)
        {
            this._書き込み.WriteLine(p_ID);
            this._書き込み.WriteLine(p_配列);
            this._書き込み.WriteLine("+");
            this._書き込み.WriteLine(p_クオリティ);
        }

        public void Dispose()
        {
            this._書き込み?.Dispose();
        }
    }
}
