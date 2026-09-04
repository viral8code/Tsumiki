namespace Tsumiki.IO
{
    internal class FastaWriter(string p_ファイル名) : IDisposable
    {
        private readonly StreamWriter _書き込み = new(p_ファイル名);

        public void V_書き込み(object p_配列ID, string p_配列)
        {
            this._書き込み.Write(">");
            this._書き込み.WriteLine(p_配列ID);
            this._書き込み.WriteLine(p_配列);
        }

        public void Dispose()
        {
            this._書き込み?.Dispose();
        }
    }
}
