namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// グラフ成果物を識別する revision
    /// </summary>
    /// <param name="A_番号"></param>
    /// <param name="A_内容ハッシュ"></param>
    internal readonly record struct グラフrevision(long A_番号, string A_内容ハッシュ)
    {
        #region 公開メソッド

        /// <summary>
        /// revision の値を検証する
        /// </summary>
        public void V_検証()
        {
            if (this.A_番号 < 0L)
            {
                throw new ArgumentOutOfRangeException(nameof(this.A_番号), "Graph revision must not be negative");
            }

            if (string.IsNullOrWhiteSpace(this.A_内容ハッシュ) || this.A_内容ハッシュ.Length != 64 || !this.A_内容ハッシュ.All(Uri.IsHexDigit))
            {
                throw new ArgumentException("Graph revision hash must be a SHA-256 hexadecimal string", nameof(this.A_内容ハッシュ));
            }
        }

        #endregion
    }
}
