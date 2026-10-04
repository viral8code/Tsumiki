using System.Text;
using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTA を 1 配列ずつ読み込む
    /// </summary>
    /// <param name="p_パス"></param>
    internal class FastaReader(string p_パス) : SequenceFileReaderBase(p_パス)
    {
        #region 定数

        /// <summary>
        /// 見出し行の最初の文字
        /// </summary>
        private const char C_見出しの印 = '>';

        #endregion

        #region 公開メソッド

        /// <summary>
        /// FASTA を 1 回で全件読み込む
        /// </summary>
        /// <param name="p_パス"></param>
        /// <returns></returns>
        public static List<(string A_ID, string A_配列)> Get_全エントリ(string p_パス)
        {
            List<(string, string)> l_結果 = [];
            using var l_読み込み = new FastaReader(p_パス);
            while (l_読み込み.Has続き())
            {
                var l_エントリ = l_読み込み.Get_次の配列();
                l_結果.Add((l_エントリ.A_ID.TrimStart('>'), l_エントリ.A_配列));
            }

            return l_結果;
        }

        /// <summary>
        /// 次の 1 配列を読み込んで返す (配列が複数行に折り返されていれば、次の見出し行の手前までをつなぐ)
        /// </summary>
        /// <returns>読み込んだ配列</returns>
        public 配列エントリ Get_次の配列()
        {
            try
            {
                var l_ID = this.Get_次の行();
                var l_配列 = this.Is見出しの前か終わり(C_見出しの印) ? string.Empty : this.Get_次の行();
                if (this.Is見出しの前か終わり(C_見出しの印))
                {
                    return new 配列エントリ(l_ID, l_配列);
                }

                var l_つないだ配列 = new StringBuilder(l_配列);
                while (!this.Is見出しの前か終わり(C_見出しの印))
                {
                    _ = l_つないだ配列.Append(this.Get_次の行_生()?.Trim());
                }

                return new 配列エントリ(l_ID, l_つないだ配列.ToString());
            }
            catch (Exception l_例外)
            {
                Logger.V_出力_警告(Logger.Get_メソッド名(), l_例外);
                throw;
            }
        }

        #endregion
    }
}
