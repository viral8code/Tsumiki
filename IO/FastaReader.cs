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
        #region 公開メソッド

        /// <summary>
        /// FASTA を 1 回で全件読み込む
        /// </summary>
        /// <param name="p_パス"></param>
        /// <remarks>
        /// ID の先頭 '>' は取り除く
        /// </remarks>
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
        /// 次の 1 配列を読み込んで返す
        /// </summary>
        /// <returns>読み込んだ配列</returns>
        public 配列エントリ Get_次の配列()
        {
            try
            {
                var l_ID = this.Get_次の行();
                var l_配列 = this.Get_次の行();

                return new 配列エントリ(l_ID, l_配列);
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
