using Tsumiki.Common;
using Tsumiki.Model.Foundation;

namespace Tsumiki.IO
{
    internal class FastaReader(string p_パス) : SequenceFileReaderBase(p_パス)
    {
        /// <summary>
        /// FASTA を1回で全件読み込む。ID の先頭 '>' は取り除く。
        /// </summary>
        public static List<(string A_ID, string A_配列)> Get_全エントリ(string p_パス)
        {
            List<(string, string)> l_結果 = [];
            using var l_読み込み = new FastaReader(p_パス);
            while (l_読み込み.Get_続きがあるか())
            {
                var l_エントリ = l_読み込み.Get_次の配列();
                l_結果.Add((l_エントリ.A_ID.TrimStart('>'), l_エントリ.A_配列));
            }
            return l_結果;
        }

        public 配列エントリ Get_次の配列()
        {
            try
            {
                var l_ID = this.Get_次の行();
                var l_配列 = this.Get_次の行();

                return new 配列エントリ(l_ID, l_配列);
            }
            catch (Exception ex)
            {
                Logger.V_出力_警告(Logger.Get_メソッド名(), ex);
                throw;
            }
        }
    }
}
