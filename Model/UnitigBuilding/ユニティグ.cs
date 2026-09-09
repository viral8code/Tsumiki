using System.Text;

namespace Tsumiki.Model.UnitigBuilding
{
    /// <summary>
    /// 分岐を持たない 1 本の配列
    /// </summary>
    /// <param name="p_ID">ユニティグ ID</param>
    /// <param name="p_配列">ユニティグの配列</param>
    internal class ユニティグ(object p_ID, string p_配列)
    {
        /// <summary>
        /// ユニティグ ID
        /// </summary>
        public readonly string A_ID = p_ID?.ToString() ?? string.Empty;

        /// <summary>
        /// 配列
        /// </summary>
        public readonly string A_配列 = p_配列;

        /// <summary>
        /// (オーバーライド) ユニティグを FASTA の 2 行として返す
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            var l_文字列 = new StringBuilder();
            _ = l_文字列.AppendLine($"ID: {this.A_ID}")
                .AppendLine($"Seq: {this.A_配列}");
            return l_文字列.ToString();
        }
    }
}
