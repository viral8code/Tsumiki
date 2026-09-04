using System.Text;

namespace Tsumiki.Model
{
    internal class ユニティグ(object p_ID, string p_配列)
    {
        public readonly string A_ID = p_ID?.ToString() ?? string.Empty;
        public readonly string A_配列 = p_配列;

        public override string ToString()
        {
            var l_文字列 = new StringBuilder();
            _ = l_文字列.AppendLine($"ID: {this.A_ID}")
                .AppendLine($"Seq: {this.A_配列}");
            return l_文字列.ToString();
        }
    }
}
