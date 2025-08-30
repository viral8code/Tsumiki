using System.Text;

namespace Tsumiki.Model
{
    internal class Unitig(object id, string seq)
    {
        public readonly string Id = id?.ToString() ?? string.Empty;
        public readonly string Sequence = seq;

        public override string ToString()
        {
            var sb = new StringBuilder();
            _ = sb.AppendLine($"ID: {this.Id}")
                .AppendLine($"Seq: {Sequence}");
            return sb.ToString();
        }
    }
}
