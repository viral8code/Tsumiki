using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Common
{
    internal static class ConfigurationManager
    {
        public static Parameters A_実行時引数 { get; set; } = new();

        public static TrustedKmerIndex A_kmerインデックス { get; set; } = null!;
    }
}
