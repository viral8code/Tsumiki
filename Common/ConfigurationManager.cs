using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Common
{
    internal static class ConfigurationManager
    {
        public static Parameters Arguments { get; set; } = new();

        public static TrustedKmerIndex TrustedKmerIndex { get; set; } = null!;
    }
}
