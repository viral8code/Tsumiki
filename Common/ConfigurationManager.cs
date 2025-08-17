using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Common
{
    internal static class ConfigurationManager
    {
        public static Parameters Arguments { get; set; } = new();

        public static CountingBloomFilter BloomFilter { get; set; } = new(0);
    }
}
