using System.Text;
using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    internal class UnitigMaker(CountingBloomFilter bloomFilter)
    {
        private readonly CountingBloomFilter set = bloomFilter;

        public Unitig MakeUnitig(string kmer)
        {
            var sb = new StringBuilder(kmer[..^1]);
            var now = kmer;
            HashSet<string> visited = [];
            while (visited.Add(now))
            {
                _ = sb.Append(now[^1]);
                string? nextKmer = null;
                for (byte i = 1; i <= 4; i++)
                {
                    var next = now[1..] + Util.ByteToBase(i);
                    if (this.set.Contains(next))
                    {
                        if (nextKmer != null)
                        {
                            break;
                        }
                        nextKmer = next;
                    }
                }
                if (nextKmer == null)
                {
                    break;
                }
                now = nextKmer;
            }
            return new Unitig(kmer.GetHashCode(), sb.ToString());
        }
    }
}
