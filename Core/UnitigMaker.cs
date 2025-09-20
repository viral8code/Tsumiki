using System.Runtime.InteropServices;
using System.Text;
using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    internal class UnitigMaker(CountingBloomFilter bloomFilter)
    {
        private readonly CountingBloomFilter set = bloomFilter;

        public Unitig MakeUnitig(Span<byte> bytes)
        {
            List<byte> list = [.. bytes[..^1]];
            var now = string.Join(string.Empty, bytes.ToArray().Select(Util.ByteToBaseString));
            HashSet<string> visited = [];
            while (visited.Add(now))
            {
                list.Add((byte)Util.GetNucleotideIDs(now[^1])[0]);
                string? nextKmer = null;
                list.Add(0);
                for (byte i = 1; i <= 4; i++)
                {
                    list[^1] = i;
                    if (this.set.Contains(CollectionsMarshal.AsSpan(list)[(list.Count - ConfigurationManager.Arguments.Kmer)..]))
                    {
                        if (nextKmer != null)
                        {
                            nextKmer = null;
                            break;
                        }
                        nextKmer = now[1..] + Util.ByteToBaseString(i);
                    }
                }
                list.RemoveAt(list.Count - 1);
                if (nextKmer == null)
                {
                    break;
                }
                now = nextKmer;
            }
            return new Unitig(list.GetHashCode(), string.Join(string.Empty, list.Select(Util.ByteToBaseString)));
        }
    }
}
