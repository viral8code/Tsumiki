using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    internal class UnitigMaker(TrustedKmerIndex bloomFilter)
    {
        private readonly TrustedKmerIndex set = bloomFilter;

        public Unitig MakeUnitig(Span<byte> bytes)
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            List<byte> list = [.. bytes[..^1]];
            var now = string.Join(string.Empty, bytes.ToArray().Select(Util.ByteToBaseString));
            HashSet<string> visited = [];
            while (visited.Add(now))
            {
                list.Add((byte)Util.GetNucleotideIDs(now[^1])[0]);
                string? nextKmer = null;
                byte nextBase = 0;
                list.Add(0);
                for (byte i = 1; i <= 4; i++)
                {
                    list[^1] = i;
                    if (this.set.Contains(CollectionsMarshal.AsSpan(list)[(list.Count - kmerLength)..]))
                    {
                        if (nextKmer != null)
                        {
                            nextKmer = null;
                            break;
                        }
                        nextKmer = now[1..] + Util.ByteToBaseString(i);
                        nextBase = i;
                    }
                }

                if (nextKmer != null)
                {
                    // unitig の定義は「内部のすべての節点が入次数1かつ出次数1である
                    // 極大パス」。この walk は出次数の条件(候補がちょうど1つ)しか
                    // 見ておらず、入次数の条件を見ていなかった。次のk-merの入次数が
                    // 2以上(=別の経路もそこへ合流している)なら、そこからは別の
                    // unitig が始まるべきで、ここで止めなければならない。
                    // これを怠ると、合流後の共有配列を複数のunitigが重複して持ち、
                    // 実データで全unitigのk-mer延べ数が実際のゲノム内容の
                    // 1.43倍にまで膨らんでいた(=同じ配列の重複出力)。さらに
                    // ContigMaker 側では、複数unitigに跨って現れるk-merが
                    // AmbiguousKmer として大量にマッピング対象から除外され
                    // (実データで1100万件超)、隣接情報が大きく損なわれていた。
                    list[^1] = nextBase;
                    if (this.set.CountInEdges(CollectionsMarshal.AsSpan(list)[(list.Count - kmerLength)..]) != 1)
                    {
                        nextKmer = null;
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
