using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    internal class UnitigMaker(TrustedKmerIndex bloomFilter)
    {
        private readonly TrustedKmerIndex set = bloomFilter;

        // 循環検出用。walk 中に同じ k-mer へ戻ってきたら打ち切る。
        // k<=64 なら 2bit パックした UInt128 をキーにでき、1歩ごとの
        // 文字列生成(k=63 では 1 歩あたり 63 文字の string を複数回確保)を
        // 完全に避けられる。k>64 の場合のみ従来どおり文字列で判定する。
        //
        // 呼び出しごとに作り直すのではなく使い回して確保を減らす
        // (GraphSimplifier は反復のたびに全 unitig を作り直すため、
        //  MakeUnitig は1回の実行で数百万回呼ばれる)。
        private readonly HashSet<UInt128> visitedPacked = [];

        private readonly HashSet<string> visitedText = [];

        /// <summary>
        /// k-mer(塩基ID 1-4、長さ64以下)を2bit/塩基で UInt128 にパックする。
        /// 向き依存の値(逆相補への正規化はしない)。循環検出は
        /// 「同じ向きで同じ k-mer に戻ったか」で判定する必要があるため。
        /// </summary>
        private static UInt128 Pack(ReadOnlySpan<byte> kmer)
        {
            UInt128 value = 0;
            foreach (var b in kmer)
            {
                value = (value << 2) | (UInt128)(b - 1);
            }
            return value;
        }

        public Unitig MakeUnitig(Span<byte> bytes)
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            var usePackedVisited = kmerLength <= 64;

            this.visitedPacked.Clear();
            this.visitedText.Clear();

            // list の末尾 kmerLength 塩基が常に「現在の k-mer」になる。
            List<byte> list = [.. bytes];

            while (true)
            {
                var current = CollectionsMarshal.AsSpan(list)[(list.Count - kmerLength)..];

                var isNew = usePackedVisited
                    ? this.visitedPacked.Add(Pack(current))
                    : this.visitedText.Add(string.Join(string.Empty, current.ToArray().Select(Util.ByteToBaseString)));
                if (!isNew)
                {
                    // 循環。従来実装は「この k-mer の最後の1塩基を付ける前」に
                    // 打ち切っていたため、同じ配列になるよう1塩基取り除く。
                    list.RemoveAt(list.Count - 1);
                    break;
                }

                // 次の1塩基を決める。候補がちょうど1つ(出次数1)でなければ
                // ここが unitig の終端。
                list.Add(0);
                byte nextBase = 0;
                var candidateCount = 0;
                for (byte i = Consts.NucleotideID.A; i <= Consts.NucleotideID.T; i++)
                {
                    list[^1] = i;
                    if (this.set.Contains(CollectionsMarshal.AsSpan(list)[(list.Count - kmerLength)..]))
                    {
                        candidateCount++;
                        if (candidateCount > 1)
                        {
                            break;
                        }
                        nextBase = i;
                    }
                }

                if (candidateCount != 1)
                {
                    list.RemoveAt(list.Count - 1);
                    break;
                }

                // unitig の定義は「内部のすべての節点が入次数1かつ出次数1である
                // 極大パス」。出次数の条件(候補がちょうど1つ)だけでなく、
                // 入次数の条件も見る必要がある。次のk-merの入次数が2以上
                // (=別の経路もそこへ合流している)なら、そこからは別の
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
                    list.RemoveAt(list.Count - 1);
                    break;
                }
            }

            return new Unitig(list.GetHashCode(), string.Join(string.Empty, list.Select(Util.ByteToBaseString)));
        }
    }
}
