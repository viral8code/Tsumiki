using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    internal class UnitigMaker(TrustedKmerIndex p_kmerインデックス)
    {
        private readonly TrustedKmerIndex _kmerインデックス = p_kmerインデックス;

        // 循環検出用。walk 中に同じ k-mer へ戻ってきたら打ち切る。
        // k<=64 なら 2bit パックした UInt128 をキーにでき、1歩ごとの
        // 文字列生成(k=63 では 1 歩あたり 63 文字の string を複数回確保)を
        // 完全に避けられる。k>64 の場合のみ従来どおり文字列で判定する。
        //
        // 呼び出しごとに作り直すのではなく使い回して確保を減らす
        // (GraphSimplifier は反復のたびに全 unitig を作り直すため、
        //  この処理は1回の実行で数百万回呼ばれる)。
        private readonly HashSet<UInt128> _訪問済み_パック = [];

        private readonly HashSet<string> _訪問済み_文字列 = [];

        /// <summary>
        /// k-mer(塩基ID 1-4、長さ64以下)を2bit/塩基で UInt128 にパックする。
        /// 向き依存の値(逆相補への正規化はしない)。循環検出は
        /// 「同じ向きで同じ k-mer に戻ったか」で判定する必要があるため。
        /// </summary>
        private static UInt128 Get_パック(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_値 = 0;
            foreach (var l_塩基ID in p_kmer)
            {
                l_値 = (l_値 << 2) | (UInt128)(l_塩基ID - 1);
            }
            return l_値;
        }

        public ユニティグ Get_ユニティグ(Span<byte> p_開始kmer)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            var l_パック経路を使うか = l_k長 <= 64;

            this._訪問済み_パック.Clear();
            this._訪問済み_文字列.Clear();

            // 末尾 k長 塩基が常に「現在の k-mer」になる。
            List<byte> l_配列 = [.. p_開始kmer];

            while (true)
            {
                var l_現在のkmer = CollectionsMarshal.AsSpan(l_配列)[(l_配列.Count - l_k長)..];

                var l_未訪問か = l_パック経路を使うか
                    ? this._訪問済み_パック.Add(Get_パック(l_現在のkmer))
                    : this._訪問済み_文字列.Add(string.Join(string.Empty, l_現在のkmer.ToArray().Select(Util.V_変換_塩基文字)));
                if (!l_未訪問か)
                {
                    // 循環。従来実装は「この k-mer の最後の1塩基を付ける前」に
                    // 打ち切っていたため、同じ配列になるよう1塩基取り除く。
                    l_配列.RemoveAt(l_配列.Count - 1);
                    break;
                }

                // 次の1塩基を決める。候補がちょうど1つ(出次数1)でなければ
                // ここが unitig の終端。
                l_配列.Add(0);
                byte l_次の塩基 = 0;
                var l_候補数 = 0;
                for (byte i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
                {
                    l_配列[^1] = i;
                    if (this._kmerインデックス.Get_含まれるか(CollectionsMarshal.AsSpan(l_配列)[(l_配列.Count - l_k長)..]))
                    {
                        l_候補数++;
                        if (l_候補数 > 1)
                        {
                            break;
                        }
                        l_次の塩基 = i;
                    }
                }

                if (l_候補数 != 1)
                {
                    l_配列.RemoveAt(l_配列.Count - 1);
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
                // 曖昧kmerとして大量にマッピング対象から除外され
                // (実データで1100万件超)、隣接情報が大きく損なわれていた。
                l_配列[^1] = l_次の塩基;
                if (this._kmerインデックス.Get_入次数(CollectionsMarshal.AsSpan(l_配列)[(l_配列.Count - l_k長)..]) != 1)
                {
                    l_配列.RemoveAt(l_配列.Count - 1);
                    break;
                }
            }

            return new ユニティグ(l_配列.GetHashCode(), string.Join(string.Empty, l_配列.Select(Util.V_変換_塩基文字)));
        }
    }
}
