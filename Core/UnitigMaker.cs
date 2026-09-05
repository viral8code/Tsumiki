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
        // k<=64 では UInt128 をキーにして1歩ごとの文字列生成を避ける。
        // 呼び出しごとに作り直さず使い回すのは、この処理が1回の実行で
        // 数百万回呼ばれるため。
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

        /// <summary>
        /// 各開始 k-mer からの walk を並列に実行し、開始 k-mer と同じ順で結果を返す。
        ///
        /// walk はカットオフ後の読み取り専用な k-mer 集合しか触らないので互いに独立。
        /// UnitigMaker 自身は呼び出しごとにクリアする訪問済み集合を持つため、
        /// ワーカーごとに1つ用意する。
        ///
        /// 重複排除は呼び出し側が元の順序で行う。どちらの向きが先に登録されるかで
        /// 採用される表現が変わるため、ここで並列に潰すと結果が実行ごとに変わる。
        /// </summary>
        public static string[] Get_walk結果(
            TrustedKmerIndex p_kmerインデックス, IReadOnlyList<byte[]> p_開始kmer)
        {
            var l_結果 = new string[p_開始kmer.Count];
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            _ = Parallel.For(
                0,
                p_開始kmer.Count,
                new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 },
                () => new UnitigMaker(p_kmerインデックス),
                (i, _, l_構築) =>
                {
                    l_結果[i] = l_構築.Get_ユニティグ(p_開始kmer[i]).A_配列;
                    return l_構築;
                },
                _ => { });

            return l_結果;
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

                // unitig は「内部のすべての節点が入次数1かつ出次数1である極大パス」。
                // 出次数だけでなく入次数も見る必要がある。次の k-mer の入次数が
                // 2以上なら別の経路が合流しており、そこからは別の unitig が始まる。
                // 怠ると合流後の共有配列を複数の unitig が重複して持ち、
                // さらにその k-mer が曖昧としてマッピング対象から外れる。
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
