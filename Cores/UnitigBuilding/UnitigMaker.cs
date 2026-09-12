using System.Runtime.InteropServices;
using Tsumiki.Commons;
using Tsumiki.Models.UnitigBuilding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// 信頼できる k-mer 集合の中を walk して、分岐を持たない 1 本の配列を作る
    /// </summary>
    /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
    internal class UnitigMaker(TrustedKmerIndex p_kmerインデックス)
    {
        #region 内部変数

        /// <summary>
        /// kmer インデックス
        /// </summary>
        private readonly TrustedKmerIndex _kmerインデックス = p_kmerインデックス;

        /// <summary>
        /// 訪問済み パック
        /// </summary>
        private readonly HashSet<UInt128> _訪問済み_パック = [];

        /// <summary>
        /// 訪問済み 文字列
        /// </summary>
        private readonly HashSet<string> _訪問済み_文字列 = [];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 各開始 k-mer からの walk を並列に実行し、開始 k-mer と同じ順で結果を返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_開始kmer"></param>
        /// <returns></returns>
        public static string[] Get_walk結果(TrustedKmerIndex p_kmerインデックス, IReadOnlyList<byte[]> p_開始kmer)
        {
            var l_結果 = new string[p_開始kmer.Count];
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            _ = Parallel.For(0, p_開始kmer.Count, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, () => new 走査状態(p_kmerインデックス), (i, _, l_状態) =>
                {
                    l_結果[i] = l_状態.Get_配列(p_開始kmer[i]);
                    return l_状態;
                }, _ => { });

            return l_結果;
        }

        /// <summary>
        /// 開始 k-mer から walk して unitig を返す
        /// </summary>
        /// <param name="p_開始kmer">walk を始める k-mer</param>
        /// <returns>組み上がった unitig</returns>
        public ユニティグ Get_ユニティグ(Span<byte> p_開始kmer)
        {
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            var l_Isパック経路使用 = l_k長 <= 64;

            this._訪問済み_パック.Clear();
            this._訪問済み_文字列.Clear();

            // 末尾 k 長 塩基が常に「現在の k-mer」になる
            List<byte> l_配列 = [.. p_開始kmer];

            while (true)
            {
                var l_現在のkmer = CollectionsMarshal.AsSpan(l_配列)[(l_配列.Count - l_k長)..];

                var l_Is未訪問 = l_Isパック経路使用
                    ? this._訪問済み_パック.Add(TryGet_パック(l_現在のkmer))
                    : this._訪問済み_文字列.Add(string.Join(string.Empty, l_現在のkmer.ToArray().Select(Util.V_変換_塩基文字)));
                if (!l_Is未訪問)
                {
                    // 循環
                    // 従来実装は「この k-mer の最後の 1 塩基を付ける前」に
                    // 打ち切っていたため、同じ配列になるよう 1 塩基取り除く
                    l_配列.RemoveAt(l_配列.Count - 1);
                    break;
                }

                // 次の 1 塩基を決める
                // 候補がちょうど 1 つ (出次数 1) でなければ
                // ここが unitig の終端
                l_配列.Add(0);
                byte l_次の塩基 = 0;
                var l_候補数 = 0;
                for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
                {
                    l_配列[^1] = i;
                    if (this._kmerインデックス.Haskmer(CollectionsMarshal.AsSpan(l_配列)[(l_配列.Count - l_k長)..]))
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

                // unitig は「内部のすべての節点が入次数 1 かつ出次数 1 である極大パス」
                // 出次数だけでなく入次数も見る必要がある
                // 次の k-mer の入次数が
                // 2 以上なら別の経路が合流しており、そこからは別の unitig が始まる
                // 怠ると合流後の共有配列を複数の unitig が重複して持ち、
                // さらにその k-mer が曖昧としてマッピング対象から外れる
                l_配列[^1] = l_次の塩基;
                if (this._kmerインデックス.Get_入次数(CollectionsMarshal.AsSpan(l_配列)[(l_配列.Count - l_k長)..]) != 1)
                {
                    l_配列.RemoveAt(l_配列.Count - 1);
                    break;
                }
            }

            return new ユニティグ(l_配列.GetHashCode(), string.Join(string.Empty, l_配列.Select(Util.V_変換_塩基文字)));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k-mer (塩基 ID 1-4、長さ 64 以下) を 2 bit/塩基で UInt128 にパックする
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <remarks>
        /// 向き依存の値 (逆相補への正規化はしない) <br/>
        /// 循環検出は「同じ向きで同じ k-mer に戻ったか」で判定する必要があるため
        /// </remarks>
        /// <returns></returns>
        private static UInt128 TryGet_パック(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_値 = 0;
            foreach (var l_塩基ID in p_kmer)
            {
                l_値 = (l_値 << 2) | (UInt128)(l_塩基ID - 1);
            }
            return l_値;
        }

        #endregion
    }
}
