using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// リファレンス無しでアセンブリの良さを測る。
    /// k が違えば k-mer 集合の意味も変わるため、比較には固定した
    /// アンカー k の集合を物差しとして使う。
    /// </summary>
    internal static class AssemblyScorer
    {
        /// <summary>
        /// 連続性の統計に含める配列の最小長。k-mer の集計側には掛けない。
        /// 短い配列に入っていても「出せている」ことに変わりはないため。
        /// </summary>
        private const int 連続性統計の最小長 = 500;

        /// <summary>
        /// p_FASTAパス のアセンブリを、アンカー k-mer 集合に対して評価する。
        /// アンカー k が 64 を超える場合(2bit パックが UInt128 に収まらない)は
        /// 評価できないため null を返す。
        /// </summary>
        public static アセンブリ評価? Get_評価(
            string p_FASTAパス,
            TrustedKmerIndex p_アンカーインデックス,
            int p_アンカーk長,
            double p_単一コピー基準値,
            long p_推定ゲノムサイズ)
        {
            if (p_アンカーk長 > 64 || p_単一コピー基準値 <= 0)
            {
                return null;
            }

            var l_観測 = Get_出現回数(p_FASTAパス, p_アンカーk長, out var l_長さ一覧, out var l_総延長);

            long l_期待延べ数 = 0;
            long l_欠損延べ数 = 0;
            long l_過剰延べ数 = 0;

            foreach (var l_kmer in p_アンカーインデックス.Get_信頼kmer一覧())
            {
                var l_正規形 = KmerPacking.Get_正規化パック(l_kmer);
                var l_カバレッジ = p_アンカーインデックス.Get_カバレッジ(l_kmer);
                var l_期待コピー数 = Math.Max(1, (int)Math.Round(l_カバレッジ / p_単一コピー基準値));
                var l_出現数 = l_観測.GetValueOrDefault(l_正規形);

                l_期待延べ数 += l_期待コピー数;
                if (l_出現数 < l_期待コピー数)
                {
                    l_欠損延べ数 += l_期待コピー数 - l_出現数;
                }
                else
                {
                    l_過剰延べ数 += l_出現数 - l_期待コピー数;
                }
            }

            var l_統計対象 = l_長さ一覧.Where(x => x >= 連続性統計の最小長).ToList();

            return new アセンブリ評価(
                A_期待延べ数: l_期待延べ数,
                A_欠損延べ数: l_欠損延べ数,
                A_過剰延べ数: l_過剰延べ数,
                A_総延長: l_統計対象.Sum(),
                A_本数: l_統計対象.Count,
                A_NG50: Get_NG50(l_統計対象, p_推定ゲノムサイズ, l_総延長));
        }

        /// <summary>
        /// アセンブリ中に各アンカー k-mer が何回現れるかを数える。逆相補は同一視する。
        /// </summary>
        private static Dictionary<UInt128, int> Get_出現回数(
            string p_FASTAパス, int p_アンカーk長, out List<int> p_長さ一覧, out long p_総延長)
        {
            Dictionary<UInt128, int> l_観測 = [];
            p_長さ一覧 = [];
            p_総延長 = 0;

            using var l_読み込み = new FastaReader(p_FASTAパス);
            while (l_読み込み.Get_続きがあるか())
            {
                var l_配列 = l_読み込み.Get_次の配列().A_配列;
                p_長さ一覧.Add(l_配列.Length);
                p_総延長 += l_配列.Length;

                for (var i = 0; i + p_アンカーk長 <= l_配列.Length; i++)
                {
                    if (KmerPacking.Get_正規化パック(l_配列, i, p_アンカーk長, out var l_正規形))
                    {
                        l_観測[l_正規形] = l_観測.GetValueOrDefault(l_正規形) + 1;
                    }
                }
            }
            return l_観測;
        }

        /// <summary>
        /// NG50。素の N50 は自分の総延長を分母にするため、配列を落として
        /// 短くなったアセンブリほど有利になり k を跨いだ比較に使えない。
        /// ゲノムサイズが分からない場合は総延長で代用する(=素の N50)。
        /// </summary>
        private static long Get_NG50(List<int> p_長さ一覧, long p_推定ゲノムサイズ, long p_総延長)
        {
            var l_分母 = p_推定ゲノムサイズ > 0 ? p_推定ゲノムサイズ : p_総延長;
            if (l_分母 <= 0)
            {
                return 0;
            }

            long l_累積 = 0;
            foreach (var l_長さ in p_長さ一覧.OrderByDescending(x => x))
            {
                l_累積 += l_長さ;
                if (l_累積 * 2 >= l_分母)
                {
                    return l_長さ;
                }
            }
            return 0;
        }
    }
}
