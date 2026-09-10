using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// リファレンス無しでアセンブリの良さを測る
    /// </summary>
    /// <remarks>
    /// k が違えば k-mer 集合の意味も変わるため、比較には固定した
    /// アンカー k の集合を物差しとして使う
    /// </remarks>
    internal static class AssemblyScorer
    {
        #region 定数

        /// <summary>
        /// 評価に含める配列の最小長
        /// </summary>
        /// <remarks>
        /// abyss-fac の既定と同じ 500 bp<br/>
        /// これより短い断片は、そこに配列が入っていても下流で使いようがない<br/>
        /// k-mer の集計にも掛けるのが要点で、掛けないと「短い破片を大量に
        /// 出しただけ」のアセンブリが完全性で得をする<br/>
        /// 実データでは
        /// k=21 の 3,421 本のうち 500 bp 以上は 322 本しかなく、残りが
        /// 完全性を底上げして低い k を有利にしていた
        /// </remarks>
        private const int 評価に含める最小長 = 500;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_FASTAパス のアセンブリを、アンカー k-mer 集合に対して評価する
        /// </summary>
        /// <param name="p_FASTAパス"></param>
        /// <param name="p_アンカーインデックス"></param>
        /// <param name="p_アンカーk長"></param>
        /// <param name="p_単一コピー基準値"></param>
        /// <param name="p_推定ゲノムサイズ"></param>
        /// <remarks>
        /// アンカー k が 64 を超える場合 (2 bit パックが UInt128 に収まらない) は
        /// 評価できないため null を返す
        /// </remarks>
        public static アセンブリ評価? Get_評価(string p_FASTAパス, TrustedKmerIndex p_アンカーインデックス, int p_アンカーk長, double p_単一コピー基準値, long p_推定ゲノムサイズ)
        {
            if (p_アンカーk長 > 64 || p_単一コピー基準値 <= 0)
            {
                return null;
            }

            var l_観測 = Get_出現回数(p_FASTAパス, p_アンカーk長, out var l_長さ一覧, out var l_総延長, out var l_環状本数, out var l_環状延長);

            var l_期待延べ数 = 0L;
            var l_欠損延べ数 = 0L;
            var l_過剰延べ数 = 0L;

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

            var l_統計対象 = l_長さ一覧.Where(x => x >= 評価に含める最小長).ToList();

            return new アセンブリ評価(
                A_期待延べ数: l_期待延べ数,
                A_欠損延べ数: l_欠損延べ数,
                A_過剰延べ数: l_過剰延べ数,
                A_総延長: l_統計対象.Sum(),
                A_本数: l_統計対象.Count,
                A_NG50: Get_NG50(l_統計対象, p_推定ゲノムサイズ, l_総延長),
                A_環状本数: l_環状本数,
                A_環状化率: p_推定ゲノムサイズ > 0 ? (double)l_環状延長 / p_推定ゲノムサイズ : 0);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// アセンブリ中に各アンカー k-mer が何回現れるかを数える
        /// </summary>
        /// <param name="p_FASTAパス"></param>
        /// <param name="p_アンカーk長"></param>
        /// <param name="p_長さ一覧"></param>
        /// <param name="p_総延長"></param>
        /// <param name="p_環状本数"></param>
        /// <param name="p_環状延長"></param>
        /// <remarks>
        /// 逆相補は同一視する<br/>
        /// 併せて、環状に閉じた配列 (名前に環状の目印が付いたもの) の
        /// 本数と総延長も集計する<br/>
        /// 短すぎる配列は数えない<br/>
        /// 環状の目印は既にそれより長い閉路にしか
        /// 付かないため、環状の集計がこの足切りで漏れることはない
        /// </remarks>
        private static Dictionary<UInt128, int> Get_出現回数(string p_FASTAパス, int p_アンカーk長, out List<int> p_長さ一覧, out long p_総延長, out int p_環状本数, out long p_環状延長)
        {
            Dictionary<UInt128, int> l_観測 = [];
            p_長さ一覧 = [];
            p_総延長 = 0;
            p_環状本数 = 0;
            p_環状延長 = 0;

            using var l_読み込み = new FastaReader(p_FASTAパス);
            while (l_読み込み.Get_続きがあるか())
            {
                var l_エントリ = l_読み込み.Get_次の配列();
                var l_配列 = l_エントリ.A_配列;
                if (l_配列.Length < 評価に含める最小長)
                {
                    continue;
                }
                p_長さ一覧.Add(l_配列.Length);
                p_総延長 += l_配列.Length;

                if (l_エントリ.A_ID.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase))
                {
                    p_環状本数++;
                    p_環状延長 += l_配列.Length;
                }

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
        /// NG50
        /// </summary>
        /// <param name="p_長さ一覧"></param>
        /// <param name="p_推定ゲノムサイズ"></param>
        /// <param name="p_総延長"></param>
        /// <remarks>
        /// 素の N50 は自分の総延長を分母にするため、配列を落として
        /// 短くなったアセンブリほど有利になり k を跨いだ比較に使えない<br/>
        /// ゲノムサイズが分からない場合は総延長で代用する (=素の N50)
        /// </remarks>
        private static long Get_NG50(List<int> p_長さ一覧, long p_推定ゲノムサイズ, long p_総延長)
        {
            var l_分母 = p_推定ゲノムサイズ > 0 ? p_推定ゲノムサイズ : p_総延長;
            if (l_分母 <= 0)
            {
                return 0;
            }

            var l_累積 = 0L;
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

        #endregion
    }
}
