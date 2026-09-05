using Tsumiki.Common;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// de Bruijn グラフの簡略化。2種類のアーティファクトを除去する。
    ///
    /// 1. tip 除去: 短い行き止まりの unitig を丸ごと除去する。行き止まりは
    ///    定義上どこにも合流しないので、他の経路が使う配列を壊さない。
    ///
    /// 2. 低カバレッジ端のトリミング: 合流点を特定せず、各 unitig の両端から
    ///    カバレッジが基準値比で著しく低い k-mer が続く間だけ剥がす。
    ///
    ///    unitig 全体の平均で判定してはいけない。SNP 様の短い分岐では共有部分の
    ///    高カバレッジに平均が引きずられて検出できず、仮に検出できても unitig 全体を
    ///    除去すると合流後の共有配列まで消して別の経路を壊す。
    ///    エラー由来の分岐は合流点までの区間だけが低カバレッジなので、そこだけ剥がす。
    /// </summary>
    internal static class GraphSimplifier
    {
        /// <summary>
        /// tip 除去と低カバレッジ端のトリミングを反復し、簡略化後の unitig 開始点を返す。
        /// 除去のたびに unitig を再構築して次数とカバレッジを評価し直すため反復する。
        ///
        /// 長さ閾値は tip 除去にのみ適用する。低カバレッジ端のトリミングは
        /// 合流先の長さに依存せず判定できるため、長さによらず全件に適用する。
        /// カバレッジの基準値は長さ加重中央値を使う。単純平均や単純中央値だと
        /// 本数の多い短い断片に引きずられ、主経路の水準から外れる。
        /// </summary>
        public static List<byte[]> V_除去_tip(
            TrustedKmerIndex p_kmerインデックス,
            int p_k長,
            int? p_tip長閾値 = null,
            int p_最大反復数 = 30,
            double p_低カバレッジ比 = 0.2)
        {
            var l_tip長閾値 = p_tip長閾値 ?? (10 * p_k長);
            var l_開始kmer = p_kmerインデックス.Get_開始kmer一覧();

            for (var l_反復 = 1; l_反復 <= p_最大反復数; l_反復++)
            {
                var l_ユニティグ群 = Get_ユニティグ群(p_kmerインデックス, l_開始kmer);
                var l_基準値 = Get_長さ加重中央カバレッジ(p_kmerインデックス, l_ユニティグ群, p_k長);
                var l_低カバレッジ閾値 = l_基準値 * p_低カバレッジ比;

                var l_除去tip数 = 0;
                var l_剥がしたkmer数 = 0;
                var l_トリミングしたunitig数 = 0;
                foreach (var l_ユニティグ in l_ユニティグ群)
                {
                    var l_塩基列 = Get_塩基列(l_ユニティグ);
                    if (l_塩基列.Length < p_k長)
                    {
                        continue;
                    }

                    if (l_ユニティグ.Length < l_tip長閾値)
                    {
                        var l_先頭次数 = p_kmerインデックス.Get_入次数(l_塩基列.AsSpan(0, p_k長));
                        var l_末尾次数 = p_kmerインデックス.Get_出次数(l_塩基列.AsSpan(l_塩基列.Length - p_k長, p_k長));

                        // 片方の端が行き止まり(そちら向きに続きがない)であれば tip とみなし、
                        // どこにも合流しないため丸ごと除去してよい。
                        // 両端とも行き止まりの場合(=孤立した短い断片)も対象に含む。
                        if (l_先頭次数 == 0 || l_末尾次数 == 0)
                        {
                            V_除去_ユニティグ全体(p_kmerインデックス, l_塩基列, p_k長);
                            l_除去tip数++;
                            continue;
                        }
                    }

                    if (l_基準値 <= 0)
                    {
                        continue;
                    }

                    var l_剥がした数 = Get_剥がした数_低カバレッジ端(p_kmerインデックス, l_塩基列, p_k長, l_低カバレッジ閾値);
                    if (l_剥がした数 > 0)
                    {
                        l_剥がしたkmer数 += l_剥がした数;
                        l_トリミングしたunitig数++;
                    }
                }

                Console.WriteLine($"[GraphSimplifier] Iteration {l_反復}: examined {l_ユニティグ群.Count} unitig(s) " +
                    $"(tip threshold < {l_tip長閾値}bp, coverage baseline {l_基準値:0.#}), " +
                    $"removed {l_除去tip数} tip(s), trimmed {l_剥がしたkmer数} low-coverage k-mer(s) from {l_トリミングしたunitig数} unitig edge(s).");

                if (l_除去tip数 == 0 && l_剥がしたkmer数 == 0)
                {
                    return l_開始kmer;
                }

                // k-mer集合が縮小されたため、開始点を再検出してから次の反復へ。
                l_開始kmer = p_kmerインデックス.Get_開始kmer一覧();
            }

            return l_開始kmer;
        }

        private static List<string> Get_ユニティグ群(
            TrustedKmerIndex p_kmerインデックス, List<byte[]> p_開始kmer)
        {
            var l_walk結果 = UnitigMaker.Get_walk結果(p_kmerインデックス, p_開始kmer);

            List<string> l_ユニティグ群 = [];
            HashSet<string> l_既出 = [];
            foreach (var l_配列 in l_walk結果)
            {
                if (l_既出.Contains(l_配列) || l_既出.Contains(Util.V_逆相補(l_配列)))
                {
                    continue;
                }
                _ = l_既出.Add(l_配列);
                _ = l_既出.Add(Util.V_逆相補(l_配列));
                l_ユニティグ群.Add(l_配列);
            }
            return l_ユニティグ群;
        }

        /// <summary>
        /// unitig の両端から、カバレッジが閾値未満の k-mer が続く間だけ除去する。
        /// 先頭側と末尾側で除去範囲が重ならないよう互いの残り長で制限する。
        /// </summary>
        private static int Get_剥がした数_低カバレッジ端(
            TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長, double p_閾値)
        {
            var l_kmer数 = p_塩基列.Length - p_k長 + 1;
            if (l_kmer数 <= 0)
            {
                return 0;
            }

            var l_除去数 = 0;

            var l_先頭から = 0;
            while (l_先頭から < l_kmer数 && p_kmerインデックス.Get_カバレッジ(p_塩基列.AsSpan(l_先頭から, p_k長)) < p_閾値)
            {
                l_先頭から++;
            }

            var l_末尾から = 0;
            while (l_末尾から < l_kmer数 - l_先頭から && p_kmerインデックス.Get_カバレッジ(p_塩基列.AsSpan(l_kmer数 - 1 - l_末尾から, p_k長)) < p_閾値)
            {
                l_末尾から++;
            }

            for (var i = 0; i < l_先頭から; i++)
            {
                p_kmerインデックス.V_除去(p_塩基列.AsSpan(i, p_k長));
                l_除去数++;
            }
            for (var i = 0; i < l_末尾から; i++)
            {
                p_kmerインデックス.V_除去(p_塩基列.AsSpan(l_kmer数 - 1 - i, p_k長));
                l_除去数++;
            }

            return l_除去数;
        }

        /// <summary>
        /// unitigを構成する全k-merのカバレッジの単純平均。
        /// </summary>
        private static double Get_平均カバレッジ(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長)
        {
            ulong l_合計 = 0;
            var l_件数 = 0;
            for (var i = 0; i + p_k長 <= p_塩基列.Length; i++)
            {
                l_合計 += p_kmerインデックス.Get_カバレッジ(p_塩基列.AsSpan(i, p_k長));
                l_件数++;
            }
            return l_件数 == 0 ? 0 : (double)l_合計 / l_件数;
        }

        /// <summary>
        /// 全unitigの平均カバレッジの長さ加重中央値。多数を占めうる短い
        /// 断片(エラー由来のtip/バブル候補そのもの)に引きずられず、
        /// ゲノムの大部分を占める正しい主経路のカバレッジ水準を推定するため、
        /// 単純平均・単純中央値ではなく塩基数で重み付けした中央値を使う。
        /// </summary>
        private static double Get_長さ加重中央カバレッジ(
            TrustedKmerIndex p_kmerインデックス, List<string> p_ユニティグ群, int p_k長)
        {
            if (p_ユニティグ群.Count == 0)
            {
                return 0;
            }

            var l_組 = p_ユニティグ群
                .Select(x => (A_長さ: (long)x.Length, A_カバレッジ: Get_平均カバレッジ(p_kmerインデックス, Get_塩基列(x), p_k長)))
                .OrderBy(x => x.A_カバレッジ)
                .ToList();
            var l_総延長 = l_組.Sum(x => x.A_長さ);
            if (l_総延長 == 0)
            {
                return 0;
            }

            var l_半分 = l_総延長 / 2.0;
            long l_累積 = 0;
            foreach (var (l_長さ, l_カバレッジ) in l_組)
            {
                l_累積 += l_長さ;
                if (l_累積 >= l_半分)
                {
                    return l_カバレッジ;
                }
            }
            return l_組[^1].A_カバレッジ;
        }

        private static byte[] Get_塩基列(string p_配列)
        {
            var l_塩基列 = new byte[p_配列.Length];
            for (var i = 0; i < p_配列.Length; i++)
            {
                l_塩基列[i] = Util.Get_塩基ID(p_配列[i]);
            }
            return l_塩基列;
        }

        private static void V_除去_ユニティグ全体(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長)
        {
            for (var i = 0; i + p_k長 <= p_塩基列.Length; i++)
            {
                p_kmerインデックス.V_除去(p_塩基列.AsSpan(i, p_k長));
            }
        }
    }
}
