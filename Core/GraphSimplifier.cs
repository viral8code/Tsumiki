using Tsumiki.Common;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// de Bruijnグラフ(TrustedKmerIndexが保持する厳密なk-mer集合)に対する
    /// グラフ簡略化。2種類のアーティファクト除去を行う:
    ///
    /// 1. tip 除去: 短い行き止まりの unitig(片端の次数が0)を丸ごと除去する。
    ///    行き止まりは定義上どこにも合流しないため、丸ごと除去しても他の
    ///    経路が使う配列を壊す心配がない。
    ///
    /// 2. 低カバレッジ端のトリミング: 真のバブル除去(2経路の合流点を
    ///    厳密に特定し、高頻度側だけを残す)ではなく簡易版。合流点を
    ///    明示的には特定せず、代わりに各unitigの両端から「カバレッジが
    ///    基準値(ゲノム全体の典型的な単一コピー相当の深度)に比べて
    ///    著しく低いk-merが続く間」だけを1つずつ剥がしていく。
    ///
    ///    重要な設計上の注意: 当初はunitig全体の平均カバレッジで判定
    ///    しようとしたが、SNP様の短い分岐(数塩基だけ異なりすぐ長い
    ///    共有配列に合流する)では、共有部分(高カバレッジ)に平均が
    ///    引きずられて低カバレッジ側を検出できないことがテストで判明した。
    ///    さらに、検出できてもunitig全体を丸ごと除去すると、合流後の
    ///    共有配列(高カバレッジ側の経路も使っている)まで消してしまい
    ///    別の経路を破壊してしまう。そのため「端から1kmerずつ、カバレッジが
    ///    基準値比で低い間だけ剥がす」方式にした。エラー由来の分岐は
    ///    合流点までの区間(=k-1個のk-mer)が本来低カバレッジのはずなので、
    ///    この区間だけを正確に剥がし、合流後の共有配列には手を付けない。
    /// </summary>
    internal static class GraphSimplifier
    {
        /// <summary>
        /// tip 除去と低カバレッジ端のトリミングを反復的に行い、
        /// 簡略化後の「unitig開始点」リストを返す。除去のたびにunitigを
        /// ゼロから再構築して次数・カバレッジを再評価するため、最大
        /// p_最大反復数 回まで反復する(それ以上変化がなくなった時点で早期終了)。
        ///
        /// 既定の長さ閾値(k*10)は、合成データ(300kbゲノム・1%エラー率・
        /// エラー訂正後)での実測に基づく: k*2(Velvet等でよく使われる値)
        /// では訂正しきれず残った少数のエラー由来の分岐がまだ長めの
        /// tip として残ってしまい、k*10まで広げてようやく大部分を吸収できた
        /// (収束までに約13反復かかったため既定の最大反復数も余裕を見て
        /// 30とした)。この長さ閾値は tip 除去(行き止まりの丸ごと除去)
        /// にのみ適用し、低カバレッジ端のトリミングは(合流先の長さに
        /// 依存せず正しく判定できるため)unitigの長さによらず全件に適用する。
        ///
        /// 既定のカバレッジ閾値(基準値の20%)は、SPAdes等が「誤った接続」
        /// 除去に使う値と同程度の、一般的に保守的とされる水準。基準値は
        /// 全unitigの平均カバレッジの長さ加重中央値(短い断片が多数を占めても、
        /// ゲノムの大部分を占める正しい主経路のカバレッジに引きずられにくくするため)。
        /// </summary>
        public static List<byte[]> V_除去_tip(
            TrustedKmerIndex p_kmerインデックス,
            int p_k長,
            int? p_tip長閾値 = null,
            int p_最大反復数 = 30,
            double p_低カバレッジ比 = 0.2)
        {
            var l_tip長閾値 = p_tip長閾値 ?? (10 * p_k長);
            var l_ユニティグ構築 = new UnitigMaker(p_kmerインデックス);
            var l_開始kmer = p_kmerインデックス.Get_開始kmer一覧();

            for (var l_反復 = 1; l_反復 <= p_最大反復数; l_反復++)
            {
                var l_ユニティグ群 = Get_ユニティグ群(l_ユニティグ構築, l_開始kmer);
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

        private static List<string> Get_ユニティグ群(UnitigMaker p_ユニティグ構築, List<byte[]> p_開始kmer)
        {
            List<string> l_ユニティグ群 = [];
            HashSet<string> l_既出 = [];
            foreach (var l_kmer in p_開始kmer)
            {
                var l_ユニティグ = p_ユニティグ構築.Get_ユニティグ(l_kmer);
                if (l_既出.Contains(l_ユニティグ.A_配列) || l_既出.Contains(Util.V_逆相補(l_ユニティグ.A_配列)))
                {
                    continue;
                }
                _ = l_既出.Add(l_ユニティグ.A_配列);
                _ = l_既出.Add(Util.V_逆相補(l_ユニティグ.A_配列));
                l_ユニティグ群.Add(l_ユニティグ.A_配列);
            }
            return l_ユニティグ群;
        }

        /// <summary>
        /// unitigの両端から、カバレッジが閾値未満のk-merが続く間だけ
        /// 1つずつ信頼できる集合から除去する(合流点に達して初めて
        /// 閾値以上のk-merが現れたら、そこで止めて先には進まない)。
        /// 先頭側と末尾側で除去範囲が重ならないよう、互いの残り長で制限する。
        /// 戻り値: 除去したk-merの総数。
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
