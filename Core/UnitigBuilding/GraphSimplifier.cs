using Tsumiki.Common;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Core.UnitigBuilding
{
    /// <summary>
    /// de Bruijn グラフの簡略化<br/>
    /// 2 種類のアーティファクトを除去する<br/>
    /// 1. tip 除去: 短い行き止まりの unitig を丸ごと除去する<br/>
    /// 行き止まりは
    /// 定義上どこにも合流しないので、他の経路が使う配列を壊さない<br/>
    /// 2. 低カバレッジ端のトリミング: 合流点を特定せず、各 unitig の両端から
    /// カバレッジが基準値比で著しく低い k-mer が続く間だけ剥がす<br/>
    /// unitig 全体の平均で判定してはいけない<br/>
    /// SNP 様の短い分岐では共有部分の
    /// 高カバレッジに平均が引きずられて検出できず、仮に検出できても unitig 全体を
    /// 除去すると合流後の共有配列まで消して別の経路を壊す<br/>
    /// エラー由来の分岐は合流点までの区間だけが低カバレッジなので、そこだけ剥がす
    /// </summary>
    internal static class GraphSimplifier
    {
        /// <summary>
        /// tip 除去と低カバレッジ端のトリミングを反復し、簡略化後の unitig 開始点を返す<br/>
        /// 除去のたびに unitig を再構築して次数とカバレッジを評価し直すため反復する<br/>
        /// 長さ閾値は tip 除去にのみ適用する<br/>
        /// 低カバレッジ端のトリミングは
        /// 合流先の長さに依存せず判定できるため、長さによらず全件に適用する<br/>
        /// カバレッジの基準値は長さ加重中央値を使う<br/>
        /// 単純平均や単純中央値だと
        /// 本数の多い短い断片に引きずられ、主経路の水準から外れる
        /// </summary>
        public static List<byte[]> V_除去_tip(
            TrustedKmerIndex p_kmerインデックス,
            int p_k長,
            int? p_リード長 = null,
            int? p_tip長閾値 = null,
            int p_最大反復数 = 30,
            double p_低カバレッジ比 = 0.2,
            double p_tipカバレッジ比 = Consts.tipとみなすカバレッジ比)
        {
            // k がリード長の半分を超えると、k を基準にした閾値は実配列まで
            // 巻き込むほど長くなるため min(k, リード長/2) を基準に取る
            var l_基準長 = p_リード長 is { } l_リード長 ? Math.Min(p_k長, l_リード長 / 2) : p_k長;
            var l_tip長閾値 = p_tip長閾値 ?? Math.Max(10 * l_基準長, p_リード長 ?? 0);
            var l_開始kmer = p_kmerインデックス.Get_開始kmer一覧();

            for (var l_反復 = 1; l_反復 <= p_最大反復数; l_反復++)
            {
                var l_ユニティグ群 = Get_ユニティグ情報(p_kmerインデックス, Get_ユニティグ群(p_kmerインデックス, l_開始kmer), p_k長);
                var l_基準値 = Get_長さ加重中央カバレッジ(l_ユニティグ群);
                var l_低カバレッジ閾値 = l_基準値 * p_低カバレッジ比;

                var l_除去tip数 = 0;
                var l_剥がしたkmer数 = 0;
                var l_トリミングしたunitig数 = 0;
                foreach (var (l_塩基列, l_平均カバレッジ) in l_ユニティグ群)
                {
                    if (l_塩基列.Length < p_k長)
                    {
                        continue;
                    }

                    if (l_塩基列.Length < l_tip長閾値)
                    {
                        var l_先頭次数 = p_kmerインデックス.Get_入次数(l_塩基列.AsSpan(0, p_k長));
                        var l_末尾次数 = p_kmerインデックス.Get_出次数(l_塩基列.AsSpan(l_塩基列.Length - p_k長, p_k長));

                        // 片方の端が行き止まり (そちら向きに続きがない) であれば tip の候補
                        // 両端とも行き止まりの場合 (=孤立した短い断片) も対象に含む
                        //
                        // ただし行き止まりであること自体は誤りの証拠にならない
                        // 実ゲノムでも
                        // カバレッジが切れた箇所では両端が行き止まりになる
                        // エラー由来の枝は
                        // カバレッジがカットオフ付近に留まるので、基準値と比べて明らかに
                        // 低いものだけを除去する
                        // これを見ないと実配列まで消してゲノム被覆率を落とす
                        if (l_先頭次数 == 0 || l_末尾次数 == 0)
                        {
                            // k-mer スペクトルの混合モデルが「無条件に信頼してよい」と
                            // 判定したカバレッジ以上なら、行き止まりに見えても除去しない
                            // 比較対象の基準値は反復のたびに evolving graph から
                            // 再計算するため、低カバレッジのノイズに引きずられて
                            // ぶれることがある
                            // 信頼下限はそれとは独立に、カットオフと
                            // 同じモデルから求めた絶対値なので、そのぶれを踏まない
                            // 安全弁になる
                            var l_信頼下限 = ConfigurationManager.A_スペクトルモデル?.A_信頼下限;
                            var l_無条件に信頼できるか = l_信頼下限 is { } l_下限 && l_平均カバレッジ >= l_下限;

                            if (!l_無条件に信頼できるか
                                && (l_基準値 <= 0 || l_平均カバレッジ < l_基準値 * p_tipカバレッジ比))
                            {
                                V_除去_ユニティグ全体(p_kmerインデックス, l_塩基列, p_k長);
                                l_除去tip数++;
                                continue;
                            }
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

                Logger.V_出力(
                    メッセージID.グラフ単純化の反復, l_反復, l_ユニティグ群.Length, l_tip長閾値, l_基準値, l_除去tip数, l_剥がしたkmer数, l_トリミングしたunitig数);

                if (l_除去tip数 == 0 && l_剥がしたkmer数 == 0)
                {
                    return l_開始kmer;
                }

                // k-mer 集合が縮小されたため、開始点を再検出してから次の反復へ
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
        /// unitig の両端から、カバレッジが閾値未満の k-mer が続く間だけ除去する<br/>
        /// 先頭側と末尾側で除去範囲が重ならないよう互いの残り長で制限する
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
        /// unitig を構成する全 k-mer のカバレッジの単純平均
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
        /// 各 unitig の塩基列と平均カバレッジ<br/>
        /// 基準値の算出と tip 判定の
        /// 両方が同じ値を使うため、まとめて 1 回だけ求める<br/>
        /// 全 unitig の全 k-mer を引くので反復のたびに数百万回のハッシュ引きになる<br/>
        /// 読み取りのみなので並列に行う
        /// </summary>
        private static (byte[] A_塩基列, double A_平均カバレッジ)[] Get_ユニティグ情報(
            TrustedKmerIndex p_kmerインデックス, List<string> p_ユニティグ群, int p_k長)
        {
            return [.. p_ユニティグ群
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数))
                .Select(x =>
                {
                    var l_塩基列 = Util.V_変換_塩基列(x);
                    return (l_塩基列, Get_平均カバレッジ(p_kmerインデックス, l_塩基列, p_k長));
                })];
        }

        /// <summary>
        /// 全 unitig の平均カバレッジの長さ加重中央値<br/>
        /// 多数を占めうる短い
        /// 断片 (エラー由来の tip/バブル候補そのもの) に引きずられず、
        /// ゲノムの大部分を占める正しい主経路のカバレッジ水準を推定するため、
        /// 単純平均・単純中央値ではなく塩基数で重み付けした中央値を使う
        /// </summary>
        private static double Get_長さ加重中央カバレッジ(
            (byte[] A_塩基列, double A_平均カバレッジ)[] p_ユニティグ群)
        {
            return StatsUtil.Get_長さ加重中央値(
                p_ユニティグ群.Select(x => ((long)x.A_塩基列.Length, x.A_平均カバレッジ)));
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
