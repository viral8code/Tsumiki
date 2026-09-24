using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// de Bruijn グラフの簡略化
    /// </summary>
    internal static class GraphSimplifier
    {
        #region 定数

        /// <summary>
        /// tip とみなすカバレッジ比
        /// </summary>
        private const double tipとみなすカバレッジ比 = 0.5D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// tip 除去と低カバレッジ端のトリミングを反復し、簡略化後の unitig 開始点を返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <param name="p_tip長閾値"></param>
        /// <param name="p_最大反復数"></param>
        /// <param name="p_低カバレッジ比"></param>
        /// <param name="p_tipカバレッジ比"></param>
        /// <param name="p_Is低カバレッジ端トリミング"></param>
        /// <returns></returns>
        public static List<byte[]> V_除去_tip(TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長 = null, int? p_tip長閾値 = null, int p_最大反復数 = 30, double p_低カバレッジ比 = 0.2D, double p_tipカバレッジ比 = tipとみなすカバレッジ比, bool p_Is低カバレッジ端トリミング = true)
        {
            var l_基準長 = p_リード長 is { } l_リード長 ? Math.Min(p_k長, l_リード長 / 2) : p_k長;
            var l_tip長閾値 = p_tip長閾値 ?? Math.Max(10 * l_基準長, p_リード長 ?? 0);
            var l_開始kmer = p_kmerインデックス.Get_開始kmer一覧();

            for (var l_反復 = 1; l_反復 <= p_最大反復数; l_反復++)
            {
                var l_unitig群 = Get_Unitig情報(p_kmerインデックス, Get_Unitig群(p_kmerインデックス, l_開始kmer), p_k長);
                var l_基準値 = Get_長さ加重中央カバレッジ(l_unitig群);

                var l_除去tip数 = 0;
                var l_剥がしたkmer数 = 0;
                var l_トリミングしたunitig数 = 0;
                foreach (var (l_塩基列, l_平均カバレッジ) in l_unitig群)
                {
                    if (l_塩基列.Length < p_k長)
                    {
                        continue;
                    }

                    var l_先頭次数 = p_kmerインデックス.Get_入次数(l_塩基列.AsSpan(0, p_k長));
                    var l_末尾次数 = p_kmerインデックス.Get_出次数(l_塩基列.AsSpan(l_塩基列.Length - p_k長, p_k長));

                    if (l_塩基列.Length < l_tip長閾値 && (l_先頭次数 == 0 || l_末尾次数 == 0))
                    {
                        var l_信頼下限 = ConfigurationManager.A_スペクトルモデル?.A_信頼下限;
                        var l_Is無条件信頼 = l_信頼下限 is { } l_下限 && l_平均カバレッジ >= l_下限;

                        var l_比較基準 = Get_tip比較基準(p_kmerインデックス, l_塩基列, p_k長, l_先頭次数, l_末尾次数, l_基準値);
                        if (!l_Is無条件信頼 && l_比較基準 is { } l_基準 && (l_基準 <= 0D || l_平均カバレッジ < l_基準 * p_tipカバレッジ比))
                        {
                            V_除去_Unitig全体(p_kmerインデックス, l_塩基列, p_k長);
                            l_除去tip数++;
                            continue;
                        }
                    }

                    var l_剥がした数 = p_Is低カバレッジ端トリミング ? Get_低カバレッジ端除去数(p_kmerインデックス, l_塩基列, p_k長, p_低カバレッジ比, l_先頭次数, l_末尾次数) : 0;
                    if (l_剥がした数 > 0)
                    {
                        l_剥がしたkmer数 += l_剥がした数;
                        l_トリミングしたunitig数++;
                    }
                }

                Logger.V_出力(メッセージID.グラフ単純化の反復, l_反復, l_unitig群.Length, l_tip長閾値, l_基準値, l_除去tip数, l_剥がしたkmer数, l_トリミングしたunitig数);

                if (l_除去tip数 == 0 && l_剥がしたkmer数 == 0)
                {
                    return l_開始kmer;
                }

                l_開始kmer = p_kmerインデックス.Get_開始kmer一覧();
            }

            return l_開始kmer;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 開始 k-mer から walk して、いまの unitig を列挙する
        /// </summary>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <param name="p_開始kmer">walk を始める k-mer</param>
        /// <returns>unitig の配列</returns>
        private static List<string> Get_Unitig群(TrustedKmerIndex p_kmerインデックス, List<byte[]> p_開始kmer)
        {
            var l_walk結果 = UnitigMaker.Get_walk結果(p_kmerインデックス, p_開始kmer);

            List<string> l_unitig群 = [];
            HashSet<string> l_既出 = [];
            foreach (var l_配列 in l_walk結果)
            {
                if (l_既出.Contains(l_配列) || l_既出.Contains(Util.V_逆相補(l_配列)))
                {
                    continue;
                }
                _ = l_既出.Add(l_配列);
                _ = l_既出.Add(Util.V_逆相補(l_配列));
                l_unitig群.Add(l_配列);
            }
            return l_unitig群;
        }

        /// <summary>
        /// tip の平均カバレッジと比べる相手
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_塩基列"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_先頭次数"></param>
        /// <param name="p_末尾次数"></param>
        /// <param name="p_基準値">全体の基準値</param>
        /// <returns>比べる相手のカバレッジ、除去の対象にしない場合は null</returns>
        private static double? Get_tip比較基準(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長, int p_先頭次数, int p_末尾次数, double p_基準値)
        {
            if (p_先頭次数 == 0 && p_末尾次数 == 0)
            {
                return p_基準値;
            }

            var l_接合側kmer = p_先頭次数 == 0 ? p_塩基列.AsSpan(p_塩基列.Length - p_k長, p_k長).ToArray() : Get_逆相補(p_塩基列.AsSpan(0, p_k長));
            var l_対抗 = Get_対抗カバレッジ(p_kmerインデックス, l_接合側kmer);
            return l_対抗 > 0UL ? l_対抗 : null;
        }

        /// <summary>
        /// 端の k-mer と同じ接合点へ合流する、別の枝の k-mer の最大カバレッジ
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_端kmer">接合点へ向かう向きにした端の k-mer</param>
        /// <returns>競合する枝が無ければ 0</returns>
        private static ulong Get_対抗カバレッジ(TrustedKmerIndex p_kmerインデックス, byte[] p_端kmer)
        {
            if (p_kmerインデックス.Get_出次数(p_端kmer) == 0)
            {
                return 0UL;
            }

            var l_兄弟 = (byte[])p_端kmer.Clone();
            var l_対抗 = 0UL;
            for (var l_塩基 = Consts.塩基ID.A; l_塩基 <= Consts.塩基ID.T; l_塩基++)
            {
                if (l_塩基 == p_端kmer[0])
                {
                    continue;
                }
                l_兄弟[0] = l_塩基;
                l_対抗 = Math.Max(l_対抗, p_kmerインデックス.Get_カバレッジ(l_兄弟));
            }
            return l_対抗;
        }

        /// <summary>
        /// unitig の両端から、カバレッジが閾値未満の k-mer が続く間だけ除去する
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_塩基列"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_低カバレッジ比"></param>
        /// <param name="p_先頭次数"></param>
        /// <param name="p_末尾次数"></param>
        /// <returns></returns>
        private static int Get_低カバレッジ端除去数(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長, double p_低カバレッジ比, int p_先頭次数, int p_末尾次数)
        {
            var l_kmer数 = p_塩基列.Length - p_k長 + 1;
            if (l_kmer数 <= 0)
            {
                return 0;
            }

            double? l_自身の中央値 = null;
            var l_先頭閾値 = p_低カバレッジ比 * (p_先頭次数 == 0 ? l_自身の中央値 ??= Get_中央カバレッジ(p_kmerインデックス, p_塩基列, p_k長) : Get_対抗カバレッジ(p_kmerインデックス, Get_逆相補(p_塩基列.AsSpan(0, p_k長))));
            var l_末尾閾値 = p_低カバレッジ比 * (p_末尾次数 == 0 ? l_自身の中央値 ?? Get_中央カバレッジ(p_kmerインデックス, p_塩基列, p_k長) : Get_対抗カバレッジ(p_kmerインデックス, p_塩基列.AsSpan(p_塩基列.Length - p_k長, p_k長).ToArray()));

            var l_先頭から = 0;
            while (l_先頭から < l_kmer数 && p_kmerインデックス.Get_カバレッジ(p_塩基列.AsSpan(l_先頭から, p_k長)) < l_先頭閾値)
            {
                l_先頭から++;
            }

            var l_末尾から = 0;
            while (l_末尾から < l_kmer数 - l_先頭から && p_kmerインデックス.Get_カバレッジ(p_塩基列.AsSpan(l_kmer数 - 1 - l_末尾から, p_k長)) < l_末尾閾値)
            {
                l_末尾から++;
            }

            for (var i = 0; i < l_先頭から; i++)
            {
                p_kmerインデックス.V_除去(p_塩基列.AsSpan(i, p_k長));
            }
            for (var i = 0; i < l_末尾から; i++)
            {
                p_kmerインデックス.V_除去(p_塩基列.AsSpan(l_kmer数 - 1 - i, p_k長));
            }

            return l_先頭から + l_末尾から;
        }

        /// <summary>
        /// unitig を構成する全 k-mer のカバレッジの中央値
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_塩基列"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        private static double Get_中央カバレッジ(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長)
        {
            var l_カバレッジ = new List<double>(p_塩基列.Length - p_k長 + 1);
            for (var i = 0; i + p_k長 <= p_塩基列.Length; i++)
            {
                l_カバレッジ.Add(p_kmerインデックス.Get_カバレッジ(p_塩基列.AsSpan(i, p_k長)));
            }
            return l_カバレッジ.Count == 0 ? 0D : StatsUtil.Get_中央値(l_カバレッジ);
        }

        /// <summary>
        /// unitig を構成する全 k-mer のカバレッジの単純平均
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_塩基列"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        private static double Get_平均カバレッジ(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長)
        {
            var l_合計 = 0UL;
            var l_件数 = 0;
            for (var i = 0; i + p_k長 <= p_塩基列.Length; i++)
            {
                l_合計 += p_kmerインデックス.Get_カバレッジ(p_塩基列.AsSpan(i, p_k長));
                l_件数++;
            }
            return l_件数 == 0 ? 0D : (double)l_合計 / l_件数;
        }

        /// <summary>
        /// 各 unitig の塩基列と平均カバレッジ
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_unitig群"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        private static (byte[] A_塩基列, double A_平均カバレッジ)[] Get_Unitig情報(TrustedKmerIndex p_kmerインデックス, List<string> p_unitig群, int p_k長)
        {
            return [.. p_unitig群
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
        /// 全 unitig の平均カバレッジの長さ加重中央値
        /// </summary>
        /// <param name="p_unitig群"></param>
        /// <returns></returns>
        private static double Get_長さ加重中央カバレッジ((byte[] A_塩基列, double A_平均カバレッジ)[] p_unitig群)
        {
            return StatsUtil.Get_長さ加重中央値(p_unitig群.Select(x => ((long)x.A_塩基列.Length, x.A_平均カバレッジ)));
        }

        /// <summary>
        /// unitig を構成する k-mer をすべて集合から外す
        /// </summary>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <param name="p_塩基列">unitig の塩基 ID 列</param>
        /// <param name="p_k長">k 長</param>
        private static void V_除去_Unitig全体(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長)
        {
            for (var i = 0; i + p_k長 <= p_塩基列.Length; i++)
            {
                p_kmerインデックス.V_除去(p_塩基列.AsSpan(i, p_k長));
            }
        }

        /// <summary>
        /// k-mer の逆相補
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        /// <returns></returns>
        private static byte[] Get_逆相補(ReadOnlySpan<byte> p_kmer)
        {
            var l_結果 = new byte[p_kmer.Length];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                l_結果[i] = (byte)(5 - p_kmer[p_kmer.Length - 1 - i]);
            }
            return l_結果;
        }

        #endregion
    }
}
