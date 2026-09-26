using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.UnitigBuilding;
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
        private const double C_tipとみなすカバレッジ比 = 0.5D;

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
        public static List<byte[]> V_除去_tip(TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長 = null, int? p_tip長閾値 = null, int p_最大反復数 = 30, double p_低カバレッジ比 = 0.2D, double p_tipカバレッジ比 = C_tipとみなすカバレッジ比, bool p_Is低カバレッジ端トリミング = true)
        {
            var l_基準長 = p_リード長 is { } l_リード長 ? Math.Min(p_k長, l_リード長 / 2) : p_k長;
            var l_tip長閾値 = p_tip長閾値 ?? Math.Max(10 * l_基準長, p_リード長 ?? 0);
            var l_開始kmer = p_kmerインデックス.Get_開始kmer一覧();

            for (var l_反復 = 1; l_反復 <= p_最大反復数; l_反復++)
            {
                var l_unitig群 = Get_Unitig情報(p_kmerインデックス, Get_Unitig群(p_kmerインデックス, l_開始kmer), p_k長);
                var l_基準値 = Get_長さ加重中央カバレッジ(l_unitig群);
                var l_信頼下限 = ConfigurationManager.A_スペクトルモデル?.A_信頼下限;

                var l_判定群 = l_unitig群
                    .AsParallel()
                    .AsOrdered()
                    .WithDegreeOfParallelism(Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数))
                    .Select(x => Get_判定(p_kmerインデックス, x.A_塩基列, x.A_平均カバレッジ, p_k長, l_tip長閾値, l_基準値, l_信頼下限, p_低カバレッジ比, p_tipカバレッジ比, p_Is低カバレッジ端トリミング))
                    .OfType<整理判定>()
                    .ToArray();
                var l_見送り = Get_見送る判定(p_kmerインデックス, l_判定群, p_k長);

                var l_除去tip数 = 0;
                var l_剥がしたkmer数 = 0;
                var l_トリミングしたunitig数 = 0;
                for (var i = 0; i < l_判定群.Length; i++)
                {
                    if (l_見送り.Contains(i))
                    {
                        continue;
                    }

                    var l_判定 = l_判定群[i];
                    if (l_判定.A_Is除去)
                    {
                        V_除去_Unitig全体(p_kmerインデックス, l_判定.A_塩基列, p_k長);
                        l_除去tip数++;
                        continue;
                    }

                    V_除去_端(p_kmerインデックス, l_判定, p_k長);
                    l_剥がしたkmer数 += l_判定.A_先頭から + l_判定.A_末尾から;
                    l_トリミングしたunitig数++;
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
        /// 反復の始めのグラフだけを見て、1 本の unitig を除くか、端を削るかを決める
        /// </summary>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <param name="p_塩基列">unitig の塩基 ID 列</param>
        /// <param name="p_平均カバレッジ">unitig の平均カバレッジ</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_tip長閾値">これより短い行き止まりを tip の候補にする</param>
        /// <param name="p_基準値">全体の基準カバレッジ</param>
        /// <param name="p_信頼下限">これ以上の平均カバレッジなら tip として除かない</param>
        /// <param name="p_低カバレッジ比">端を削る閾値の比</param>
        /// <param name="p_tipカバレッジ比">tip とみなすカバレッジ比</param>
        /// <param name="p_Is低カバレッジ端トリミング">端を削るか</param>
        /// <returns>判定、何もしないなら null</returns>
        private static 整理判定? Get_判定(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, double p_平均カバレッジ, int p_k長, int p_tip長閾値, double p_基準値, double? p_信頼下限, double p_低カバレッジ比, double p_tipカバレッジ比, bool p_Is低カバレッジ端トリミング)
        {
            if (p_塩基列.Length < p_k長)
            {
                return null;
            }

            var l_先頭次数 = p_kmerインデックス.Get_入次数(p_塩基列.AsSpan(0, p_k長));
            var l_末尾次数 = p_kmerインデックス.Get_出次数(p_塩基列.AsSpan(p_塩基列.Length - p_k長, p_k長));
            var l_kmer数 = p_塩基列.Length - p_k長 + 1;

            if (p_塩基列.Length < p_tip長閾値 && (l_先頭次数 == 0 || l_末尾次数 == 0))
            {
                var l_Is無条件信頼 = p_信頼下限 is { } l_下限 && p_平均カバレッジ >= l_下限;
                var l_比較基準 = Get_tip比較基準(p_kmerインデックス, p_塩基列, p_k長, l_先頭次数, l_末尾次数, p_基準値);
                if (!l_Is無条件信頼 && l_比較基準 is { } l_基準 && (l_基準 <= 0D || p_平均カバレッジ < l_基準 * p_tipカバレッジ比))
                {
                    return new 整理判定(p_塩基列, p_平均カバレッジ, true, l_kmer数, 0, Get_合流側kmer群(p_塩基列, p_k長, l_先頭次数 > 0, l_末尾次数 > 0));
                }
            }

            if (!p_Is低カバレッジ端トリミング)
            {
                return null;
            }

            var (l_先頭から, l_末尾から) = Get_低カバレッジ端の数(p_kmerインデックス, p_塩基列, p_k長, p_低カバレッジ比, l_先頭次数, l_末尾次数);
            return l_先頭から + l_末尾から == 0
                ? null
                : new 整理判定(p_塩基列, p_平均カバレッジ, false, l_先頭から, l_末尾から, Get_合流側kmer群(p_塩基列, p_k長, l_先頭から > 0 && l_先頭次数 > 0, l_末尾から > 0 && l_末尾次数 > 0));
        }

        /// <summary>
        /// 削る側の端のうち、合流点に接する端の k-mer を、合流点へ向かう向きで返す
        /// </summary>
        /// <param name="p_塩基列">unitig の塩基 ID 列</param>
        /// <param name="p_k長">k 長</param>
        /// <param name="p_Is先頭">先頭の端を含めるか</param>
        /// <param name="p_Is末尾">末尾の端を含めるか</param>
        /// <returns>合流点へ向かう向きの端 k-mer</returns>
        private static byte[][] Get_合流側kmer群(byte[] p_塩基列, int p_k長, bool p_Is先頭, bool p_Is末尾)
        {
            List<byte[]> l_結果 = [];
            if (p_Is先頭)
            {
                l_結果.Add(Get_逆相補(p_塩基列.AsSpan(0, p_k長)));
            }

            if (p_Is末尾)
            {
                l_結果.Add(p_塩基列.AsSpan(p_塩基列.Length - p_k長, p_k長).ToArray());
            }

            return [.. l_結果];
        }

        /// <summary>
        /// 判定をそのまま当てると合流点の枝がすべて消える所で、残す 1 本の判定を選ぶ
        /// </summary>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <param name="p_判定群">反復の始めのグラフで下した判定</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>当てずに見送る判定の番号</returns>
        private static HashSet<int> Get_見送る判定(TrustedKmerIndex p_kmerインデックス, 整理判定[] p_判定群, int p_k長)
        {
            HashSet<string> l_除くkmer = [];
            Dictionary<string, List<int>> l_合流点別 = [];
            for (var i = 0; i < p_判定群.Length; i++)
            {
                var l_判定 = p_判定群[i];
                foreach (var l_位置 in Get_削る位置(l_判定, p_k長))
                {
                    _ = l_除くkmer.Add(Get_正規キー(l_判定.A_塩基列.AsSpan(l_位置.A_開始, l_位置.A_長さ)));
                }

                foreach (var l_端 in l_判定.A_合流側kmer群)
                {
                    var l_合流点 = Convert.ToHexString(l_端.AsSpan(1));
                    if (!l_合流点別.TryGetValue(l_合流点, out var l_番号群))
                    {
                        l_番号群 = [];
                        l_合流点別[l_合流点] = l_番号群;
                    }

                    l_番号群.Add(i);
                }
            }

            HashSet<int> l_見送り = [];
            foreach (var (l_合流点, l_番号群) in l_合流点別)
            {
                var l_枝 = Convert.FromHexString(l_合流点);
                var l_兄弟 = new byte[l_枝.Length + 1];
                l_枝.CopyTo(l_兄弟, 1);
                var l_Is全枝消失 = true;
                for (var l_塩基 = Consts.塩基ID.A; l_塩基 <= Consts.塩基ID.T && l_Is全枝消失; l_塩基++)
                {
                    l_兄弟[0] = l_塩基;
                    l_Is全枝消失 = !p_kmerインデックス.Haskmer(l_兄弟) || l_除くkmer.Contains(Get_正規キー(l_兄弟));
                }

                if (l_Is全枝消失)
                {
                    _ = l_見送り.Add(l_番号群
                        .OrderByDescending(x => p_判定群[x].A_平均カバレッジ)
                        .ThenByDescending(x => p_判定群[x].A_塩基列.Length)
                        .ThenBy(x => Get_正規キー(p_判定群[x].A_塩基列), StringComparer.Ordinal)
                        .First());
                }
            }

            return l_見送り;
        }

        /// <summary>
        /// 判定で削る k-mer の範囲
        /// </summary>
        /// <param name="p_判定">判定</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>削る k-mer の開始位置と長さ</returns>
        private static IEnumerable<(int A_開始, int A_長さ)> Get_削る位置(整理判定 p_判定, int p_k長)
        {
            var l_kmer数 = p_判定.A_塩基列.Length - p_k長 + 1;
            for (var i = 0; i < p_判定.A_先頭から; i++)
            {
                yield return (i, p_k長);
            }

            for (var i = 0; i < p_判定.A_末尾から; i++)
            {
                yield return (l_kmer数 - 1 - i, p_k長);
            }
        }

        /// <summary>
        /// 判定どおりに unitig の両端の k-mer を集合から外す
        /// </summary>
        /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
        /// <param name="p_判定">判定</param>
        /// <param name="p_k長">k 長</param>
        private static void V_除去_端(TrustedKmerIndex p_kmerインデックス, 整理判定 p_判定, int p_k長)
        {
            foreach (var (l_開始, l_長さ) in Get_削る位置(p_判定, p_k長))
            {
                p_kmerインデックス.V_除去(p_判定.A_塩基列.AsSpan(l_開始, l_長さ));
            }
        }

        /// <summary>
        /// 向きによらない配列のキー
        /// </summary>
        /// <param name="p_配列">塩基 ID 列</param>
        /// <returns>配列と逆相補のうち小さい方の 16 進表記</returns>
        private static string Get_正規キー(ReadOnlySpan<byte> p_配列)
        {
            var l_逆相補 = Get_逆相補(p_配列);
            return Convert.ToHexString(p_配列.SequenceCompareTo(l_逆相補) <= 0 ? p_配列 : l_逆相補);
        }

        /// <summary>
        /// unitig の両端から、カバレッジが閾値未満の k-mer が続く数を数える
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_塩基列"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_低カバレッジ比"></param>
        /// <param name="p_先頭次数"></param>
        /// <param name="p_末尾次数"></param>
        /// <returns>先頭から削る数と末尾から削る数</returns>
        private static (int A_先頭から, int A_末尾から) Get_低カバレッジ端の数(TrustedKmerIndex p_kmerインデックス, byte[] p_塩基列, int p_k長, double p_低カバレッジ比, int p_先頭次数, int p_末尾次数)
        {
            var l_kmer数 = p_塩基列.Length - p_k長 + 1;
            if (l_kmer数 <= 0)
            {
                return (0, 0);
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

            return (l_先頭から, l_末尾から);
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
