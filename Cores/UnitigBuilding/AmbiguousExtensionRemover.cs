using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// 控えに同水準の兄弟がいる一本道を、偶然そう見えているだけとみなして断つ
    /// </summary>
    internal static class AmbiguousExtensionRemover
    {
        #region 定数

        /// <summary>
        /// 控えの兄弟に対して、残った続きに要求する優勢の比
        /// </summary>
        private const double 優勢とみなす比 = 1.5D;

        #endregion

        #region テストメソッド

        /// <summary>
        /// 断った k-mer の数を返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <returns>除去した k-mer の件数</returns>
        public static int Get_除去数(TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            using var l_計測 = new StageTimer($"ambiguous-extension k={p_k長}");
            if (p_kmerインデックス.A_控えkmer数 == 0)
            {
                return 0;
            }

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            var l_候補 = p_kmerインデックス.Get_信頼kmer一覧()
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(l_スレッド数)
                .SelectMany(x => Get_断つべき後続(p_kmerインデックス, x))
                .ToList();

            var l_除去数 = 0;
            foreach (var l_後続 in l_候補)
            {
                if (p_kmerインデックス.Get_カバレッジ(l_後続) == 0UL)
                {
                    continue;
                }
                p_kmerインデックス.V_除去(l_後続);
                l_除去数++;
            }

            Logger.V_出力(メッセージID.曖昧な一本道の除去, l_除去数, l_候補.Count);
            return l_除去数;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k-mer の両方の向きについて、断つべき後続を返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        private static IEnumerable<byte[]> Get_断つべき後続(TrustedKmerIndex p_kmerインデックス, byte[] p_kmer)
        {
            if (Get_断つべき後続_片向き(p_kmerインデックス, p_kmer) is { } l_順鎖)
            {
                yield return l_順鎖;
            }
            if (Get_断つべき後続_片向き(p_kmerインデックス, Util.V_逆相補(p_kmer).ToArray()) is { } l_逆鎖)
            {
                yield return l_逆鎖;
            }
        }

        /// <summary>
        /// 一本道に見える後続のうち、控えに同水準の兄弟がいるものを返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_kmer"></param>
        /// <returns>断つべき後続、無ければ null</returns>
        private static byte[]? Get_断つべき後続_片向き(TrustedKmerIndex p_kmerインデックス, byte[] p_kmer)
        {
            var l_候補 = new byte[p_kmer.Length];
            p_kmer.AsSpan(1).CopyTo(l_候補);

            byte[]? l_後続 = null;
            var l_後続のカバレッジ = 0UL;
            var l_控えの最多 = 0UL;
            for (var l_塩基 = Consts.塩基ID.A; l_塩基 <= Consts.塩基ID.T; l_塩基++)
            {
                l_候補[^1] = l_塩基;
                var l_カバレッジ = p_kmerインデックス.Get_カバレッジ(l_候補);
                if (l_カバレッジ > 0UL)
                {
                    if (l_後続 is not null)
                    {
                        return null;
                    }
                    l_後続 = l_候補.ToArray();
                    l_後続のカバレッジ = l_カバレッジ;
                    continue;
                }
                l_控えの最多 = Math.Max(l_控えの最多, p_kmerインデックス.Get_控えカバレッジ(l_候補));
            }

            return l_後続 is null || l_控えの最多 == 0UL || l_後続のカバレッジ >= l_控えの最多 * 優勢とみなす比 ? null : p_kmerインデックス.Get_入次数(l_後続) == 1 ? l_後続 : null;
        }

        #endregion

    }
}
