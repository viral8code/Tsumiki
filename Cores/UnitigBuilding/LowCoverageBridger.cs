using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// カットオフ未満で控えた k-mer を辿って、信頼できる k-mer 集合の行き止まりを繋ぎ直す
    /// </summary>
    /// <remarks>
    /// GC に偏ったライブラリでは AT に富む領域のカバレッジが全体のカットオフを割り、実在する配列の途中でグラフが途切れる<br/>
    /// 続きが控えの k-mer だけで一本道に続き、再び信頼できる k-mer に合流する場合に限って足すことで、エラー由来の枝を持ち込まない
    /// </remarks>
    internal static class LowCoverageBridger
    {
        #region 定数

        /// <summary>
        /// 控えに残す出現回数の下限
        /// </summary>
        /// <remarks>
        /// 1 回しか見ていない k-mer はエラーと区別できない
        /// </remarks>
        public const ulong 控えの最小出現回数 = 2UL;

        /// <summary>
        /// 行き止まりの k-mer の出現回数に対して、続きに要求する出現回数の比
        /// </summary>
        /// <remarks>
        /// 隣り合う k-mer は k-1 塩基を共有するのでカバレッジは急には変わらず、行き止まりから桁違いに落ちる続きはエラー由来
        /// </remarks>
        private const double 行き止まりに対する最小比 = 0.1D;

        /// <summary>
        /// 続きの候補が複数あるとき、最多の候補に要求する次点との比
        /// </summary>
        private const double 優勢とみなす比 = 3D;

        /// <summary>
        /// 辿ってよい k-mer 数のリード長に対する倍率
        /// </summary>
        /// <remarks>
        /// 長く途切れた区間ほど、控えの k-mer だけで正しい一本道を選べている根拠が薄くなる
        /// </remarks>
        private const int 歩数上限のリード長倍率 = 4;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 行き止まりを架橋し、足した k-mer の数を返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 終わったら控えは手放す
        /// </remarks>
        /// <returns></returns>
        public static int Get_架橋kmer数(TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長)
        {
            using var l_計測 = new StageTimer($"low-coverage-bridge k={p_k長}");
            if (p_kmerインデックス.A_控えkmer数 == 0)
            {
                p_kmerインデックス.V_解放_控え();
                return 0;
            }

            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_歩数上限 = 歩数上限のリード長倍率 * Math.Max(p_リード長 ?? 0, p_k長);

            // 探索は読み取りだけなので並列に行い、集合への追加は全経路が出揃ってから決まった順で行う
            // 同じ k-mer を複数の経路が足すと先に足した側の出現回数が残るため、順序が変わると結果も変わる
            var l_行き止まり = p_kmerインデックス.Get_信頼kmer一覧()
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(l_スレッド数)
                .SelectMany(x => Get_行き止まりの向き(p_kmerインデックス, x))
                .ToList();

            var l_経路群 = l_行き止まり
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(l_スレッド数)
                .Select(x => Get_架橋経路(p_kmerインデックス, x, p_k長, l_歩数上限))
                .Where(x => x is not null)
                .ToList();

            var l_追加数 = 0;
            foreach (var l_経路 in l_経路群)
            {
                foreach (var (l_kmer, l_出現回数) in l_経路!)
                {
                    if (p_kmerインデックス.Try追加_信頼kmer(l_kmer, l_出現回数))
                    {
                        l_追加数++;
                    }
                }
            }

            Logger.V_出力(メッセージID.低カバレッジ架橋結果, l_経路群.Count, l_行き止まり.Count, l_追加数);
            p_kmerインデックス.V_解放_控え();
            return l_追加数;
        }

        /// <summary>
        /// 行き止まりから控えの k-mer を辿り、信頼できる k-mer に合流するまでの経路を返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_行き止まり">出次数 0 の向きにした k-mer</param>
        /// <param name="p_k長"></param>
        /// <param name="p_歩数上限"></param>
        /// <returns>足すべき k-mer と出現回数、合流できなければ null</returns>
        internal static List<(byte[] A_kmer, ulong A_出現回数)>? Get_架橋経路(TrustedKmerIndex p_kmerインデックス, byte[] p_行き止まり, int p_k長, int p_歩数上限)
        {
            var l_最小出現回数 = Math.Max(控えの最小出現回数, (ulong)Math.Ceiling(p_kmerインデックス.Get_カバレッジ(p_行き止まり) * 行き止まりに対する最小比));
            var l_現在 = (byte[])p_行き止まり.Clone();
            var l_候補 = new byte[p_k長];
            List<(byte[] A_kmer, ulong A_出現回数)> l_経路 = [];
            HashSet<UInt128> l_通過済み = [];

            for (var l_歩数 = 0; l_歩数 < p_歩数上限; l_歩数++)
            {
                l_現在.AsSpan(1).CopyTo(l_候補);
                var l_信頼数 = 0;
                var l_最多 = 0UL;
                var l_次点 = 0UL;
                var l_最多の塩基 = Consts.塩基ID.A;
                for (var l_塩基 = Consts.塩基ID.A; l_塩基 <= Consts.塩基ID.T; l_塩基++)
                {
                    l_候補[^1] = l_塩基;
                    if (p_kmerインデックス.Haskmer(l_候補))
                    {
                        l_信頼数++;
                        continue;
                    }

                    var l_出現回数 = p_kmerインデックス.Get_控えカバレッジ(l_候補);
                    if (l_出現回数 > l_最多)
                    {
                        l_次点 = l_最多;
                        l_最多 = l_出現回数;
                        l_最多の塩基 = l_塩基;
                    }
                    else if (l_出現回数 > l_次点)
                    {
                        l_次点 = l_出現回数;
                    }
                }

                // 行き止まりの直後に信頼できる k-mer は無い (出次数 0) ので、合流した経路は必ず控えを 1 つ以上含む
                if (l_信頼数 > 0)
                {
                    return l_信頼数 == 1 ? l_経路 : null;
                }

                if (l_最多 < l_最小出現回数)
                {
                    return null;
                }

                if (l_次点 > 0UL && l_最多 < l_次点 * 優勢とみなす比)
                {
                    return null;
                }

                l_候補[^1] = l_最多の塩基;
                if (!l_通過済み.Add(KmerPacking.TryGet_正規化キー(l_候補)))
                {
                    return null;
                }
                l_経路.Add((l_候補.ToArray(), l_最多));
                (l_現在, l_候補) = (l_候補, l_現在);
            }
            return null;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k-mer の両方の向きのうち、出次数 0 のものを返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        private static IEnumerable<byte[]> Get_行き止まりの向き(TrustedKmerIndex p_kmerインデックス, byte[] p_kmer)
        {
            if (p_kmerインデックス.Get_出次数(p_kmer) == 0)
            {
                yield return p_kmer;
            }
            var l_逆相補 = Util.V_逆相補(p_kmer).ToArray();
            if (p_kmerインデックス.Get_出次数(l_逆相補) == 0)
            {
                yield return l_逆相補;
            }
        }

        #endregion
    }
}
