namespace Tsumiki.Commons
{
    /// <summary>
    /// 中央値・分位点・N50・長さ加重中央値など、複数箇所で必要になる
    /// 分布の要約統計をまとめる
    /// </summary>
    /// <remarks>
    /// 「整列してから累積和が半分を超えた点を採る」
    /// という同じ骨格の実装がファイルごとに個別に書かれていたのを 1 箇所にする
    /// </remarks>
    internal static class StatsUtil
    {
        #region 公開メソッド

        /// <summary>
        /// 中央値を返す
        /// </summary>
        /// <param name="p_値一覧">元の値</param>
        /// <returns>中央値、値が無ければ 0</returns>
        public static int Get_中央値(IReadOnlyCollection<int> p_値一覧)
        {
            var l_整列済み = p_値一覧.Order().ToList();
            var l_中央 = l_整列済み.Count / 2;
            return l_整列済み.Count % 2 == 0 ? (l_整列済み[l_中央 - 1] + l_整列済み[l_中央]) / 2 : l_整列済み[l_中央];
        }

        /// <summary>
        /// 中央値を返す
        /// </summary>
        /// <param name="p_値一覧">元の値</param>
        /// <returns>中央値、値が無ければ 0</returns>
        public static double Get_中央値(IReadOnlyCollection<double> p_値一覧)
        {
            var l_整列済み = p_値一覧.Order().ToList();
            var l_中央 = l_整列済み.Count / 2;
            return l_整列済み.Count % 2 == 0 ? (l_整列済み[l_中央 - 1] + l_整列済み[l_中央]) / 2.0 : l_整列済み[l_中央];
        }

        /// <summary>
        /// 整列済みの一覧から分位点の値を取り出す
        /// </summary>
        /// <param name="p_整列済み"></param>
        /// <param name="p_分位"></param>
        public static int Get_分位点(IReadOnlyList<int> p_整列済み, double p_分位)
        {
            return p_整列済み[Math.Clamp((int)(p_分位 * (p_整列済み.Count - 1)), 0, p_整列済み.Count - 1)];
        }

        /// <summary>
        /// 長さで重み付けした値の中央値
        /// </summary>
        /// <param name="p_組"></param>
        /// <remarks>
        /// 累積長が総延長の半分を超えた点の
        /// 値を採る<br/>
        /// 短い断片が本数で多数を占めていても、実際の塩基の
        /// 大部分が属する水準を代表させたい場面 (単一コピー領域の
        /// カバレッジ基準値など) で使う
        /// </remarks>
        public static double Get_長さ加重中央値(IEnumerable<(long A_長さ, double A_値)> p_組)
        {
            var l_整列済み = p_組.OrderBy(x => x.A_値).ToList();
            var l_総延長 = l_整列済み.Sum(x => x.A_長さ);
            if (l_整列済み.Count == 0 || l_総延長 == 0)
            {
                return 0;
            }

            var l_半分 = l_総延長 / 2.0;
            var l_累積 = 0L;
            foreach (var (l_長さ, l_値) in l_整列済み)
            {
                l_累積 += l_長さ;
                if (l_累積 >= l_半分)
                {
                    return l_値;
                }
            }
            return l_整列済み[^1].A_値;
        }

        /// <summary>
        /// 長さ一覧から N50/L50 を求める
        /// </summary>
        /// <param name="p_長さ一覧"></param>
        /// <remarks>
        /// N50 は「この長さ以上の配列だけで
        /// 総延長の半分に達する」最小の長さ、L50 はそのために必要な本数
        /// </remarks>
        public static (long A_N50, int A_L50) Get_N50(IReadOnlyCollection<long> p_長さ一覧)
        {
            if (p_長さ一覧.Count == 0)
            {
                return (0, 0);
            }

            var l_降順 = p_長さ一覧.OrderByDescending(x => x).ToList();
            var l_半分 = l_降順.Sum() / 2.0;
            var l_累積 = 0L;
            for (var i = 0; i < l_降順.Count; i++)
            {
                l_累積 += l_降順[i];
                if (l_累積 >= l_半分)
                {
                    return (l_降順[i], i + 1);
                }
            }
            return (l_降順[^1], l_降順.Count);
        }

        #endregion
    }
}
