namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer の出現回数(count)ごとの分布(count -> 何種類のユニークk-merがその
    /// countを持つか)から、エラー由来の低頻度k-merと真のゲノム由来k-merを
    /// 分ける「谷」を推定する。実コインデータでは典型的に count=1近辺に
    /// エラー由来の大きな山があり、その後カバレッジ相当のcountに真のゲノム
    /// 由来の山がある(二峰性分布)。谷の位置がカットオフの目安になる。
    /// </summary>
    internal static class KmerHistogram
    {
        /// <summary>
        /// count=1から順に頻度を見ていき、単調減少から増加に転じる直前の count を
        /// 谷(推奨カットオフ)として返す。単調減少のまま終わった場合や
        /// ヒストグラムが空の場合は null を返す。
        /// </summary>
        public static ulong? SuggestCutoff(IReadOnlyDictionary<ulong, long> histogram, ulong maxCountToScan = 1000)
        {
            if (histogram.Count == 0)
            {
                return null;
            }

            var maxKey = Math.Min(histogram.Keys.Max(), maxCountToScan);
            if (maxKey < 2)
            {
                return null;
            }

            long? previous = null;
            for (var count = 1UL; count <= maxKey; count++)
            {
                var freq = histogram.GetValueOrDefault(count, 0L);
                if (previous is { } prevFreq && freq > prevFreq)
                {
                    return count - 1;
                }
                previous = freq;
            }
            return null;
        }

        /// <summary>
        /// count=1からmaxCountまでのヒストグラムを1行にまとめたサマリ文字列を作る
        /// (ログ表示用)。
        /// </summary>
        public static string FormatSummary(IReadOnlyDictionary<ulong, long> histogram, ulong maxCount = 20)
        {
            var parts = new List<string>();
            var top = Math.Min(maxCount, histogram.Count == 0 ? 0 : histogram.Keys.Max());
            for (var count = 1UL; count <= top; count++)
            {
                parts.Add($"{count}:{histogram.GetValueOrDefault(count, 0L)}");
            }
            return string.Join(", ", parts);
        }
    }
}
