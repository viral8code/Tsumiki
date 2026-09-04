namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer の出現回数ごとの分布(出現回数 -> 何種類のユニークk-merがその
    /// 回数を持つか)から、エラー由来の低頻度k-merと真のゲノム由来k-merを
    /// 分ける「谷」を推定する。実データでは典型的に出現回数1近辺に
    /// エラー由来の大きな山があり、その後カバレッジ相当の回数に真のゲノム
    /// 由来の山がある(二峰性分布)。谷の位置がカットオフの目安になる。
    /// </summary>
    internal static class KmerHistogram
    {
        /// <summary>
        /// 出現回数1から順に頻度を見ていき、単調減少から増加に転じる直前の回数を
        /// 谷(推奨カットオフ)として返す。単調減少のまま終わった場合や
        /// ヒストグラムが空の場合は null を返す。
        /// </summary>
        public static ulong? Get_推奨カットオフ(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限 = 1000)
        {
            if (p_ヒストグラム.Count == 0)
            {
                return null;
            }

            var l_最大出現回数 = Math.Min(p_ヒストグラム.Keys.Max(), p_走査上限);
            if (l_最大出現回数 < 2)
            {
                return null;
            }

            long? l_直前の頻度 = null;
            for (var l_出現回数 = 1UL; l_出現回数 <= l_最大出現回数; l_出現回数++)
            {
                var l_頻度 = p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L);
                if (l_直前の頻度 is { } l_前回 && l_頻度 > l_前回)
                {
                    return l_出現回数 - 1;
                }
                l_直前の頻度 = l_頻度;
            }
            return null;
        }

        /// <summary>
        /// 出現回数1から上限までのヒストグラムを1行にまとめた要約文字列を作る
        /// (ログ表示用)。
        /// </summary>
        public static string Get_要約(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_表示上限 = 20)
        {
            var l_項目 = new List<string>();
            var l_上限 = Math.Min(p_表示上限, p_ヒストグラム.Count == 0 ? 0 : p_ヒストグラム.Keys.Max());
            for (var l_出現回数 = 1UL; l_出現回数 <= l_上限; l_出現回数++)
            {
                l_項目.Add($"{l_出現回数}:{p_ヒストグラム.GetValueOrDefault(l_出現回数, 0L)}");
            }
            return string.Join(", ", l_項目);
        }
    }
}
