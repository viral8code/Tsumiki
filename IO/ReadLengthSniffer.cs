namespace Tsumiki.IO
{
    /// <summary>
    /// リードファイルの先頭を標本抽出して代表的なリード長を求める
    /// </summary>
    /// <remarks>
    /// k 長の自動選択に使う
    /// </remarks>
    internal static class ReadLengthSniffer
    {
        #region 公開メソッド

        /// <summary>
        /// 標本の中央値をリード長とする
        /// </summary>
        /// <param name="p_ファイルパス"></param>
        /// <param name="p_標本上限"></param>
        /// <remarks>
        /// トリミング済みデータでは長さがばらつくため、平均や最大値より中央値のほうが実態に近い<br/>
        /// 1 リードも読めなければ null
        /// </remarks>
        /// <returns></returns>
        public static int? Get_代表リード長(string p_ファイルパス, int p_標本上限 = 20_000)
        {
            var l_長さ標本 = new List<int>();
            using (var l_読み込み = new FastqReader(p_ファイルパス))
            {
                while (l_長さ標本.Count < p_標本上限 && l_読み込み.Get_続きがあるか())
                {
                    l_長さ標本.Add(l_読み込み.Get_次のリード_軽量().A_生リード!.Length);
                }
            }

            if (l_長さ標本.Count == 0)
            {
                return null;
            }

            l_長さ標本.Sort();
            return l_長さ標本[l_長さ標本.Count / 2];
        }

        /// <summary>
        /// 両ファイルの代表リード長のうち短いほう
        /// </summary>
        /// <param name="p_リード1のパス"></param>
        /// <param name="p_リード2のパス"></param>
        /// <param name="p_標本上限"></param>
        /// <remarks>
        /// 長いほうに合わせると短い側のリードが丸ごと使えなくなりうる
        /// </remarks>
        /// <returns></returns>
        public static int? Get_代表リード長(string p_リード1のパス, string? p_リード2のパス, int p_標本上限 = 20_000)
        {
            var l_リード長1 = Get_代表リード長(p_リード1のパス, p_標本上限);
            if (string.IsNullOrWhiteSpace(p_リード2のパス))
            {
                return l_リード長1;
            }

            var l_リード長2 = Get_代表リード長(p_リード2のパス, p_標本上限);
            return l_リード長1 is not { } l_長さ1 ? l_リード長2 : l_リード長2 is not { } l_長さ2 ? l_リード長1 : Math.Min(l_長さ1, l_長さ2);
        }

        #endregion
    }
}
