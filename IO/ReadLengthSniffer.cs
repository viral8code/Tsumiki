namespace Tsumiki.IO
{
    /// <summary>
    /// リードファイルの先頭を標本抽出して代表的なリード長を求める。
    ///
    /// k 長の自動選択に使う。リード長はライブラリの世代によって 75bp から
    /// 300bp まで大きく変わるうえ、適正な k はリード長にほぼ比例するため、
    /// 固定の既定値ではどのデータにも合わない。
    /// </summary>
    internal static class ReadLengthSniffer
    {
        /// <summary>
        /// 標本の中央値をリード長とする。トリミング済みのデータではリード長が
        /// ばらつくが、平均や最大値と違って中央値なら「大多数のリードが
        /// これ以上の長さを持つ」という保証に近い値になる。
        /// 1リードも読めなかった場合は null。
        /// </summary>
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
        /// ペアエンドの両ファイルから代表リード長を求め、短いほうを返す。
        /// k は「どちらのリードからも k-mer が取れる」必要があるため、
        /// 長いほうに合わせると短い側のリードが丸ごと使えなくなりうる。
        /// </summary>
        public static int? Get_代表リード長(string p_リード1のパス, string? p_リード2のパス, int p_標本上限 = 20_000)
        {
            var l_リード長1 = Get_代表リード長(p_リード1のパス, p_標本上限);
            if (string.IsNullOrWhiteSpace(p_リード2のパス))
            {
                return l_リード長1;
            }

            var l_リード長2 = Get_代表リード長(p_リード2のパス, p_標本上限);
            if (l_リード長1 is not { } l_長さ1)
            {
                return l_リード長2;
            }
            if (l_リード長2 is not { } l_長さ2)
            {
                return l_リード長1;
            }
            return Math.Min(l_長さ1, l_長さ2);
        }
    }
}
