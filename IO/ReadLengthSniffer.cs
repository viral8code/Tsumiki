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
                while (l_長さ標本.Count < p_標本上限 && l_読み込み.Has続き())
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
        /// リード長ごとの本数を全リードから数える
        /// </summary>
        /// <param name="p_ファイルパス"></param>
        /// <remarks>
        /// 先頭だけを見ると、長さの違うライブラリを連結したファイルで後ろ側が見えない
        /// </remarks>
        /// <returns></returns>
        public static Dictionary<int, long> Get_リード長分布(string p_ファイルパス)
        {
            Dictionary<int, long> l_分布 = [];
            if (string.IsNullOrWhiteSpace(p_ファイルパス) || !File.Exists(p_ファイルパス))
            {
                return l_分布;
            }

            using var l_読み込み = new FastqReader(p_ファイルパス);
            while (l_読み込み.Has続き())
            {
                var l_長さ = l_読み込み.Get_次のリード_軽量().A_生リード!.Length;
                l_分布[l_長さ] = l_分布.GetValueOrDefault(l_長さ) + 1;
            }
            return l_分布;
        }

        /// <summary>
        /// k の梯子の上限を決めるリード長
        /// </summary>
        /// <param name="p_分布">リード長ごとの本数</param>
        /// <param name="p_必要な塩基の割合">この割合以上の塩基を供給する長さまでを採る</param>
        /// <remarks>
        /// 長いリードが少数でも混ざっていれば高い k へ届くが、そこでカバレッジが痩せては意味がない<br/>
        /// 「その長さ以上のリードが全体の塩基のうち何割を出しているか」で線を引く<br/>
        /// 長さが 1 種類しかない (単一ライブラリの) ファイルでは、その長さがそのまま返る
        /// </remarks>
        /// <returns>1 リードも無ければ null</returns>
        public static int? Get_梯子上限のリード長(IReadOnlyDictionary<int, long> p_分布, double p_必要な塩基の割合 = 0.1D)
        {
            if (p_分布.Count == 0)
            {
                return null;
            }

            var l_総塩基 = p_分布.Sum(x => (double)x.Key * x.Value);
            if (l_総塩基 <= 0D)
            {
                return null;
            }

            var l_長い順 = p_分布.Keys.OrderByDescending(x => x).ToList();
            var l_累積 = 0D;
            foreach (var l_長さ in l_長い順)
            {
                l_累積 += (double)l_長さ * p_分布[l_長さ];
                if (l_累積 / l_総塩基 >= p_必要な塩基の割合)
                {
                    return l_長さ;
                }
            }
            return l_長い順[^1];
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

        /// <summary>
        /// 全ライブラリを合わせた分布から、k の梯子の上限を決めるリード長を返す
        /// </summary>
        /// <param name="p_ライブラリ群">ライブラリごとのリードの組</param>
        /// <param name="p_分布">数え上げた分布の書き留め先</param>
        /// <returns></returns>
        public static int? Get_梯子上限のリード長(IEnumerable<(string A_リード1, string A_リード2)> p_ライブラリ群, out Dictionary<int, long> p_分布)
        {
            p_分布 = [];
            foreach (var (A_リード1, A_リード2) in p_ライブラリ群)
            {
                foreach (var l_パス in new[] { A_リード1, A_リード2 })
                {
                    foreach (var (l_長さ, l_本数) in Get_リード長分布(l_パス))
                    {
                        p_分布[l_長さ] = p_分布.GetValueOrDefault(l_長さ) + l_本数;
                    }
                }
            }
            return Get_梯子上限のリード長(p_分布);
        }

        #endregion
    }
}
