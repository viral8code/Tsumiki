using System.Text;
using Tsumiki.Commons;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 配列の端どうしを、跨ぐリード 2 本以上で繋ぐ<br/>
    /// 端の相手は、リードの続きに現れる一意な配列 (先頭の一意な 25-mer) で見つけ、互いに相手がただ 1 つに決まる組だけを Scaffolder.Get_リードで埋めた繋ぎ目 で埋める<br/>
    /// 配列 i の順向きを頂点 2i、逆相補を頂点 2i+1 とする
    /// </summary>
    internal static class ReadEndJoiner
    {
        #region 定数

        /// <summary>
        /// 端の錨の長さ (ReadMinimizerIndex の最短の問い合わせ長)
        /// </summary>
        private const int C_錨長 = ReadMinimizerIndex.C_最短の問い合わせ長;

        /// <summary>
        /// 一意かどうかを見る k-mer の長さ (Scaffolder.C_右の錨長)
        /// </summary>
        private const int C_一意の長さ = Scaffolder.C_右の錨長;

        /// <summary>
        /// 相手の配列の先頭のうち、一意な 25-mer を索引に入れる範囲の長さ
        /// </summary>
        private const int C_相手を探す先頭の長さ = 200;

        /// <summary>
        /// 相手を決める票数の下限、および埋めるのに要るリード数
        /// </summary>
        private const int C_要るリード数 = 2;

        /// <summary>
        /// 首位の票が次点の何倍以上か
        /// </summary>
        private const int C_首位の優勢比 = 3;

        /// <summary>
        /// Get_前後群 に渡す、続きを見るリード数の上限
        /// </summary>
        private const int C_続きを見るリード数の上限 = 64;

        /// <summary>
        /// これより短い配列は繋がない (両端を削っても配列が残るように)
        /// </summary>
        private const int C_対象の最短長 = 150;

        /// <summary>
        /// 繋ぎ目で左の末尾・右の先頭から削ってみる長さの上限<br/>
        /// 高い k で k-mer が欠けて、k−1 未満で重なる断片どうしも、重なりごと見つけられるようにする
        /// </summary>
        private const int C_削る上限 = 100;

        /// <summary>
        /// 端が反復かを確かめる、末尾の錨の位置 (末尾から錨の終わりまでの距離)<br/>
        /// ゲノムで 2 コピー以上ある配列の端どうしを繋ぐと、間のコピーを飛ばすため
        /// </summary>
        private static readonly int[] C_反復を見る位置群 = [0, 50, 100, 150];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 互いに一意な端の組を、跨ぐリードの続きで埋め、鎖にして配列を繋ぐ
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安 (Scaffolder.Get_一意の出現数) 、分からなければ 0</param>
        /// <param name="p_候補数">互いに一意な端の組 (埋める前) の数</param>
        /// <param name="p_繋いだ数">連結に使った辺の数</param>
        /// <returns>
        /// 繋いだ後の配列の ID と配列<br/>
        /// 鎖に入らなかった配列は元の順と ID のまま
        /// </returns>
        public static List<(string A_ID, string A_配列)> Get_繋いだ配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, ReadMinimizerIndex p_索引, int p_リード長, int p_一意の出現数, out int p_候補数, out int p_繋いだ数)
        {
            var l_配列数 = p_配列群.Count;
            var l_頂点 = new string[2 * l_配列数];
            var l_対象 = new bool[l_配列数];
            for (var i = 0; i < l_配列数; i++)
            {
                l_頂点[2 * i] = p_配列群[i].A_配列;
                l_頂点[(2 * i) + 1] = Util.V_逆相補_曖昧塩基あり(p_配列群[i].A_配列);
                l_対象[i] = p_配列群[i].A_配列.Length >= C_対象の最短長 && !p_配列群[i].A_ID.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase);
            }

            var l_登録 = Get_一意な先頭の登録(l_頂点, l_対象);
            var l_相手 = Get_相手(l_頂点, l_対象, p_索引, p_リード長, l_登録);
            var l_出る = new (int A_先, int A_左から削る長さ, string A_埋める配列, int A_右から削る長さ)?[l_頂点.Length];
            var l_入る = new int[l_頂点.Length];
            Array.Fill(l_入る, -1);
            var l_候補数 = 0;
            for (var v = 0; v < l_頂点.Length; v++)
            {
                var w = l_相手[v];
                if (w < 0 || l_相手[w ^ 1] != (v ^ 1) || v >= (w ^ 1))
                {
                    continue;
                }

                l_候補数++;
                if (Is端が反復(p_索引, l_頂点[v], p_一意の出現数) || Is端が反復(p_索引, l_頂点[w ^ 1], p_一意の出現数))
                {
                    continue;
                }

                if (Scaffolder.Get_リードで埋めた繋ぎ目(p_索引, new StringBuilder(l_頂点[v]), l_頂点[w], p_リード長, p_一意の出現数, C_要るリード数, C_削る上限) is not { } l_辺)
                {
                    continue;
                }

                l_出る[v] = (w, l_辺.A_左から削る長さ, l_辺.A_埋める配列, l_辺.A_右から削る長さ);
                l_出る[w ^ 1] = (v ^ 1, l_辺.A_右から削る長さ, Util.V_逆相補_曖昧塩基あり(l_辺.A_埋める配列), l_辺.A_左から削る長さ);
                l_入る[w] = v;
                l_入る[v ^ 1] = w ^ 1;
            }

            p_候補数 = l_候補数;
            return Get_鎖にして書き出す(p_配列群, l_頂点, l_出る, l_入る, out p_繋いだ数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 向き付き配列の右端の付近 (C_反復を見る位置群) の錨のどれかが反復の中 (<see cref = "Scaffolder.Is反復の中"/>) か
        /// </summary>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_配列">向き付き配列</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安、分からなければ 0</param>
        /// <returns>反復の中なら true</returns>
        private static bool Is端が反復(ReadMinimizerIndex p_索引, string p_配列, int p_一意の出現数)
        {
            foreach (var l_距離 in C_反復を見る位置群)
            {
                var l_位置 = p_配列.Length - C_錨長 - l_距離;
                if (l_位置 < 0 || p_配列.AsSpan(l_位置, C_錨長).IndexOfAnyExcept(Consts.塩基文字) >= 0)
                {
                    continue;
                }

                if (Scaffolder.Is反復の中(p_索引, p_配列, l_位置, p_一意の出現数))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// 配列の 25-mer の A/C/G/T の塩基値 (A=0, C=1, G=2, T=3)、それ以外は -1
        /// </summary>
        /// <param name="p_塩基">塩基</param>
        /// <returns>塩基値、または -1</returns>
        private static int Get_塩基値(char p_塩基)
        {
            return p_塩基 switch
            {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => -1,
            };
        }

        /// <summary>
        /// 25-mer の正準形 (その配列と逆相補の小さいほう) を返す<br/>
        /// A/C/G/T 以外を含めば false
        /// </summary>
        /// <param name="p_断片">25-mer</param>
        /// <param name="p_正準">正準形の 2 ビット表現</param>
        /// <returns>正準形を求めたか</returns>
        private static bool Is成功_正準形(ReadOnlySpan<char> p_断片, out ulong p_正準)
        {
            var l_順 = 0UL;
            var l_逆 = 0UL;
            foreach (var l_塩基 in p_断片)
            {
                var l_値 = Get_塩基値(l_塩基);
                if (l_値 < 0)
                {
                    p_正準 = 0UL;
                    return false;
                }

                l_順 = (l_順 << 2) | (uint)l_値;
                l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_値) << (2 * (C_一意の長さ - 1)));
            }

            p_正準 = Math.Min(l_順, l_逆);
            return true;
        }

        /// <summary>
        /// 配列の 25-mer の正準形を、数えたい正準形の表 (数は 0 で入っている) に数え足す<br/>
        /// ローリングで 1 回だけ走る
        /// </summary>
        /// <param name="p_配列">順向きの配列 (対象外も含む)</param>
        /// <param name="p_数">正準形ごとの出現数 (表に無いものは数えない)</param>
        private static void V_数える(string p_配列, Dictionary<ulong, int> p_数)
        {
            const ulong l_マスク = (1UL << (2 * C_一意の長さ)) - 1UL;
            var l_順 = 0UL;
            var l_逆 = 0UL;
            var l_連続 = 0;
            foreach (var l_塩基 in p_配列)
            {
                var l_値 = Get_塩基値(l_塩基);
                if (l_値 < 0)
                {
                    l_連続 = 0;
                    l_順 = 0UL;
                    l_逆 = 0UL;
                    continue;
                }

                l_順 = ((l_順 << 2) | (uint)l_値) & l_マスク;
                l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_値) << (2 * (C_一意の長さ - 1)));
                l_連続 = Math.Min(l_連続 + 1, C_一意の長さ);
                if (l_連続 < C_一意の長さ)
                {
                    continue;
                }

                var l_正準 = Math.Min(l_順, l_逆);
                if (p_数.TryGetValue(l_正準, out var l_数))
                {
                    p_数[l_正準] = l_数 + 1;
                }
            }
        }

        /// <summary>
        /// 対象の頂点の先頭の一意な 25-mer を、向きのままの文字列をキーにして登録する<br/>
        /// 同じキーが別の頂点から来たものは使わない (登録から外し、以後も登録しない)
        /// </summary>
        /// <param name="p_頂点">向き付きの配列 (2i が順向き、2i+1 が逆相補)</param>
        /// <param name="p_対象">配列ごとの対象かどうか</param>
        /// <returns>キーから頂点と位置への表</returns>
        private static Dictionary<string, (int A_頂点, int A_位置)> Get_一意な先頭の登録(string[] p_頂点, bool[] p_対象)
        {
            var l_数 = new Dictionary<ulong, int>();
            for (var v = 0; v < p_頂点.Length; v++)
            {
                if (!p_対象[v / 2])
                {
                    continue;
                }

                var l_最後 = Math.Min(C_相手を探す先頭の長さ, p_頂点[v].Length - C_一意の長さ);
                for (var p = 0; p <= l_最後; p++)
                {
                    if (Is成功_正準形(p_頂点[v].AsSpan(p, C_一意の長さ), out var l_正準))
                    {
                        _ = l_数.TryAdd(l_正準, 0);
                    }
                }
            }

            for (var i = 0; i < p_頂点.Length; i += 2)
            {
                V_数える(p_頂点[i], l_数);
            }

            var l_登録 = new Dictionary<string, (int A_頂点, int A_位置)>();
            var l_禁止 = new HashSet<string>();
            for (var v = 0; v < p_頂点.Length; v++)
            {
                if (!p_対象[v / 2])
                {
                    continue;
                }

                var l_最後 = Math.Min(C_相手を探す先頭の長さ, p_頂点[v].Length - C_一意の長さ);
                for (var p = 0; p <= l_最後; p++)
                {
                    if (!Is成功_正準形(p_頂点[v].AsSpan(p, C_一意の長さ), out var l_正準) || l_数[l_正準] != 1)
                    {
                        continue;
                    }

                    var l_キー = p_頂点[v].Substring(p, C_一意の長さ);
                    if (l_禁止.Contains(l_キー))
                    {
                        continue;
                    }

                    if (!l_登録.TryGetValue(l_キー, out var l_既存))
                    {
                        l_登録[l_キー] = (v, p);
                    }
                    else if (l_既存.A_頂点 != v)
                    {
                        _ = l_登録.Remove(l_キー);
                        _ = l_禁止.Add(l_キー);
                    }
                }
            }

            return l_登録;
        }

        /// <summary>
        /// 続きを、左から 25-mer ずつ見て、最初に登録と当たった頂点を返す
        /// </summary>
        /// <param name="p_続き">リードの続き</param>
        /// <param name="p_登録">一意な先頭の登録</param>
        /// <returns>当たった頂点、当たらなければ null</returns>
        private static int? Get_最初に当たった頂点(string p_続き, Dictionary<string, (int A_頂点, int A_位置)> p_登録)
        {
            for (var p = 0; p + C_一意の長さ <= p_続き.Length; p++)
            {
                if (p_登録.TryGetValue(p_続き.Substring(p, C_一意の長さ), out var l_当たり))
                {
                    return l_当たり.A_頂点;
                }
            }

            return null;
        }

        /// <summary>
        /// 対象の各頂点について、末尾の錨を含むリードの続きから、相手の頂点を決める<br/>
        /// 票が首位 (C_要るリード数 以上、次点の C_首位の優勢比 倍以上) のときだけ相手とする
        /// </summary>
        /// <param name="p_頂点">向き付きの配列</param>
        /// <param name="p_対象">配列ごとの対象かどうか</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_登録">一意な先頭の登録</param>
        /// <returns>頂点ごとの相手の頂点、無ければ -1</returns>
        private static int[] Get_相手(string[] p_頂点, bool[] p_対象, ReadMinimizerIndex p_索引, int p_リード長, Dictionary<string, (int A_頂点, int A_位置)> p_登録)
        {
            var l_相手 = new int[p_頂点.Length];
            Array.Fill(l_相手, -1);
            for (var v = 0; v < p_頂点.Length; v++)
            {
                if (!p_対象[v / 2])
                {
                    continue;
                }

                var l_錨 = p_頂点[v].AsSpan(p_頂点[v].Length - C_錨長);
                if (l_錨.IndexOfAnyExcept(Consts.塩基文字) >= 0)
                {
                    continue;
                }

                var l_票 = new Dictionary<int, int>();
                foreach (var (_, l_続き) in p_索引.Get_前後群(l_錨, p_リード長, C_続きを見るリード数の上限))
                {
                    if (Get_最初に当たった頂点(l_続き, p_登録) is { } l_w && l_w != v && l_w != (v ^ 1))
                    {
                        l_票[l_w] = l_票.GetValueOrDefault(l_w) + 1;
                    }
                }

                var l_首位 = -1;
                var l_首位の票 = 0;
                var l_次点の票 = 0;
                foreach (var (l_w, l_数) in l_票)
                {
                    if (l_数 > l_首位の票)
                    {
                        l_次点の票 = l_首位の票;
                        l_首位の票 = l_数;
                        l_首位 = l_w;
                    }
                    else if (l_数 == l_首位の票)
                    {
                        l_次点の票 = l_数;
                    }
                    else if (l_数 > l_次点の票)
                    {
                        l_次点の票 = l_数;
                    }
                }

                if (l_首位 >= 0 && l_首位の票 >= C_要るリード数 && l_首位の票 >= C_首位の優勢比 * l_次点の票)
                {
                    l_相手[v] = l_首位;
                }
            }

            return l_相手;
        }

        /// <summary>
        /// 確定した辺を辿って鎖を作り、配列を書き出す<br/>
        /// 鎖の始まりは、入る辺の無い頂点 (入る辺を逆に辿って行き着く頂点) とする<br/>
        /// 環状の鎖は 2i (配列 i の順向き) から始め、戻る辺は使わない<br/>
        /// 使った配列は二度使わず、鎖に入らなかった配列は元のまま出す
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_頂点">向き付きの配列</param>
        /// <param name="p_出る">頂点ごとの確定した出る辺</param>
        /// <param name="p_入る">頂点ごとの入る辺の元の頂点、無ければ -1</param>
        /// <param name="p_繋いだ数">連結に使った辺の数</param>
        /// <returns>配列の出力 (元の番号の順)</returns>
        private static List<(string A_ID, string A_配列)> Get_鎖にして書き出す(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, string[] p_頂点, (int A_先, int A_左から削る長さ, string A_埋める配列, int A_右から削る長さ)?[] p_出る, int[] p_入る, out int p_繋いだ数)
        {
            var l_配列数 = p_配列群.Count;
            var l_使った = new bool[l_配列数];
            var l_結果 = new List<(int A_キー, string A_ID, string A_配列)>();
            p_繋いだ数 = 0;
            for (var i = 0; i < l_配列数; i++)
            {
                if (l_使った[i])
                {
                    continue;
                }

                var l_頭 = 2 * i;
                for (var l_歩 = 0; l_歩 <= p_頂点.Length && p_入る[l_頭] >= 0; l_歩++)
                {
                    var l_前 = p_入る[l_頭];
                    if (l_前 == 2 * i)
                    {
                        l_頭 = 2 * i;
                        break;
                    }

                    if (l_使った[l_前 / 2] || l_前 / 2 == i)
                    {
                        break;
                    }

                    l_頭 = l_前;
                }

                var l_鎖 = new StringBuilder(p_頂点[l_頭]);
                var l_今 = l_頭;
                var l_前からの削る長さ = 0;
                var l_連結数 = 0;
                while (p_出る[l_今] is { } l_辺)
                {
                    if (l_前からの削る長さ + l_辺.A_左から削る長さ >= p_頂点[l_今].Length)
                    {
                        break;
                    }

                    var l_先 = l_辺.A_先;
                    if (l_先 / 2 == l_頭 / 2 || l_使った[l_先 / 2])
                    {
                        break;
                    }

                    if (l_連結数 == 0)
                    {
                        l_使った[l_頭 / 2] = true;
                    }

                    l_鎖.Length -= l_辺.A_左から削る長さ;
                    _ = l_鎖.Append(l_辺.A_埋める配列);
                    _ = l_鎖.Append(p_頂点[l_先], l_辺.A_右から削る長さ, p_頂点[l_先].Length - l_辺.A_右から削る長さ);
                    l_使った[l_先 / 2] = true;
                    l_連結数++;
                    l_前からの削る長さ = l_辺.A_右から削る長さ;
                    l_今 = l_先;
                }

                if (l_連結数 == 0)
                {
                    continue;
                }

                p_繋いだ数 += l_連結数;
                l_結果.Add((l_頭 / 2, p_配列群[l_頭 / 2].A_ID, l_鎖.ToString()));
            }

            for (var i = 0; i < l_配列数; i++)
            {
                if (!l_使った[i])
                {
                    l_結果.Add((i, p_配列群[i].A_ID, p_配列群[i].A_配列));
                }
            }

            return [.. l_結果.OrderBy(x => x.A_キー).Select(x => (x.A_ID, x.A_配列))];
        }

        #endregion
    }
}
