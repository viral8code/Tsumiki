using System.Globalization;
using System.Text;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evidence;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// 全リードのペアで端どうしの相手を決め、その間を局所組み立てで一意に埋められたときだけ繋ぐ<br/>
    /// 配列 i の順向きを頂点 2i、逆相補を頂点 2i+1 とする<br/>
    /// 辺 u→x は、u の右端の後ろに x を (x の向きのまま) 続けること
    /// </summary>
    internal static class PairBridgeFiller
    {
        #region 定数

        /// <summary>
        /// 対象とする配列の最短長
        /// </summary>
        private const int C_対象の最短長 = 500;

        /// <summary>
        /// 種を取る、頂点の右端の窓の長さ (配列の末尾から)、および票を数える距離の上限
        /// </summary>
        private const int C_窓長 = 1000;

        /// <summary>
        /// 種の k-mer の長さ
        /// </summary>
        private const int C_種長 = 25;

        /// <summary>
        /// 挿入長を決めるのに要る標本数の下限
        /// </summary>
        private const int C_挿入長の最少標本数 = 100;

        /// <summary>
        /// 挿入長の標本として認める上限
        /// </summary>
        private const int C_挿入長の上限 = 2000;

        /// <summary>
        /// 相手を決める票数の下限
        /// </summary>
        private const int C_要るペア数 = 5;

        /// <summary>
        /// 首位の票が次点の何倍以上か
        /// </summary>
        private const int C_首位の優勢比 = 3;

        /// <summary>
        /// 局所組み立てに使う、端のアンカーの長さ
        /// </summary>
        private const int C_アンカー長 = 300;

        /// <summary>
        /// 局所組み立ての文脈長
        /// </summary>
        private const int C_局所k = 63;

        /// <summary>
        /// 局所組み立ての長さの余裕幅
        /// </summary>
        private const int C_長さの余裕幅 = 100;

        /// <summary>
        /// 1 つの辺に集める局所リードの上限
        /// </summary>
        private const int C_局所リード数の上限 = 4_000;

        /// <summary>
        /// 埋められない繋ぎ目を未確認の印 (N) で繋ぐのに要る、首位の票数の下限<br/>
        /// 埋めて繋ぐより強い証拠があるときだけにする
        /// </summary>
        private const int C_印で繋ぐのに要るペア数 = 10;

        /// <summary>
        /// 未確認の繋ぎ目を N 列で繋ぐときの、N の最小数
        /// </summary>
        private const int C_N列の最小数 = 10;

        /// <summary>
        /// 反復の端を調べる長さ (u の末尾、x の先頭)
        /// </summary>
        private const int C_反復を調べる長さ = 500;

        /// <summary>
        /// 反復の端と見なすまでに、両方に共有してよい 25-mer の数
        /// </summary>
        private const int C_共有を許す数 = 2;

        /// <summary>
        /// 配列の一致による重なりを探すときの、u の末尾の窓の長さ
        /// </summary>
        private const int C_重なりを探す窓長 = 301;

        /// <summary>
        /// 配列の一致による重なりを探すときの k 長 (重なりの上限は k - 1)
        /// </summary>
        private const int C_重なりを探すk長 = 135;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 全リードのペアで端どうしを結び、局所組み立てで間を埋めて配列を繋ぐ
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_ライブラリ群">ペアの順リードと逆リードのファイルのパス (逆リードが空のものは飛ばす)</param>
        /// <param name="p_統計">橋渡しの統計</param>
        /// <returns>
        /// 繋いだ後の配列の ID と配列<br/>
        /// 繋がらなかった配列は元の順と ID のまま
        /// </returns>
        public static List<(string A_ID, string A_配列)> Get_繋いだ配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IReadOnlyList<(string A_順リード, string A_逆リード)> p_ライブラリ群, out 橋渡し統計 p_統計)
        {
            var l_頂点 = Get_頂点群(p_配列群, out var l_対象);
            var l_表 = Get_種の索引(l_頂点, l_対象);
            var l_並列数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);

            var (l_標本, l_票) = Get_1回目の集計(l_頂点, l_表, p_ライブラリ群, l_並列数);
            if (l_標本.Count < C_挿入長の最少標本数)
            {
                p_統計 = default;
                return p_配列群.ToList();
            }

            var l_挿入長 = StatsUtil.Get_中央値(l_標本);
            var l_候補全体 = Get_候補の辺(l_票, l_挿入長);
            if (l_候補全体.Count == 0)
            {
                p_統計 = new 橋渡し統計(l_挿入長, 0, 0, 0, 0, 0, 0, 0);
                return p_配列群.ToList();
            }

            var l_候補 = new List<候補の辺>();
            var l_重なり群 = new List<int?>();
            var l_反復の端で除いた数 = 0;
            foreach (var l_辺 in l_候補全体)
            {
                var l_u文字列 = l_頂点[l_辺.A_u];
                var l_x文字列 = l_頂点[l_辺.A_x];
                var l_重なり = Get_重なりの長さ(l_u文字列, l_x文字列);
                if (Is反復の端(l_u文字列, l_x文字列, l_重なり))
                {
                    l_反復の端で除いた数++;
                    continue;
                }

                l_候補.Add(l_辺);
                l_重なり群.Add(l_重なり);
            }

            if (l_候補.Count == 0)
            {
                p_統計 = new 橋渡し統計(l_挿入長, l_候補全体.Count, 0, 0, 0, 0, l_反復の端で除いた数, 0);
                return p_配列群.ToList();
            }

            var l_読み = Get_局所リード群(l_頂点, l_表, p_ライブラリ群, l_候補, l_並列数);
            var (l_結果, l_判定) = Get_局所結果群(l_頂点, l_候補, l_読み);

            var l_出る = new (int A_先, int A_左から削る長さ, string A_埋める配列, int A_右から削る長さ)?[l_頂点.Length];
            var l_入る = new int[l_頂点.Length];
            Array.Fill(l_入る, -1);
            var l_埋めた数 = 0;
            var l_一意でない数 = 0;
            var l_届かない数 = 0;
            var l_印の数 = 0;
            var l_重なりで畳んだ数 = 0;
            for (var i = 0; i < l_候補.Count; i++)
            {
                var l_u = l_候補[i].A_u;
                var l_x = l_候補[i].A_x;
                var l_埋める配列 = l_結果[i];
                if (l_埋める配列 is not null)
                {
                    l_埋めた数++;
                    l_出る[l_u] = (l_x, 0, l_埋める配列, 0);
                    l_出る[l_x ^ 1] = (l_u ^ 1, 0, Util.V_逆相補_曖昧塩基あり(l_埋める配列), 0);
                    l_入る[l_x] = l_u;
                    l_入る[l_u ^ 1] = l_x ^ 1;
                    continue;
                }

                if (l_判定[i] == ギャップ充填判定.一意でない)
                {
                    l_一意でない数++;
                }
                else
                {
                    l_届かない数++;
                }

                if (l_重なり群[i] is { } l_重なり長)
                {
                    l_重なりで畳んだ数++;
                    l_出る[l_u] = (l_x, 0, "", l_重なり長);
                    l_出る[l_x ^ 1] = (l_u ^ 1, l_重なり長, "", 0);
                    l_入る[l_x] = l_u;
                    l_入る[l_u ^ 1] = l_x ^ 1;
                    continue;
                }

                if (l_候補[i].A_票数 >= C_印で繋ぐのに要るペア数)
                {
                    var l_N列 = new string('N', Math.Max(C_N列の最小数, l_候補[i].A_間の長さ));
                    l_印の数++;
                    l_出る[l_u] = (l_x, 0, l_N列, 0);
                    l_出る[l_x ^ 1] = (l_u ^ 1, 0, l_N列, 0);
                    l_入る[l_x] = l_u;
                    l_入る[l_u ^ 1] = l_x ^ 1;
                }
            }

            var l_書き出し = ReadEndJoiner.Get_鎖にして書き出す(p_配列群, l_頂点, l_出る, l_入る, out var l_繋いだ数);
            p_統計 = new 橋渡し統計(l_挿入長, l_候補全体.Count, l_埋めた数, l_一意でない数, l_届かない数, l_印の数, l_反復の端で除いた数, l_重なりで畳んだ数);
            return l_書き出し;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 配列から頂点の文字列を作る (2i が順向き、2i+1 が逆相補)
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_対象">配列ごとに、橋渡しの対象かどうか</param>
        /// <returns>頂点の文字列</returns>
        internal static string[] Get_頂点群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, out bool[] p_対象)
        {
            var l_配列数 = p_配列群.Count;
            var l_頂点 = new string[2 * l_配列数];
            p_対象 = new bool[l_配列数];
            for (var i = 0; i < l_配列数; i++)
            {
                l_頂点[2 * i] = p_配列群[i].A_配列;
                l_頂点[(2 * i) + 1] = Util.V_逆相補_曖昧塩基あり(p_配列群[i].A_配列);
                p_対象[i] = p_配列群[i].A_配列.Length >= C_対象の最短長 && !p_配列群[i].A_ID.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase);
            }

            return l_頂点;
        }

        /// <summary>
        /// 対象の頂点の右端の窓から、25-mer の正準形を種として索引に入れる<br/>
        /// 全配列 (順向き) の中で 2 回以上現れる正準形は捨てる
        /// </summary>
        /// <param name="p_頂点">頂点の文字列</param>
        /// <param name="p_対象">配列ごとに、橋渡しの対象かどうか</param>
        /// <returns>正準形から (頂点, 位置, 向きの印) への表</returns>
        internal static Dictionary<ulong, (int A_頂点, int A_位置, bool A_順向き)> Get_種の索引(string[] p_頂点, bool[] p_対象)
        {
            var l_表 = new Dictionary<ulong, (int A_頂点, int A_位置, bool A_順向き)>();
            for (var v = 0; v < p_頂点.Length; v++)
            {
                if (!p_対象[v / 2])
                {
                    continue;
                }

                var l_配列 = p_頂点[v];
                var l_始め = Math.Max(0, l_配列.Length - C_窓長);
                for (var p = l_始め; p + C_種長 <= l_配列.Length; p++)
                {
                    if (Is成功_キー(l_配列.AsSpan(p, C_種長), out var l_キー, out var l_順向き))
                    {
                        _ = l_表.TryAdd(l_キー, (v, p, l_順向き));
                    }
                }
            }

            var l_数 = new Dictionary<ulong, int>(l_表.Count);
            foreach (var l_キー in l_表.Keys)
            {
                l_数[l_キー] = 0;
            }

            for (var i = 0; i < p_頂点.Length; i += 2)
            {
                V_数える(p_頂点[i], l_数);
            }

            foreach (var l_キー in l_数.Where(x => x.Value >= 2).Select(x => x.Key).ToList())
            {
                _ = l_表.Remove(l_キー);
            }

            return l_表;
        }

        /// <summary>
        /// リードの 25-mer を先頭から見て、最初に索引に当たったものから、リードの 5' 端の位置を求める
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_表">種の索引</param>
        /// <param name="p_頂点">頂点の文字列</param>
        /// <returns>
        /// (頂点 u, 頂点 u の文字列でのリードの 5' 端の位置 q)<br/>
        /// q は負や長さ以上になることもある。当たりが無ければ null
        /// </returns>
        internal static (int A_頂点, int A_位置)? Get_当たり(string p_リード, IReadOnlyDictionary<ulong, (int A_頂点, int A_位置, bool A_順向き)> p_表, string[] p_頂点)
        {
            const ulong l_マスク = (1UL << (2 * C_種長)) - 1UL;
            var l_順 = 0UL;
            var l_逆 = 0UL;
            var l_連続 = 0;
            for (var i = 0; i < p_リード.Length; i++)
            {
                var l_値 = ReadEndJoiner.Get_塩基値(p_リード[i]);
                if (l_値 < 0)
                {
                    l_連続 = 0;
                    l_順 = 0UL;
                    l_逆 = 0UL;
                    continue;
                }

                l_順 = ((l_順 << 2) | (uint)l_値) & l_マスク;
                l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_値) << (2 * (C_種長 - 1)));
                l_連続 = Math.Min(l_連続 + 1, C_種長);
                if (l_連続 < C_種長)
                {
                    continue;
                }

                var l_正準 = Math.Min(l_順, l_逆);
                if (!p_表.TryGetValue(l_正準, out var l_当たり))
                {
                    continue;
                }

                var l_r = i - C_種長 + 1;
                var l_読みの向き = l_順 == l_正準;
                if (l_読みの向き == l_当たり.A_順向き)
                {
                    return (l_当たり.A_頂点, l_当たり.A_位置 - l_r);
                }

                var l_長さ = p_頂点[l_当たり.A_頂点].Length;
                return (l_当たり.A_頂点 ^ 1, l_長さ - l_当たり.A_位置 - C_種長 - l_r);
            }

            return null;
        }

        /// <summary>
        /// 向き付きの辺 (a→b) を正準化する。a→b と b^1→a^1 は同じ繋ぎ
        /// </summary>
        /// <param name="p_a">始点</param>
        /// <param name="p_b">終点</param>
        /// <returns>正準化した (始点, 終点)</returns>
        private static (int A_始, int A_終) Get_正準辺(int p_a, int p_b)
        {
            return p_a < (p_b ^ 1) || (p_a == (p_b ^ 1) && p_b <= (p_a ^ 1)) ? (p_a, p_b) : (p_b ^ 1, p_a ^ 1);
        }

        /// <summary>
        /// 1 回目の走査で、挿入長の標本と、辺の票 (正準辺ごとの d1 + d2 の一覧) を集める<br/>
        /// 並列に計算し、集計は読んだ順に行う
        /// </summary>
        /// <param name="p_頂点">頂点の文字列</param>
        /// <param name="p_表">種の索引</param>
        /// <param name="p_ライブラリ群">ペアのファイル</param>
        /// <param name="p_並列数">並列数</param>
        /// <returns>挿入長の標本と、正準辺ごとの票</returns>
        private static (List<int> 標本, Dictionary<(int A_始, int A_終), List<int>> 票) Get_1回目の集計(string[] p_頂点, IReadOnlyDictionary<ulong, (int A_頂点, int A_位置, bool A_順向き)> p_表, IReadOnlyList<(string A_順リード, string A_逆リード)> p_ライブラリ群, int p_並列数)
        {
            var l_標本 = new List<int>();
            var l_票 = new Dictionary<(int A_始, int A_終), List<int>>();
            var l_結果群 = Get_全ペア(p_ライブラリ群)
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(p_並列数)
                .Select(l_ペア => Get_組の結果(l_ペア.A_順リード, l_ペア.A_逆リード, p_頂点, p_表));
            foreach (var l_結果 in l_結果群)
            {
                if (l_結果.A_標本あり)
                {
                    l_標本.Add(l_結果.A_挿入長);
                }

                if (!l_結果.A_票あり)
                {
                    continue;
                }

                var l_鍵 = (l_結果.A_始, l_結果.A_終);
                if (!l_票.TryGetValue(l_鍵, out var l_一覧))
                {
                    l_一覧 = [];
                    l_票[l_鍵] = l_一覧;
                }

                l_一覧.Add(l_結果.A_和);
            }

            return (l_標本, l_票);
        }

        /// <summary>
        /// 1 組のペアから、挿入長の標本または辺の票を求める
        /// </summary>
        /// <param name="p_順リード">順リード</param>
        /// <param name="p_逆リード">逆リード</param>
        /// <param name="p_頂点">頂点の文字列</param>
        /// <param name="p_表">種の索引</param>
        /// <returns>組の結果 (どちらも無ければ既定値)</returns>
        private static 組の結果 Get_組の結果(string p_順リード, string p_逆リード, string[] p_頂点, IReadOnlyDictionary<ulong, (int A_頂点, int A_位置, bool A_順向き)> p_表)
        {
            var l_当たり1 = Get_当たり(p_順リード, p_表, p_頂点);
            var l_当たり2 = Get_当たり(p_逆リード, p_表, p_頂点);
            if (l_当たり1 is null || l_当たり2 is null)
            {
                return default;
            }

            var l_u1 = l_当たり1.Value.A_頂点;
            var l_u2 = l_当たり2.Value.A_頂点;
            var l_d1 = p_頂点[l_u1].Length - l_当たり1.Value.A_位置;
            var l_d2 = p_頂点[l_u2].Length - l_当たり2.Value.A_位置;
            if (l_u2 == (l_u1 ^ 1))
            {
                var l_断片長 = l_d1 + l_d2 - p_頂点[l_u1].Length;
                return new 組の結果(l_断片長 >= 1 && l_断片長 <= C_挿入長の上限, l_断片長, false, 0, 0, 0);
            }

            if (l_u1 / 2 != l_u2 / 2 && l_d1 > 0 && l_d1 <= C_窓長 && l_d2 > 0 && l_d2 <= C_窓長)
            {
                var (l_始, l_終) = Get_正準辺(l_u1, l_u2 ^ 1);
                return new 組の結果(false, 0, true, l_始, l_終, l_d1 + l_d2);
            }

            return default;
        }

        /// <summary>
        /// 全ライブラリのペアを、逆リードが空のものを飛ばして順に返す
        /// </summary>
        /// <param name="p_ライブラリ群">ペアのファイル</param>
        /// <returns>順リードと逆リードの組</returns>
        private static IEnumerable<(string A_順リード, string A_逆リード)> Get_全ペア(IReadOnlyList<(string A_順リード, string A_逆リード)> p_ライブラリ群)
        {
            foreach (var l_ライブラリ in p_ライブラリ群)
            {
                if (string.IsNullOrEmpty(l_ライブラリ.A_逆リード))
                {
                    continue;
                }

                foreach (var l_ペア in FastqReader.Get_ペア塩基列(l_ライブラリ.A_順リード, l_ライブラリ.A_逆リード))
                {
                    yield return l_ペア;
                }
            }
        }

        /// <summary>
        /// 票から、相互に一意で間の見込みが範囲内の候補の辺を選ぶ
        /// </summary>
        /// <param name="p_票">正準辺ごとの票 (d1 + d2 の一覧)</param>
        /// <param name="p_挿入長">挿入長の中央値</param>
        /// <returns>候補の辺 (u→x、u は u &lt; x^1 で一度だけ)</returns>
        private static List<候補の辺> Get_候補の辺(Dictionary<(int A_始, int A_終), List<int>> p_票, int p_挿入長)
        {
            var l_向き = new Dictionary<(int A_始, int A_終), List<int>>();
            foreach (var (l_鍵, l_一覧) in p_票)
            {
                l_向き[(l_鍵.A_始, l_鍵.A_終)] = l_一覧;
                l_向き[(l_鍵.A_終 ^ 1, l_鍵.A_始 ^ 1)] = l_一覧;
            }

            var l_首位 = new Dictionary<int, (int A_x, List<int> A_一覧)>();
            foreach (var l_群 in l_向き.GroupBy(x => x.Key.A_始))
            {
                var l_並び = l_群.OrderByDescending(x => x.Value.Count).ThenBy(x => x.Key.A_終).ToList();
                var l_首位の票 = l_並び[0].Value.Count;
                var l_次点の票 = l_並び.Count > 1 ? l_並び[1].Value.Count : 0;
                if (l_首位の票 >= C_要るペア数 && l_首位の票 >= C_首位の優勢比 * l_次点の票)
                {
                    l_首位[l_群.Key] = (l_並び[0].Key.A_終, l_並び[0].Value);
                }
            }

            var l_候補 = new List<候補の辺>();
            foreach (var l_u in l_首位.Keys.OrderBy(x => x))
            {
                var l_x = l_首位[l_u].A_x;
                if (l_u >= (l_x ^ 1) || l_u / 2 == l_x / 2)
                {
                    continue;
                }

                if (!l_首位.TryGetValue(l_x ^ 1, out var l_逆) || l_逆.A_x != (l_u ^ 1))
                {
                    continue;
                }

                var l_間の見込み = p_挿入長 - StatsUtil.Get_中央値(l_首位[l_u].A_一覧);
                if (l_間の見込み > Consts.ギャップ充填のギャップ長上限)
                {
                    continue;
                }

                l_候補.Add(new 候補の辺(l_u, l_x, l_間の見込み, l_首位[l_u].A_一覧.Count));
            }

            return l_候補;
        }

        /// <summary>
        /// 2 回目の走査で、候補の辺ごとに、端に当たったペアの両方のリードを集める<br/>
        /// 1 つの辺につき C_局所リード数の上限 本まで、読んだ順に集める
        /// </summary>
        /// <param name="p_頂点">頂点の文字列</param>
        /// <param name="p_表">種の索引</param>
        /// <param name="p_ライブラリ群">ペアのファイル</param>
        /// <param name="p_候補">候補の辺</param>
        /// <param name="p_並列数">並列数</param>
        /// <returns>候補の辺ごとの局所リード</returns>
        private static List<読取証拠>[] Get_局所リード群(string[] p_頂点, IReadOnlyDictionary<ulong, (int A_頂点, int A_位置, bool A_順向き)> p_表, IReadOnlyList<(string A_順リード, string A_逆リード)> p_ライブラリ群, List<候補の辺> p_候補, int p_並列数)
        {
            var l_候補の鍵 = new Dictionary<(int A_始, int A_終), int>();
            for (var i = 0; i < p_候補.Count; i++)
            {
                l_候補の鍵[(p_候補[i].A_u, p_候補[i].A_x ^ 1)] = i;
            }

            var l_読み = new List<読取証拠>[p_候補.Count];
            for (var i = 0; i < l_読み.Length; i++)
            {
                l_読み[i] = [];
            }

            var l_結果群 = Get_全ペア(p_ライブラリ群)
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(p_並列数)
                .Select((l_ペア, l_番号) => Get_該当の読み(l_ペア.A_順リード, l_ペア.A_逆リード, l_番号, p_頂点, p_表, l_候補の鍵));
            foreach (var l_群 in l_結果群)
            {
                foreach (var (l_候補番号, l_読取) in l_群)
                {
                    if (l_読み[l_候補番号].Count < C_局所リード数の上限)
                    {
                        l_読み[l_候補番号].Add(l_読取);
                    }
                }
            }

            return l_読み;
        }

        /// <summary>
        /// ペアの 2 本のリードが、候補の辺の u と x^1 の端に当たるかを調べ、当たれば局所リードを返す<br/>
        /// u に当たったリードはそのまま、x^1 に当たったリードは逆相補にして、u→x の向きで揃える
        /// </summary>
        /// <param name="p_順リード">順リード</param>
        /// <param name="p_逆リード">逆リード</param>
        /// <param name="p_番号">ペアの通し番号 (ペア ID になる)</param>
        /// <param name="p_頂点">頂点の文字列</param>
        /// <param name="p_表">種の索引</param>
        /// <param name="p_候補の鍵">(u, x^1) から候補の辺の番号への表</param>
        /// <returns>候補の辺の番号と局所リードの組</returns>
        private static (int A_候補, 読取証拠 A_読取)[] Get_該当の読み(string p_順リード, string p_逆リード, int p_番号, string[] p_頂点, IReadOnlyDictionary<ulong, (int A_頂点, int A_位置, bool A_順向き)> p_表, IReadOnlyDictionary<(int A_始, int A_終), int> p_候補の鍵)
        {
            var l_当たり1 = Get_当たり(p_順リード, p_表, p_頂点);
            var l_当たり2 = Get_当たり(p_逆リード, p_表, p_頂点);
            if (l_当たり1 is null || l_当たり2 is null)
            {
                return [];
            }

            var l_h1 = l_当たり1.Value.A_頂点;
            var l_h2 = l_当たり2.Value.A_頂点;
            var l_d1 = p_頂点[l_h1].Length - l_当たり1.Value.A_位置;
            var l_d2 = p_頂点[l_h2].Length - l_当たり2.Value.A_位置;
            if (l_d1 <= 0 || l_d1 > C_窓長 || l_d2 <= 0 || l_d2 > C_窓長)
            {
                return [];
            }

            var l_ペアID = p_番号.ToString(CultureInfo.InvariantCulture);
            if (p_候補の鍵.TryGetValue((l_h1, l_h2), out var l_候補1))
            {
                return
                [
                    (l_候補1, new 読取証拠(p_順リード, l_ペアID)),
                    (l_候補1, new 読取証拠(Util.V_逆相補_曖昧塩基あり(p_逆リード), l_ペアID)),
                ];
            }

            if (p_候補の鍵.TryGetValue((l_h2, l_h1), out var l_候補2))
            {
                return
                [
                    (l_候補2, new 読取証拠(Util.V_逆相補_曖昧塩基あり(p_順リード), l_ペアID)),
                    (l_候補2, new 読取証拠(p_逆リード, l_ペアID)),
                ];
            }

            return [];
        }

        /// <summary>
        /// 候補の辺ごとに、アンカーの間を局所組み立てで埋める<br/>
        /// 局所リードが無い辺は、埋めずに到達不能とする
        /// </summary>
        /// <param name="p_頂点">頂点の文字列</param>
        /// <param name="p_候補">候補の辺</param>
        /// <param name="p_読み">候補の辺ごとの局所リード</param>
        /// <returns>埋める配列 (埋まらなければ null) と判定の配列</returns>
        private static (string?[] 結果, ギャップ充填判定[] 判定) Get_局所結果群(string[] p_頂点, List<候補の辺> p_候補, List<読取証拠>[] p_読み)
        {
            var l_結果 = new string?[p_候補.Count];
            var l_判定 = new ギャップ充填判定[p_候補.Count];
            var l_並列設定 = new ParallelOptions
            {
                MaxDegreeOfParallelism = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数)
            };
            _ = Parallel.For(0, p_候補.Count, l_並列設定, i =>
            {
                if (p_読み[i].Count == 0)
                {
                    l_判定[i] = ギャップ充填判定.到達不能;
                    return;
                }

                var l_u文字列 = p_頂点[p_候補[i].A_u];
                var l_x文字列 = p_頂点[p_候補[i].A_x];
                var l_左 = l_u文字列.Substring(Math.Max(0, l_u文字列.Length - C_アンカー長));
                var l_右 = l_x文字列[..Math.Min(C_アンカー長, l_x文字列.Length)];
                l_結果[i] = LocalAssembler.Get_適応kの局所結果(new 局所ギャップ(0, 0, Math.Max(0, p_候補[i].A_間の長さ), l_左, l_右), p_読み[i], C_局所k, out var l_局所判定, C_長さの余裕幅);
                l_判定[i] = l_局所判定;
            });
            return (l_結果, l_判定);
        }

        /// <summary>
        /// u の末尾と x の先頭が完全に一致する重なりの長さを、配列の一致で求める (u の末尾の窓は C_重なりを探す窓長)
        /// </summary>
        /// <param name="p_u文字列">u の文字列</param>
        /// <param name="p_x文字列">x の文字列</param>
        /// <returns>重なりの長さ、決められなければ null</returns>
        private static int? Get_重なりの長さ(string p_u文字列, string p_x文字列)
        {
            var l_窓 = p_u文字列[Math.Max(0, p_u文字列.Length - C_重なりを探す窓長)..];
            return Scaffolder.Get_配列の一致で決めた重なり長(new StringBuilder(l_窓), p_x文字列, C_重なりを探すk長);
        }

        /// <summary>
        /// 候補の辺の両端が反復の端かを調べる<br/>
        /// u の末尾 C_反復を調べる長さ (重なり o があれば最後の o bp を除く) と、x の先頭 C_反復を調べる長さ で、
        /// 両方に現れる正準の 25-mer が C_共有を許す数 を超えれば反復の端とする
        /// </summary>
        /// <param name="p_u文字列">u の文字列</param>
        /// <param name="p_x文字列">x の文字列</param>
        /// <param name="p_重なり">u の末尾と x の先頭の重なりの長さ、無ければ null</param>
        /// <returns>反復の端なら true</returns>
        private static bool Is反復の端(string p_u文字列, string p_x文字列, int? p_重なり)
        {
            var l_u末尾 = p_u文字列[Math.Max(0, p_u文字列.Length - C_反復を調べる長さ)..];
            var l_除く = Math.Min(p_重なり ?? 0, l_u末尾.Length);
            l_u末尾 = l_u末尾[..(l_u末尾.Length - l_除く)];
            var l_x先頭 = p_x文字列[..Math.Min(C_反復を調べる長さ, p_x文字列.Length)];
            var l_x種 = Get_種の集合(l_x先頭);
            var l_共有 = 0;
            foreach (var l_キー in Get_種の集合(l_u末尾))
            {
                if (l_x種.Contains(l_キー))
                {
                    l_共有++;
                }
            }

            return l_共有 > C_共有を許す数;
        }

        /// <summary>
        /// 配列に現れる 25-mer の正準形の集合を求める
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <returns>正準形の集合 (A/C/G/T 以外を含む 25-mer は除く)</returns>
        private static HashSet<ulong> Get_種の集合(string p_配列)
        {
            var l_集合 = new HashSet<ulong>();
            for (var p = 0; p + C_種長 <= p_配列.Length; p++)
            {
                if (Is成功_キー(p_配列.AsSpan(p, C_種長), out var l_キー, out _))
                {
                    _ = l_集合.Add(l_キー);
                }
            }

            return l_集合;
        }

        /// <summary>
        /// 配列の 25-mer の正準形と、その向きの印を求める
        /// </summary>
        /// <param name="p_断片">25-mer</param>
        /// <param name="p_正準">正準形の 2 ビット表現</param>
        /// <param name="p_順向き">断片そのものが正準形か</param>
        /// <returns>求めたか (A/C/G/T 以外を含めば false)</returns>
        private static bool Is成功_キー(ReadOnlySpan<char> p_断片, out ulong p_正準, out bool p_順向き)
        {
            var l_順 = 0UL;
            var l_逆 = 0UL;
            foreach (var l_塩基 in p_断片)
            {
                var l_値 = ReadEndJoiner.Get_塩基値(l_塩基);
                if (l_値 < 0)
                {
                    p_正準 = 0UL;
                    p_順向き = false;
                    return false;
                }

                l_順 = (l_順 << 2) | (uint)l_値;
                l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_値) << (2 * (C_種長 - 1)));
            }

            p_正準 = Math.Min(l_順, l_逆);
            p_順向き = l_順 == p_正準;
            return true;
        }

        /// <summary>
        /// 配列の 25-mer の正準形を、数えたい正準形の表 (数は 0 で入っている) に数え足す
        /// </summary>
        /// <param name="p_配列">順向きの配列 (対象外も含む)</param>
        /// <param name="p_数">正準形ごとの出現数 (表に無いものは数えない)</param>
        private static void V_数える(string p_配列, Dictionary<ulong, int> p_数)
        {
            const ulong l_マスク = (1UL << (2 * C_種長)) - 1UL;
            var l_順 = 0UL;
            var l_逆 = 0UL;
            var l_連続 = 0;
            foreach (var l_塩基 in p_配列)
            {
                var l_値 = ReadEndJoiner.Get_塩基値(l_塩基);
                if (l_値 < 0)
                {
                    l_連続 = 0;
                    l_順 = 0UL;
                    l_逆 = 0UL;
                    continue;
                }

                l_順 = ((l_順 << 2) | (uint)l_値) & l_マスク;
                l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_値) << (2 * (C_種長 - 1)));
                l_連続 = Math.Min(l_連続 + 1, C_種長);
                if (l_連続 < C_種長)
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
        /// 1 組のペアの結果 (挿入長の標本、または辺の票)
        /// </summary>
        /// <param name="A_標本あり">挿入長の標本があるか</param>
        /// <param name="A_挿入長">挿入長の標本</param>
        /// <param name="A_票あり">辺の票があるか</param>
        /// <param name="A_始">票の正準辺の始点</param>
        /// <param name="A_終">票の正準辺の終点</param>
        /// <param name="A_和">d1 + d2</param>
        private readonly record struct 組の結果(bool A_標本あり, int A_挿入長, bool A_票あり, int A_始, int A_終, int A_和);

        /// <summary>
        /// 候補の辺 (u→x) と、u と x の端の間の見込み
        /// </summary>
        /// <param name="A_u">始点の頂点</param>
        /// <param name="A_x">終点の頂点</param>
        /// <param name="A_間の長さ">間の見込み (挿入長 − d1 − d2 の中央値)</param>
        /// <param name="A_票数">首位の辺の票数 (ペアの数)</param>
        private readonly record struct 候補の辺(int A_u, int A_x, int A_間の長さ, int A_票数);

        #endregion
    }
}
