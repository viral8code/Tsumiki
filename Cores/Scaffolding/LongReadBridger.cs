using System.Text;
using Tsumiki.Commons;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Scaffolding
{
    /// <summary>
    /// ロングリードを最終配列 (錨) に当て、反復を跨ぐ互いに一意な端の組 (橋) を決める
    /// </summary>
    internal static class LongReadBridger
    {
        #region 定数

        /// <summary>
        /// 錨にする配列の長さの下限
        /// </summary>
        private const int C_錨の最短長 = 1000;

        /// <summary>
        /// minimizer の k-mer の長さ
        /// </summary>
        private const int C_種長 = 15;

        /// <summary>
        /// minimizer を選ぶ窓の k-mer の数
        /// </summary>
        private const int C_窓の種数 = 10;

        /// <summary>
        /// 索引に残す minimizer の、錨全体での出現数の上限
        /// </summary>
        private const int C_種の出現の上限 = 8;

        /// <summary>
        /// 当て方に使うリードの長さの下限
        /// </summary>
        private const int C_リードの最短長 = 1000;

        /// <summary>
        /// 鎖で隣り合う当たりの、リード上・錨上の差の上限
        /// </summary>
        private const int C_鎖の差の上限 = 5000;

        /// <summary>
        /// 鎖で隣り合う当たりの、リード上と錨上の差のずれを許す割合
        /// </summary>
        private const double C_鎖のずれの割合 = 0.2D;

        /// <summary>
        /// 鎖で隣り合う当たりの、リード上と錨上の差のずれを許す余裕 (bp)
        /// </summary>
        private const int C_鎖のずれの余裕 = 50;

        /// <summary>
        /// 鎖を作るときに遡る当たりの数
        /// </summary>
        private const int C_鎖を遡る数 = 50;

        /// <summary>
        /// 当たりとして採るのに要る鎖の minimizer の数
        /// </summary>
        private const int C_当たりに要る種の数 = 10;

        /// <summary>
        /// 当たりとして採るのに要る錨上の幅
        /// </summary>
        private const int C_当たりに要る錨上の幅 = 500;

        /// <summary>
        /// リード上で重なる当たりのうち、種の少ないほうを捨てる重なりの割合 (短いほうの当たりの長さに対して)
        /// </summary>
        private const double C_重なりで捨てる割合 = 0.5D;

        /// <summary>
        /// リード上で重なる 2 つの当たりの種の数の比がこれ以上なら、どちらとも決められないとして両方を捨てる
        /// </summary>
        private const double C_曖昧とみなす種の比 = 0.8D;

        /// <summary>
        /// 深さを数える区画の長さ
        /// </summary>
        private const int C_区画長 = 200;

        /// <summary>
        /// 深さが基準の何倍以上なら多コピー (反復) と見なすか
        /// </summary>
        private const double C_多コピーとみなす深さの比 = 1.5D;

        /// <summary>
        /// 多コピーかどうかを判定する錨の長さの上限 (これ以上の錨はプラスミドなどの単位とみて、多コピーと見なさない)
        /// </summary>
        private const int C_多コピーを判定する錨の長さの上限 = 10_000;

        /// <summary>
        /// 橋の証拠に使う当たりの、当たりの端から錨の端までの長さの上限 (リードが錨の端を跨いだとみなせる距離)
        /// </summary>
        private const int C_端までの残りの上限 = 1000;

        /// <summary>
        /// 当たりを橋の証拠に使うのに要る、丸ごと覆う反復でない区画の数
        /// </summary>
        private const int C_証拠に要る一意な区画の数 = 3;

        /// <summary>
        /// 相手を決めるのに要る票の下限
        /// </summary>
        private const int C_要る票 = 3;

        /// <summary>
        /// 首位の票が次点の何倍以上か
        /// </summary>
        private const int C_首位の優勢比 = 3;

        /// <summary>
        /// 間を埋めるときに許す幅の、間の長さに対する割合
        /// </summary>
        private const double C_許す幅の割合 = 0.08D;

        /// <summary>
        /// 間を埋めるときに許す幅の余裕 (bp)
        /// </summary>
        private const int C_許す幅の余裕 = 60;

        /// <summary>
        /// 重なりで畳むときの重なりの長さの下限
        /// </summary>
        private const int C_重なりの最短長 = 20;

        /// <summary>
        /// 重なりの配列の先頭と末尾がこの長さ以上一致すれば (周期を持てば)、縦列反復の中の見かけの重なりとみて使わない
        /// </summary>
        private const int C_周期とみなす境界の長さ = 10;

        /// <summary>
        /// ロングリードの区間に含める、当たりの外側の長さ (u の当たりの終わりの手前、x の当たりの始まりの先)
        /// </summary>
        private const int C_区間の外側の長さ = 300;

        /// <summary>
        /// 橋ごとに残すロングリードの区間の数の上限
        /// </summary>
        private const int C_区間の数の上限 = 30;

        /// <summary>
        /// 道を選ぶのに使える区間の数の下限
        /// </summary>
        private const int C_道選びに要る区間の数 = 3;

        /// <summary>
        /// 歩き始めの位置を決めるときに、区間に整列する u の末尾の長さ
        /// </summary>
        private const int C_歩き始めの整列長 = 60;

        /// <summary>
        /// 歩き始めの整列で許す編集距離の割合
        /// </summary>
        private const double C_歩き始めの誤りの割合 = 0.25D;

        /// <summary>
        /// 分かれ道の先読みの長さ
        /// </summary>
        private const int C_先読み長 = 24;

        /// <summary>
        /// 分かれ道の候補ごとの先読みの配列の数の上限
        /// </summary>
        private const int C_先読みの数の上限 = 8;

        /// <summary>
        /// 整列の窓を、整列する長さの何倍にするか
        /// </summary>
        private const double C_窓の倍率 = 1.3D;

        /// <summary>
        /// 整列の窓に足す余裕 (bp)
        /// </summary>
        private const int C_窓の余裕 = 10;

        /// <summary>
        /// 未整列の塩基がこの長さに達したら、区間の現在位置を進める
        /// </summary>
        private const int C_まとめて整列する長さ = 100;

        /// <summary>
        /// 区間を使い続けるのに許す、整列の編集距離の割合
        /// </summary>
        private const double C_歩きの誤りの割合 = 0.3D;

        /// <summary>
        /// 分かれ道で、首位の票が次点の何倍以上か
        /// </summary>
        private const int C_道選びの優勢比 = 2;

        /// <summary>
        /// 道を選んで埋めたとき、見込みの長さから許す幅の倍率 (L2 の幅に対して)
        /// </summary>
        private const int C_道選びの幅の倍率 = 2;

        /// <summary>
        /// 埋めた配列を直すのに要る、数えた区間の数の下限
        /// </summary>
        private const int C_直すのに要る区間の数 = 5;

        /// <summary>
        /// 埋めた配列を直すとき、別の塩基を示す区間の割合の下限
        /// </summary>
        private const double C_直す塩基の割合 = 0.7D;

        /// <summary>
        /// 埋めた配列を直すとき、元の塩基を示す区間の割合の上限
        /// </summary>
        private const double C_元の塩基の割合の上限 = 0.2D;

        /// <summary>
        /// 埋めた配列を直すとき、区間に整列するときに末尾に足す x の先頭の長さ
        /// </summary>
        private const int C_直すときに足す右の長さ = 60;

        /// <summary>
        /// 直す塩基を短いリードで確かめる k-mer の、中心から片側の長さ (k-mer の長さは 2 倍 + 1)
        /// </summary>
        public const int C_確かめるkmerの片側の長さ = 20;

        /// <summary>
        /// 区間の整列で欠失を表す文字
        /// </summary>
        private const char C_欠失の文字 = '-';

        /// <summary>
        /// 経路の材料にする GFA の segment の長さの上限
        /// </summary>
        public const int C_材料にするsegmentの長さの上限 = 20_000;

        /// <summary>
        /// 経路の材料に足す、橋の両端の配列の長さ (k に足す長さ)
        /// </summary>
        private const int C_材料にする端の長さ = 300;

        /// <summary>
        /// TSV の見出し
        /// </summary>
        private const string C_TSVの見出し = "left_id\tleft_strand\tright_id\tright_strand\tvotes\tsecond_votes\tgap_median\tjoin";

        /// <summary>
        /// 重なりで畳んだ印
        /// </summary>
        private const string C_重なりの印 = "overlap:";

        /// <summary>
        /// 経路で埋めた印
        /// </summary>
        private const string C_経路の印 = "path:";

        /// <summary>
        /// 繋がなかった印
        /// </summary>
        private const string C_繋がない印 = "none:";

        /// <summary>
        /// ロングリードで道を選んで埋めた印
        /// </summary>
        private const string C_道選びの印 = "guided:";

        /// <summary>
        /// 繋がなかった理由: 重なりが見つからない
        /// </summary>
        private const string C_理由_重なりなし = "no_overlap";

        /// <summary>
        /// 繋がなかった理由: 経路が一意でない
        /// </summary>
        private const string C_理由_一意でない = "ambiguous";

        /// <summary>
        /// 繋がなかった理由: 経路が見つからない
        /// </summary>
        private const string C_理由_到達不能 = "unreachable";

        /// <summary>
        /// 繋がなかった理由: 探索を打ち切った
        /// </summary>
        private const string C_理由_打切り = "aborted";

        /// <summary>
        /// 繋がなかった理由: 経路の材料 (GFA) が無い
        /// </summary>
        private const string C_理由_グラフなし = "no_graph";

        /// <summary>
        /// 順向きの印
        /// </summary>
        private const string C_順向きの印 = "+";

        /// <summary>
        /// 逆向きの印
        /// </summary>
        private const string C_逆向きの印 = "-";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// ロングリードを錨に当て、互いに一意な端の組 (橋) を決める
        /// </summary>
        /// <param name="p_配列群">最終配列の ID と配列</param>
        /// <param name="p_ロングリード群">ロングリードの塩基列</param>
        /// <param name="p_統計">橋の組を決めた結果の統計</param>
        /// <returns>決めた橋 (u→x と x^1→u^1 のうち一方だけ)</returns>
        public static List<ロングリードの橋> Get_決めた橋(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IEnumerable<string> p_ロングリード群, out ロングリード橋渡し統計 p_統計)
        {
            return Get_決めた橋(p_配列群, p_ロングリード群, false, out p_統計, out _);
        }

        /// <summary>
        /// ロングリードを錨に当て、互いに一意な端の組 (橋) を決める。求めれば、橋ごとにそれを跨いだロングリードの区間も返す<br/>
        /// 区間は、橋を決めたあとにロングリードをもう一度読んで取り出す
        /// </summary>
        /// <param name="p_配列群">最終配列の ID と配列</param>
        /// <param name="p_ロングリード群">ロングリードの塩基列 (2 回読めること)</param>
        /// <param name="p_Is区間も求める">区間も求めるか</param>
        /// <param name="p_統計">橋の組を決めた結果の統計</param>
        /// <param name="p_区間群">橋ごとの区間 (u→x の向き、橋と同じ順)、求めなければ空</param>
        /// <returns>決めた橋 (u→x と x^1→u^1 のうち一方だけ)</returns>
        public static List<ロングリードの橋> Get_決めた橋(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IEnumerable<string> p_ロングリード群, bool p_Is区間も求める, out ロングリード橋渡し統計 p_統計, out List<List<string>> p_区間群)
        {
            var l_Is錨 = p_配列群.Select(x => x.A_配列.Length >= C_錨の最短長 && !x.A_ID.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase)).ToArray();
            var l_索引 = Get_索引(p_配列群, l_Is錨);
            var l_長さ = p_配列群.Select(x => x.A_配列.Length).ToArray();
            var l_並列数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_リード数 = 0;
            var l_当たり群 = new List<List<ロングリードの当たり>>();
            var l_リード長群 = new List<int>();
            foreach (var (l_リード長, l_当たり) in p_ロングリード群.Where(x => x.Length >= C_リードの最短長).AsParallel().AsOrdered().WithDegreeOfParallelism(l_並列数).Select(x => (x.Length, Get_当たり群(x, l_索引, l_長さ))))
            {
                l_リード数++;
                if (l_当たり.Count > 0)
                {
                    l_当たり群.Add(l_当たり);
                    l_リード長群.Add(l_リード長);
                }
            }

            var l_区画の深さ = Get_区画の深さ(l_当たり群, l_長さ, l_Is錨);
            var l_基準 = Get_基準の深さ(l_区画の深さ);
            var l_錨の深さ = l_区画の深さ.Select(Get_錨の深さ).ToArray();
            var l_Is多コピー = l_錨の深さ.Select((x, i) => l_長さ[i] < C_多コピーを判定する錨の長さの上限 && x >= C_多コピーとみなす深さの比 * l_基準).ToArray();
            var l_票 = new Dictionary<(int A_u, int A_x), List<int>>();
            bool Is使える(ロングリードの当たり x) => !l_Is多コピー[x.A_錨] && Is証拠に使える(x, l_区画の深さ[x.A_錨], Math.Max(l_基準, l_錨の深さ[x.A_錨]));
            for (var r = 0; r < l_当たり群.Count; r++)
            {
                V_数える_証拠(l_当たり群[r].Where(Is使える).ToList(), l_リード長群[r], l_長さ, l_票);
            }

            var l_橋群 = Get_互いに一意な橋(l_票);
            p_統計 = new ロングリード橋渡し統計(l_リード数, l_当たり群.Sum(x => x.Count), l_Is錨.Count(x => x), l_Is多コピー.Count(x => x), l_票.Keys.Count(x => x == Get_正準辺(x.A_u, x.A_x)), l_橋群.Count);
            p_区間群 = p_Is区間も求める ? Get_区間群(p_ロングリード群, l_索引, l_長さ, Is使える, l_橋群, l_並列数) : [];
            return l_橋群;
        }

        /// <summary>
        /// 決めた橋の間を、重なりで畳むか、グラフの経路で埋めて繋ぐ
        /// </summary>
        /// <param name="p_配列群">最終配列の ID と配列</param>
        /// <param name="p_橋群">決めた橋</param>
        /// <param name="p_経路の材料">経路を探す k-mer を取る配列 (GFA の segment など)、無ければ null</param>
        /// <param name="p_区間群">橋ごとのロングリードの区間 (u→x の向き、橋と同じ順)。経路が一意でないときに道を選ぶのと、埋めた配列を直すのに使う。無ければ null</param>
        /// <param name="p_短いリードで確かめる">直す塩基を中心にした k-mer が短いリードで確かめられるか。無ければ直さない</param>
        /// <param name="p_k長">経路を探す k-mer の長さ (選ばれた k)</param>
        /// <param name="p_繋ぎ方群">橋ごとの埋め方 (p_橋群 と同じ順)</param>
        /// <param name="p_統計">繋いだ結果の統計</param>
        /// <returns>繋いだ後の配列の ID と配列 (繋がなかった配列は元のまま)</returns>
        public static List<(string A_ID, string A_配列)> Get_繋いだ配列群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IReadOnlyList<ロングリードの橋> p_橋群, IEnumerable<string>? p_経路の材料, IReadOnlyList<IReadOnlyList<string>>? p_区間群, Func<string, bool>? p_短いリードで確かめる, int p_k長, out List<ロングリードの繋ぎ方> p_繋ぎ方群, out ロングリード繋ぎ統計 p_統計)
        {
            var l_頂点 = new string[2 * p_配列群.Count];
            for (var i = 0; i < p_配列群.Count; i++)
            {
                l_頂点[2 * i] = p_配列群[i].A_配列;
                l_頂点[(2 * i) + 1] = Util.V_逆相補_曖昧塩基あり(p_配列群[i].A_配列);
            }

            var l_集合 = p_経路の材料 is null ? null : Get_経路の材料(p_経路の材料, p_橋群, l_頂点, p_k長);
            var l_出る = new (int A_先, int A_左から削る長さ, string A_埋める配列, int A_右から削る長さ)?[l_頂点.Length];
            var l_入る = new int[l_頂点.Length];
            Array.Fill(l_入る, -1);
            p_繋ぎ方群 = [];
            var l_直した塩基数 = 0;
            for (var b = 0; b < p_橋群.Count; b++)
            {
                var l_橋 = p_橋群[b];
                var l_繋ぎ方 = Get_繋ぎ方(l_頂点[l_橋.A_u], l_頂点[l_橋.A_x], l_橋.A_間の長さ, l_集合, p_k長);
                if (l_繋ぎ方.A_理由 is C_理由_一意でない or C_理由_到達不能 or C_理由_打切り && l_集合 is not null && p_区間群 is not null && b < p_区間群.Count
                    && Get_ロングリードで選んだ経路(l_頂点[l_橋.A_u], l_頂点[l_橋.A_x], l_橋.A_間の長さ, l_集合, p_k長, p_区間群[b]) is { } l_選んだ経路)
                {
                    l_繋ぎ方 = new ロングリードの繋ぎ方(null, l_選んだ経路, null, true);
                }

                if (l_繋ぎ方.A_埋める配列 is { } l_直す前 && p_短いリードで確かめる is not null && p_区間群 is not null && b < p_区間群.Count)
                {
                    l_繋ぎ方 = l_繋ぎ方 with { A_埋める配列 = Get_直した埋める配列(l_頂点[l_橋.A_u], l_直す前, l_頂点[l_橋.A_x], p_区間群[b], p_短いリードで確かめる, out var l_直した数) };
                    l_直した塩基数 += l_直した数;
                }

                p_繋ぎ方群.Add(l_繋ぎ方);
                if (l_繋ぎ方.A_重なり is { } l_重なり)
                {
                    l_出る[l_橋.A_u] = (l_橋.A_x, 0, string.Empty, l_重なり);
                    l_出る[l_橋.A_x ^ 1] = (l_橋.A_u ^ 1, l_重なり, string.Empty, 0);
                }
                else if (l_繋ぎ方.A_埋める配列 is { } l_埋める配列)
                {
                    l_出る[l_橋.A_u] = (l_橋.A_x, 0, l_埋める配列, 0);
                    l_出る[l_橋.A_x ^ 1] = (l_橋.A_u ^ 1, 0, Util.V_逆相補_曖昧塩基あり(l_埋める配列), 0);
                }
                else
                {
                    continue;
                }

                l_入る[l_橋.A_x] = l_橋.A_u;
                l_入る[l_橋.A_u ^ 1] = l_橋.A_x ^ 1;
            }

            var l_環状群 = Get_環状の鎖群(p_配列群, l_頂点, l_出る, l_入る);
            var l_結果 = ReadEndJoiner.Get_鎖にして書き出す(p_配列群, l_頂点, l_出る, l_入る, out _);
            if (l_環状群.Count > 0)
            {
                var l_外す = l_環状群.SelectMany(x => x.A_配列の番号群).Select(x => p_配列群[x].A_ID).ToHashSet(StringComparer.Ordinal);
                var l_環状の配列 = l_環状群.ToDictionary(x => p_配列群[x.A_配列の番号群[0]].A_ID, x => (x.A_ID, x.A_配列), StringComparer.Ordinal);
                l_結果 = [.. l_結果.Where(x => !l_外す.Contains(x.A_ID) || l_環状の配列.ContainsKey(x.A_ID)).Select(x => l_環状の配列.TryGetValue(x.A_ID, out var l_環) ? l_環 : x)];
            }

            p_統計 = new ロングリード繋ぎ統計(p_繋ぎ方群.Count(x => x.A_重なり is not null), p_繋ぎ方群.Count(x => x.A_埋める配列 is not null && !x.A_Isロングリードで選んだ), p_繋ぎ方群.Count(x => x.A_Isロングリードで選んだ), p_繋ぎ方群.Count(x => x.A_理由 is not null), p_配列群.Count, l_結果.Count, l_直した塩基数, l_環状群.Count);
            return l_結果;
        }

        /// <summary>
        /// 決めた橋を TSV に書き出す
        /// </summary>
        /// <param name="p_パス">出力先</param>
        /// <param name="p_配列群">最終配列の ID と配列</param>
        /// <param name="p_橋群">決めた橋</param>
        /// <param name="p_繋ぎ方群">橋ごとの埋め方 (p_橋群 と同じ順)、無ければ null</param>
        public static void V_書き出し_橋(string p_パス, IReadOnlyList<(string A_ID, string A_配列)> p_配列群, IReadOnlyList<ロングリードの橋> p_橋群, IReadOnlyList<ロングリードの繋ぎ方>? p_繋ぎ方群 = null)
        {
            using var l_書き込み = new StreamWriter(p_パス);
            l_書き込み.WriteLine(C_TSVの見出し);
            for (var i = 0; i < p_橋群.Count; i++)
            {
                var l_橋 = p_橋群[i];
                var l_繋ぎ = p_繋ぎ方群 is null ? string.Empty : Get_繋ぎ方の印(p_繋ぎ方群[i]);
                l_書き込み.WriteLine($"{p_配列群[l_橋.A_u / 2].A_ID}\t{Get_向きの印(l_橋.A_u)}\t{p_配列群[l_橋.A_x / 2].A_ID}\t{Get_向きの印(l_橋.A_x)}\t{l_橋.A_票}\t{l_橋.A_次点の票}\t{l_橋.A_間の長さ}\t{l_繋ぎ}");
            }
        }

        /// <summary>
        /// 配列の minimizer (k = C_種長、窓 C_窓の種数) を求める
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <returns>(正準形, 位置, 配列の向きの k-mer が正準形か) の並び (位置の順)</returns>
        public static List<(ulong A_キー, int A_位置, bool A_Is順向き)> Get_minimizer群(string p_配列)
        {
            var l_k数 = p_配列.Length - C_種長 + 1;
            var l_結果 = new List<(ulong A_キー, int A_位置, bool A_Is順向き)>();
            if (l_k数 < C_窓の種数)
            {
                return l_結果;
            }

            var l_ハッシュ = new ulong[l_k数];
            var l_キー = new ulong[l_k数];
            var l_Is順向き = new bool[l_k数];
            const ulong l_マスク = (1UL << (2 * C_種長)) - 1UL;
            var l_順 = 0UL;
            var l_逆 = 0UL;
            var l_連続 = 0;
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_値 = ReadEndJoiner.Get_塩基値(p_配列[i]);
                if (l_値 < 0)
                {
                    l_連続 = 0;
                    l_順 = 0UL;
                    l_逆 = 0UL;
                }
                else
                {
                    l_順 = ((l_順 << 2) | (uint)l_値) & l_マスク;
                    l_逆 = (l_逆 >> 2) | ((ulong)(3 - l_値) << (2 * (C_種長 - 1)));
                    l_連続 = Math.Min(l_連続 + 1, C_種長);
                }

                var l_始め = i - C_種長 + 1;
                if (l_始め < 0)
                {
                    continue;
                }

                if (l_連続 < C_種長)
                {
                    l_ハッシュ[l_始め] = ulong.MaxValue;
                    continue;
                }

                l_キー[l_始め] = Math.Min(l_順, l_逆);
                l_Is順向き[l_始め] = l_順 < l_逆;
                l_ハッシュ[l_始め] = Get_混ぜた値(l_キー[l_始め]);
            }

            var l_前 = -1;
            for (var s = 0; s + C_窓の種数 <= l_k数; s++)
            {
                var l_最小 = s;
                for (var j = s + 1; j < s + C_窓の種数; j++)
                {
                    if (l_ハッシュ[j] < l_ハッシュ[l_最小])
                    {
                        l_最小 = j;
                    }
                }

                if (l_最小 != l_前 && l_ハッシュ[l_最小] != ulong.MaxValue)
                {
                    l_結果.Add((l_キー[l_最小], l_最小, l_Is順向き[l_最小]));
                    l_前 = l_最小;
                }
            }

            return l_結果;
        }

        /// <summary>
        /// 錨の minimizer の索引を作る (C_種の出現の上限 を超えて出るものは捨てる)
        /// </summary>
        /// <param name="p_配列群">最終配列の ID と配列</param>
        /// <param name="p_Is錨">配列ごとに錨かどうか</param>
        /// <returns>正準形から (錨, 位置, 錨の向きの k-mer が正準形か) の一覧への表</returns>
        public static Dictionary<ulong, List<(int A_錨, int A_位置, bool A_Is順向き)>> Get_索引(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, bool[] p_Is錨)
        {
            var l_索引 = new Dictionary<ulong, List<(int A_錨, int A_位置, bool A_Is順向き)>>();
            for (var i = 0; i < p_配列群.Count; i++)
            {
                if (!p_Is錨[i])
                {
                    continue;
                }

                foreach (var (l_キー, l_位置, l_Is順向き) in Get_minimizer群(p_配列群[i].A_配列))
                {
                    if (!l_索引.TryGetValue(l_キー, out var l_一覧))
                    {
                        l_一覧 = [];
                        l_索引[l_キー] = l_一覧;
                    }

                    l_一覧.Add((i, l_位置, l_Is順向き));
                }
            }

            foreach (var l_キー in l_索引.Where(x => x.Value.Count > C_種の出現の上限).Select(x => x.Key).ToList())
            {
                _ = l_索引.Remove(l_キー);
            }

            return l_索引;
        }

        /// <summary>
        /// 1 本のリードの当たりを求める ((錨, 向き) ごとに最良の共線的な鎖を 1 本)<br/>
        /// リード上で C_重なりで捨てる割合 以上重なる当たりは、種の多いほうだけを残す (似た配列への二次の当たりを捨てる)<br/>
        /// 重なる当たりの種の数が、残すほうの C_曖昧とみなす種の比 以上なら、どちらとも決められないので両方を捨てる
        /// </summary>
        /// <param name="p_リード">リードの配列</param>
        /// <param name="p_索引">錨の minimizer の索引</param>
        /// <param name="p_長さ">配列ごとの長さ</param>
        /// <returns>当たり (リード上の始まりの順)</returns>
        public static List<ロングリードの当たり> Get_当たり群(string p_リード, Dictionary<ulong, List<(int A_錨, int A_位置, bool A_Is順向き)>> p_索引, int[] p_長さ)
        {
            var l_群 = new Dictionary<(int A_錨, bool A_Is順向き), List<(int A_r, int A_y)>>();
            foreach (var (l_キー, l_位置, l_Is順向き) in Get_minimizer群(p_リード))
            {
                if (!p_索引.TryGetValue(l_キー, out var l_一覧))
                {
                    continue;
                }

                foreach (var (l_錨, l_錨の位置, l_錨のIs順向き) in l_一覧)
                {
                    var l_Is同じ向き = l_Is順向き == l_錨のIs順向き;
                    var l_y = l_Is同じ向き ? l_錨の位置 : p_長さ[l_錨] - C_種長 - l_錨の位置;
                    if (!l_群.TryGetValue((l_錨, l_Is同じ向き), out var l_点群))
                    {
                        l_点群 = [];
                        l_群[(l_錨, l_Is同じ向き)] = l_点群;
                    }

                    l_点群.Add((l_位置, l_y));
                }
            }

            var l_結果 = new List<ロングリードの当たり>();
            foreach (var ((l_錨, l_Is同じ向き), l_点群) in l_群)
            {
                if (Get_最良の鎖(l_点群) is not { } l_鎖)
                {
                    continue;
                }

                var l_錨始 = l_Is同じ向き ? l_鎖.A_y始 : p_長さ[l_錨] - C_種長 - l_鎖.A_y終;
                var l_錨終 = (l_Is同じ向き ? l_鎖.A_y終 : p_長さ[l_錨] - C_種長 - l_鎖.A_y始) + C_種長;
                if (l_鎖.A_数 >= C_当たりに要る種の数 && l_錨終 - l_錨始 >= C_当たりに要る錨上の幅)
                {
                    l_結果.Add(new ロングリードの当たり(l_錨, l_Is同じ向き, l_鎖.A_r始, l_鎖.A_r終 + C_種長, l_錨始, l_錨終, l_鎖.A_数));
                }
            }

            var l_残す = new List<ロングリードの当たり>();
            var l_曖昧 = new HashSet<ロングリードの当たり>();
            foreach (var l_当たり in l_結果.OrderByDescending(x => x.A_種の数).ThenBy(x => x.A_錨).ThenBy(x => x.A_Is順向き))
            {
                var l_重なる = l_残す.Where(x => Is大きく重なる(x, l_当たり)).ToList();
                if (l_重なる.Count == 0)
                {
                    l_残す.Add(l_当たり);
                    continue;
                }

                foreach (var l_残した当たり in l_重なる.Where(x => l_当たり.A_種の数 >= C_曖昧とみなす種の比 * x.A_種の数))
                {
                    _ = l_曖昧.Add(l_残した当たり);
                }
            }

            return [.. l_残す.Where(x => !l_曖昧.Contains(x)).OrderBy(x => x.A_リード始).ThenBy(x => x.A_錨)];
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 確定した辺をたどって元の頂点に戻る輪を探し、環状の配列にする。輪に入った頂点の辺は外す (鎖の書き出しで二度使わないため)<br/>
        /// 輪は頂点 2i (i の小さい順) から探し、辺を順に当ててつなぎ、最後の辺の重なり (右を削る長さ) は末尾から除く (先頭の配列の始まりを保つ)
        /// </summary>
        /// <param name="p_配列群">配列の ID と配列</param>
        /// <param name="p_頂点">向き付きの配列</param>
        /// <param name="p_出る">頂点ごとの確定した出る辺 (輪の分は消す)</param>
        /// <param name="p_入る">頂点ごとの入る辺の元の頂点 (輪の分は消す)</param>
        /// <returns>環状の配列の ID・配列と、輪に入った配列の番号 (先頭が ID の元)</returns>
        private static List<(string A_ID, string A_配列, List<int> A_配列の番号群)> Get_環状の鎖群(IReadOnlyList<(string A_ID, string A_配列)> p_配列群, string[] p_頂点, (int A_先, int A_左から削る長さ, string A_埋める配列, int A_右から削る長さ)?[] p_出る, int[] p_入る)
        {
            var l_結果 = new List<(string A_ID, string A_配列, List<int> A_配列の番号群)>();
            var l_使った = new bool[p_配列群.Count];
            for (var i = 0; i < p_配列群.Count; i++)
            {
                if (l_使った[i] || p_出る[2 * i] is null)
                {
                    continue;
                }

                var l_輪 = new List<int> { 2 * i };
                var l_今 = 2 * i;
                var l_Is輪 = false;
                while (p_出る[l_今] is { } l_辺)
                {
                    if (l_辺.A_先 == 2 * i)
                    {
                        l_Is輪 = true;
                        break;
                    }

                    if (l_輪.Any(x => x / 2 == l_辺.A_先 / 2))
                    {
                        break;
                    }

                    l_輪.Add(l_辺.A_先);
                    l_今 = l_辺.A_先;
                }

                if (!l_Is輪)
                {
                    continue;
                }

                var l_鎖 = new StringBuilder(p_頂点[l_輪[0]]);
                for (var k = 0; k < l_輪.Count; k++)
                {
                    var l_辺 = p_出る[l_輪[k]]!.Value;
                    l_鎖.Length -= Math.Min(l_鎖.Length, l_辺.A_左から削る長さ);
                    _ = l_鎖.Append(l_辺.A_埋める配列);
                    if (k + 1 < l_輪.Count)
                    {
                        _ = l_鎖.Append(p_頂点[l_辺.A_先], l_辺.A_右から削る長さ, p_頂点[l_辺.A_先].Length - l_辺.A_右から削る長さ);
                    }
                    else
                    {
                        l_鎖.Length -= Math.Min(l_鎖.Length, l_辺.A_右から削る長さ);
                    }
                }

                foreach (var l_頂点 in l_輪)
                {
                    l_使った[l_頂点 / 2] = true;
                    p_出る[l_頂点] = null;
                    p_出る[l_頂点 ^ 1] = null;
                    p_入る[l_頂点] = -1;
                    p_入る[l_頂点 ^ 1] = -1;
                }

                var l_ID = p_配列群[i].A_ID.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase) ? p_配列群[i].A_ID : $"{p_配列群[i].A_ID}_{Consts.環状の目印}";
                l_結果.Add((l_ID, l_鎖.ToString(), [.. l_輪.Select(x => x / 2)]));
            }

            return l_結果;
        }

        /// <summary>
        /// 1 つの橋の間を、重なりで畳むか、グラフの経路で埋めるかを決める
        /// </summary>
        /// <param name="p_左">左の頂点の配列 (u)</param>
        /// <param name="p_右">右の頂点の配列 (x)</param>
        /// <param name="p_間">間の長さの見込み</param>
        /// <param name="p_集合">経路の材料の k-mer 集合、無ければ null</param>
        /// <param name="p_k長">経路を探す k-mer の長さ</param>
        /// <returns>埋め方</returns>
        private static ロングリードの繋ぎ方 Get_繋ぎ方(string p_左, string p_右, int p_間, LocalKmerSet? p_集合, int p_k長)
        {
            var l_幅 = (int)(C_許す幅の割合 * Math.Abs(p_間)) + C_許す幅の余裕;
            var l_重なり群 = p_間 <= l_幅 ? Get_重なり群(p_左, p_右, p_間, l_幅) : [];
            if (l_重なり群.Count == 1 && Get_最長の境界(p_右.AsSpan(0, l_重なり群[0])) < C_周期とみなす境界の長さ)
            {
                return new ロングリードの繋ぎ方(l_重なり群[0], null, null);
            }

            if (l_重なり群.Count > 1)
            {
                return new ロングリードの繋ぎ方(null, null, C_理由_一意でない);
            }

            if (p_間 + l_幅 < 0)
            {
                return new ロングリードの繋ぎ方(null, null, C_理由_重なりなし);
            }

            if (p_集合 is null)
            {
                return new ロングリードの繋ぎ方(null, null, C_理由_グラフなし);
            }

            if (p_左.Length < p_k長 || p_右.Length < p_k長 || p_左.AsSpan(p_左.Length - p_k長).IndexOfAnyExcept(Consts.塩基文字) >= 0 || p_右.AsSpan(0, p_k長).IndexOfAnyExcept(Consts.塩基文字) >= 0)
            {
                return new ロングリードの繋ぎ方(null, null, C_理由_到達不能);
            }

            var (l_経路, l_判定) = ConstrainedPathFinder.Get_経路(Util.V_変換_塩基列(p_左[^p_k長..]), Util.V_変換_塩基列(p_右[..p_k長]), Math.Max(0, p_間 - l_幅), p_間 + l_幅, p_集合, p_k長);
            return l_判定 == ギャップ充填判定.充填済み && l_経路 is not null
                ? new ロングリードの繋ぎ方(null, l_経路, null)
                : new ロングリードの繋ぎ方(null, null, l_判定 switch
                {
                    ギャップ充填判定.一意でない => C_理由_一意でない,
                    ギャップ充填判定.探索打切り => C_理由_打切り,
                    _ => C_理由_到達不能,
                });
        }

        /// <summary>
        /// グラフ (k-mer 集合) の上を u の末尾から歩き、分かれ道ではロングリードの区間の多数決で進む向きを選んで、x の先頭の k-mer まで行く道を求める<br/>
        /// 道はグラフの k-mer だけからなる。決まらなければ null
        /// </summary>
        /// <param name="p_左">左の頂点の配列 (u)</param>
        /// <param name="p_右">右の頂点の配列 (x)</param>
        /// <param name="p_間">間の長さの見込み</param>
        /// <param name="p_集合">経路の材料の k-mer 集合</param>
        /// <param name="p_k長">k-mer の長さ</param>
        /// <param name="p_区間群">橋を跨いだロングリードの区間 (u→x の向き)</param>
        /// <returns>埋める配列 (両端の k-mer を除いた間)、決まらなければ null</returns>
        private static string? Get_ロングリードで選んだ経路(string p_左, string p_右, int p_間, LocalKmerSet p_集合, int p_k長, IReadOnlyList<string> p_区間群)
        {
            if (p_左.Length < Math.Max(p_k長, C_歩き始めの整列長) || p_右.Length < p_k長)
            {
                return null;
            }

            var l_幅 = C_道選びの幅の倍率 * ((int)(C_許す幅の割合 * Math.Abs(p_間)) + C_許す幅の余裕);
            var l_位置 = Get_歩き始めの位置群(p_左, p_区間群);

            var l_目標 = Util.V_変換_塩基列(p_右[..p_k長]);
            var l_今 = Util.V_変換_塩基列(p_左[^p_k長..]);
            var l_足した = new StringBuilder();
            var l_未整列 = 0;
            var l_作業 = new byte[p_k長];
            while (l_位置.Count >= C_道選びに要る区間の数)
            {
                if (l_足した.Length >= p_k長 && l_今.AsSpan().SequenceEqual(l_目標))
                {
                    var l_埋める長さ = l_足した.Length - p_k長;
                    return Math.Abs(l_埋める長さ - p_間) <= l_幅 ? l_足した.ToString(0, l_埋める長さ) : null;
                }

                if (l_足した.Length > p_間 + l_幅 + p_k長)
                {
                    return null;
                }

                var l_候補 = new List<byte>();
                for (var l_塩基 = Consts.塩基ID.A; l_塩基 <= Consts.塩基ID.T; l_塩基++)
                {
                    Array.Copy(l_今, 1, l_作業, 0, p_k長 - 1);
                    l_作業[p_k長 - 1] = l_塩基;
                    if (p_集合.Haskmer(l_作業))
                    {
                        l_候補.Add(l_塩基);
                    }
                }

                if (l_候補.Count == 0)
                {
                    return null;
                }

                byte l_選んだ;
                if (l_候補.Count == 1)
                {
                    l_選んだ = l_候補[0];
                }
                else
                {
                    l_位置 = Get_進めた位置(l_位置, l_足した, l_未整列);
                    l_未整列 = 0;
                    if (l_位置.Count < C_道選びに要る区間の数 || Get_多数決の塩基(l_今, l_候補, p_集合, p_k長, l_位置) is not { } l_多数)
                    {
                        return null;
                    }

                    l_選んだ = l_多数;
                }

                Array.Copy(l_今, 1, l_今, 0, p_k長 - 1);
                l_今[p_k長 - 1] = l_選んだ;
                _ = l_足した.Append(Consts.塩基文字[l_選んだ - Consts.塩基ID.A]);
                l_未整列++;
                if (l_未整列 >= C_まとめて整列する長さ)
                {
                    l_位置 = Get_進めた位置(l_位置, l_足した, l_未整列);
                    l_未整列 = 0;
                }
            }

            return null;
        }

        /// <summary>
        /// 各区間で、u の末尾 C_歩き始めの整列長 bp を区間の先頭側に半大域整列し、u の末尾に当たる区間の位置を求める (誤りの多い区間は除く)
        /// </summary>
        /// <param name="p_左">左の頂点の配列 (u)</param>
        /// <param name="p_区間群">橋を跨いだロングリードの区間 (u→x の向き)</param>
        /// <returns>区間と、u の末尾に当たる位置</returns>
        private static List<(string A_区間, int A_位置)> Get_歩き始めの位置群(string p_左, IReadOnlyList<string> p_区間群)
        {
            var l_位置 = new List<(string A_区間, int A_位置)>();
            if (p_左.Length < C_歩き始めの整列長)
            {
                return l_位置;
            }

            var l_始めの配列 = p_左[^C_歩き始めの整列長..];
            foreach (var l_区間 in p_区間群)
            {
                var l_窓の長さ = Math.Min(l_区間.Length, C_区間の外側の長さ + C_端までの残りの上限 + (2 * C_区間の外側の長さ));
                var (l_距離, l_終わり) = Get_半大域の編集距離(l_始めの配列, l_区間.AsSpan(0, l_窓の長さ));
                if (l_距離 <= C_歩き始めの誤りの割合 * C_歩き始めの整列長)
                {
                    l_位置.Add((l_区間, l_終わり));
                }
            }

            return l_位置;
        }

        /// <summary>
        /// 埋めた配列を、橋を跨いだロングリードの区間に整列して位置ごとに塩基を数え、多くの区間が別の塩基を示し、
        /// その塩基を中心にした k-mer が短いリードで確かめられる所だけを置換する (挿入・欠失は直さない)
        /// </summary>
        /// <param name="p_左">左の頂点の配列 (u)</param>
        /// <param name="p_埋める配列">埋めた配列</param>
        /// <param name="p_右">右の頂点の配列 (x)</param>
        /// <param name="p_区間群">橋を跨いだロングリードの区間 (u→x の向き)</param>
        /// <param name="p_短いリードで確かめる">k-mer が短いリードで確かめられるか</param>
        /// <param name="p_直した数">直した塩基の数</param>
        /// <returns>直した埋める配列</returns>
        private static string Get_直した埋める配列(string p_左, string p_埋める配列, string p_右, IReadOnlyList<string> p_区間群, Func<string, bool> p_短いリードで確かめる, out int p_直した数)
        {
            p_直した数 = 0;
            var l_位置群 = Get_歩き始めの位置群(p_左, p_区間群);
            if (p_埋める配列.Length == 0 || l_位置群.Count < C_直すのに要る区間の数)
            {
                return p_埋める配列;
            }

            var l_整列する配列 = p_埋める配列 + p_右[..Math.Min(C_直すときに足す右の長さ, p_右.Length)];
            var l_数 = new int[p_埋める配列.Length, Consts.塩基文字.Length + 1];
            var l_区間の数 = new int[p_埋める配列.Length];
            foreach (var (l_区間, l_始め) in l_位置群)
            {
                var l_位置 = l_始め;
                for (var p = 0; p < l_整列する配列.Length && p < p_埋める配列.Length; p += C_まとめて整列する長さ)
                {
                    var l_長さ = Math.Min(C_まとめて整列する長さ, l_整列する配列.Length - p);
                    var l_窓の長さ = Math.Min(l_区間.Length - l_位置, (int)(C_窓の倍率 * l_長さ) + C_窓の余裕);
                    if (l_窓の長さ <= 0)
                    {
                        break;
                    }

                    var (l_距離, l_終わり, l_対応) = Get_塩基の対応(l_整列する配列.AsSpan(p, l_長さ), l_区間.AsSpan(l_位置, l_窓の長さ));
                    if (l_距離 > C_歩きの誤りの割合 * l_長さ)
                    {
                        break;
                    }

                    for (var i = 0; i < l_長さ && p + i < p_埋める配列.Length; i++)
                    {
                        var l_番号 = l_対応[i] == C_欠失の文字 ? Consts.塩基文字.Length : Consts.塩基文字.IndexOf(l_対応[i]);
                        if (l_番号 >= 0)
                        {
                            l_数[p + i, l_番号]++;
                            l_区間の数[p + i]++;
                        }
                    }

                    l_位置 += l_終わり;
                }
            }

            var l_片側 = C_確かめるkmerの片側の長さ;
            var l_文脈 = (p_左[^Math.Min(l_片側, p_左.Length)..] + p_埋める配列 + p_右[..Math.Min(l_片側, p_右.Length)]).ToCharArray();
            var l_ずれ = Math.Min(l_片側, p_左.Length);
            for (var i = 0; i < p_埋める配列.Length; i++)
            {
                var l_n = l_区間の数[i];
                var l_元 = Consts.塩基文字.IndexOf(p_埋める配列[i]);
                if (l_n < C_直すのに要る区間の数 || l_元 < 0 || l_数[i, l_元] > C_元の塩基の割合の上限 * l_n)
                {
                    continue;
                }

                for (var c = 0; c < Consts.塩基文字.Length; c++)
                {
                    var l_中心 = l_ずれ + i;
                    if (c == l_元 || l_数[i, c] < C_直す塩基の割合 * l_n || l_中心 < l_片側 || l_中心 + l_片側 >= l_文脈.Length)
                    {
                        continue;
                    }

                    var l_元の文字 = l_文脈[l_中心];
                    l_文脈[l_中心] = Consts.塩基文字[c];
                    if (p_短いリードで確かめる(new string(l_文脈, l_中心 - l_片側, (2 * l_片側) + 1)))
                    {
                        p_直した数++;
                    }
                    else
                    {
                        l_文脈[l_中心] = l_元の文字;
                    }

                    break;
                }
            }

            return new string(l_文脈, l_ずれ, p_埋める配列.Length);
        }

        /// <summary>
        /// 配列を窓に半大域整列し (配列は全体、窓の両端は無料)、配列の位置ごとに対応する窓の塩基 (欠失は C_欠失の文字) を求める
        /// </summary>
        /// <param name="p_配列">整列する配列</param>
        /// <param name="p_窓">窓</param>
        /// <returns>編集距離、窓での終わりの位置 (含まない)、配列の位置ごとの対応する塩基</returns>
        private static (int A_距離, int A_終わり, char[] A_対応) Get_塩基の対応(ReadOnlySpan<char> p_配列, ReadOnlySpan<char> p_窓)
        {
            var l_m = p_配列.Length;
            var l_n = p_窓.Length;
            var l_表 = new int[l_m + 1, l_n + 1];
            for (var i = 1; i <= l_m; i++)
            {
                l_表[i, 0] = i;
                for (var j = 1; j <= l_n; j++)
                {
                    var l_一致 = l_表[i - 1, j - 1] + (p_配列[i - 1] == p_窓[j - 1] ? 0 : 1);
                    l_表[i, j] = Math.Min(l_一致, Math.Min(l_表[i - 1, j], l_表[i, j - 1]) + 1);
                }
            }

            var l_終わり = 0;
            for (var j = 1; j <= l_n; j++)
            {
                if (l_表[l_m, j] < l_表[l_m, l_終わり])
                {
                    l_終わり = j;
                }
            }

            var l_対応 = new char[l_m];
            var (l_i, l_j) = (l_m, l_終わり);
            while (l_i > 0)
            {
                if (l_j > 0 && l_表[l_i, l_j] == l_表[l_i - 1, l_j - 1] + (p_配列[l_i - 1] == p_窓[l_j - 1] ? 0 : 1))
                {
                    l_対応[l_i - 1] = p_窓[l_j - 1];
                    l_i--;
                    l_j--;
                }
                else if (l_表[l_i, l_j] == l_表[l_i - 1, l_j] + 1)
                {
                    l_対応[l_i - 1] = C_欠失の文字;
                    l_i--;
                }
                else
                {
                    l_j--;
                }
            }

            return (l_表[l_m, l_終わり], l_終わり, l_対応);
        }

        /// <summary>
        /// 足した塩基のうち未整列の分を、各区間の現在位置から始まる窓に整列して位置を進める。誤りの多い区間は捨てる
        /// </summary>
        /// <param name="p_位置">区間と現在位置</param>
        /// <param name="p_足した">足した塩基</param>
        /// <param name="p_未整列">足した塩基のうち未整列の長さ</param>
        /// <returns>進めた区間と位置</returns>
        private static List<(string A_区間, int A_位置)> Get_進めた位置(List<(string A_区間, int A_位置)> p_位置, StringBuilder p_足した, int p_未整列)
        {
            if (p_未整列 == 0)
            {
                return p_位置;
            }

            var l_配列 = p_足した.ToString(p_足した.Length - p_未整列, p_未整列);
            var l_結果 = new List<(string A_区間, int A_位置)>();
            foreach (var (l_区間, l_位置) in p_位置)
            {
                var l_窓の長さ = Math.Min(l_区間.Length - l_位置, (int)(C_窓の倍率 * p_未整列) + C_窓の余裕);
                if (l_窓の長さ <= 0)
                {
                    continue;
                }

                var (l_距離, l_終わり) = Get_半大域の編集距離(l_配列, l_区間.AsSpan(l_位置, l_窓の長さ));
                if (l_距離 <= C_歩きの誤りの割合 * p_未整列)
                {
                    l_結果.Add((l_区間, l_位置 + l_終わり));
                }
            }

            return l_結果;
        }

        /// <summary>
        /// 分かれ道で、候補の塩基ごとの先読みの配列を各区間の現在位置からの窓に整列し、編集距離が最小の候補に票を入れ、多数決の塩基を求める
        /// </summary>
        /// <param name="p_今">現在の k-mer</param>
        /// <param name="p_候補">候補の塩基</param>
        /// <param name="p_集合">k-mer 集合</param>
        /// <param name="p_k長">k-mer の長さ</param>
        /// <param name="p_位置">区間と現在位置</param>
        /// <returns>首位の票が C_要る票 以上かつ次点の C_道選びの優勢比 倍以上の塩基、無ければ null</returns>
        private static byte? Get_多数決の塩基(byte[] p_今, List<byte> p_候補, LocalKmerSet p_集合, int p_k長, List<(string A_区間, int A_位置)> p_位置)
        {
            var l_先読み = p_候補.Select(x => Get_先読み群(p_今, x, p_集合, p_k長)).ToList();
            var l_票 = new int[p_候補.Count];
            var l_窓の長さ = (int)(C_窓の倍率 * C_先読み長) + C_窓の余裕;
            foreach (var (l_区間, l_位置) in p_位置)
            {
                var l_窓 = l_区間.AsSpan(l_位置, Math.Min(l_窓の長さ, l_区間.Length - l_位置));
                var l_最小 = int.MaxValue;
                var l_最小の候補 = -1;
                var l_並ぶ = false;
                for (var c = 0; c < p_候補.Count; c++)
                {
                    var l_距離 = int.MaxValue;
                    foreach (var l_先 in l_先読み[c])
                    {
                        l_距離 = Math.Min(l_距離, Get_半大域の編集距離(l_先, l_窓).A_距離);
                    }

                    if (l_距離 < l_最小)
                    {
                        l_最小 = l_距離;
                        l_最小の候補 = c;
                        l_並ぶ = false;
                    }
                    else if (l_距離 == l_最小)
                    {
                        l_並ぶ = true;
                    }
                }

                if (!l_並ぶ && l_最小の候補 >= 0)
                {
                    l_票[l_最小の候補]++;
                }
            }

            var l_順 = Enumerable.Range(0, p_候補.Count).OrderByDescending(x => l_票[x]).ToList();
            var l_首位 = l_票[l_順[0]];
            var l_次点 = l_順.Count > 1 ? l_票[l_順[1]] : 0;
            return l_首位 >= C_要る票 && l_首位 >= C_道選びの優勢比 * l_次点 ? p_候補[l_順[0]] : null;
        }

        /// <summary>
        /// 現在の k-mer に候補の塩基を足したところから、グラフの上で到達できる長さ C_先読み長 の配列を集める (C_先読みの数の上限 本まで)
        /// </summary>
        /// <param name="p_今">現在の k-mer</param>
        /// <param name="p_塩基">候補の塩基</param>
        /// <param name="p_集合">k-mer 集合</param>
        /// <param name="p_k長">k-mer の長さ</param>
        /// <returns>先読みの配列</returns>
        private static List<string> Get_先読み群(byte[] p_今, byte p_塩基, LocalKmerSet p_集合, int p_k長)
        {
            var l_結果 = new List<string>();
            var l_始め = new byte[p_k長];
            Array.Copy(p_今, 1, l_始め, 0, p_k長 - 1);
            l_始め[p_k長 - 1] = p_塩基;
            var l_積み = new Stack<(byte[] A_kmer, string A_配列)>();
            l_積み.Push((l_始め, Consts.塩基文字[p_塩基 - Consts.塩基ID.A].ToString()));
            while (l_積み.Count > 0 && l_結果.Count < C_先読みの数の上限)
            {
                var (l_kmer, l_配列) = l_積み.Pop();
                if (l_配列.Length >= C_先読み長)
                {
                    l_結果.Add(l_配列);
                    continue;
                }

                var l_伸びた = false;
                for (var l_次 = Consts.塩基ID.T; l_次 >= Consts.塩基ID.A; l_次--)
                {
                    var l_新 = new byte[p_k長];
                    Array.Copy(l_kmer, 1, l_新, 0, p_k長 - 1);
                    l_新[p_k長 - 1] = l_次;
                    if (p_集合.Haskmer(l_新))
                    {
                        l_積み.Push((l_新, l_配列 + Consts.塩基文字[l_次 - Consts.塩基ID.A]));
                        l_伸びた = true;
                    }
                }

                if (!l_伸びた)
                {
                    l_結果.Add(l_配列);
                }
            }

            return l_結果;
        }

        /// <summary>
        /// 配列を窓に半大域整列した編集距離 (配列は全体、窓の両端は無料) と、整列の窓での終わりの位置を求める
        /// </summary>
        /// <param name="p_配列">整列する配列</param>
        /// <param name="p_窓">窓</param>
        /// <returns>編集距離と、窓での終わりの位置 (含まない)</returns>
        private static (int A_距離, int A_終わり) Get_半大域の編集距離(ReadOnlySpan<char> p_配列, ReadOnlySpan<char> p_窓)
        {
            var l_前 = new int[p_窓.Length + 1];
            var l_今 = new int[p_窓.Length + 1];
            for (var i = 1; i <= p_配列.Length; i++)
            {
                l_今[0] = i;
                for (var j = 1; j <= p_窓.Length; j++)
                {
                    var l_一致 = l_前[j - 1] + (p_配列[i - 1] == p_窓[j - 1] ? 0 : 1);
                    l_今[j] = Math.Min(l_一致, Math.Min(l_前[j], l_今[j - 1]) + 1);
                }

                (l_前, l_今) = (l_今, l_前);
            }

            var l_最小 = p_配列.Length;
            var l_終わり = 0;
            for (var j = 0; j <= p_窓.Length; j++)
            {
                if (l_前[j] < l_最小)
                {
                    l_最小 = l_前[j];
                    l_終わり = j;
                }
            }

            return (l_最小, l_終わり);
        }

        /// <summary>
        /// 左の末尾と右の先頭が完全に一致する重なりの長さを、許す幅 (−間 ± 幅) の中ですべて求める<br/>
        /// 2 通り以上あれば (縦列反復の中など) 重なりは一意に決まらない
        /// </summary>
        /// <param name="p_左">左の頂点の配列</param>
        /// <param name="p_右">右の頂点の配列</param>
        /// <param name="p_間">間の長さの見込み (負は重なり)</param>
        /// <param name="p_幅">許す幅</param>
        /// <returns>重なりの長さ (短い順)</returns>
        private static List<int> Get_重なり群(string p_左, string p_右, int p_間, int p_幅)
        {
            var l_下限 = Math.Max(C_重なりの最短長, -p_間 - p_幅);
            var l_上限 = Math.Min(Math.Min(p_左.Length, p_右.Length) - 1, -p_間 + p_幅);
            var l_結果 = new List<int>();
            for (var o = l_下限; o <= l_上限; o++)
            {
                if (p_左.AsSpan(p_左.Length - o).SequenceEqual(p_右.AsSpan(0, o)))
                {
                    l_結果.Add(o);
                }
            }

            return l_結果;
        }

        /// <summary>
        /// 配列の先頭と末尾が一致する最長の長さ (配列全体を除く境界、KMP の失敗関数の最後の値) を求める
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <returns>最長の境界の長さ</returns>
        private static int Get_最長の境界(ReadOnlySpan<char> p_配列)
        {
            var l_失敗 = new int[p_配列.Length];
            for (var i = 1; i < p_配列.Length; i++)
            {
                var l_j = l_失敗[i - 1];
                while (l_j > 0 && p_配列[i] != p_配列[l_j])
                {
                    l_j = l_失敗[l_j - 1];
                }

                l_失敗[i] = p_配列[i] == p_配列[l_j] ? l_j + 1 : 0;
            }

            return p_配列.Length == 0 ? 0 : l_失敗[^1];
        }

        /// <summary>
        /// 経路の材料の配列と、橋の両端の配列の末尾・先頭から、k-mer 集合を作る
        /// </summary>
        /// <param name="p_材料">経路の材料の配列</param>
        /// <param name="p_橋群">決めた橋</param>
        /// <param name="p_頂点">向き付きの配列</param>
        /// <param name="p_k長">k-mer の長さ</param>
        /// <returns>k-mer 集合</returns>
        private static LocalKmerSet Get_経路の材料(IEnumerable<string> p_材料, IReadOnlyList<ロングリードの橋> p_橋群, string[] p_頂点, int p_k長)
        {
            var l_集合 = new LocalKmerSet(p_k長);
            foreach (var l_配列 in p_材料)
            {
                V_登録_kmer(l_集合, l_配列, p_k長);
            }

            var l_端の長さ = p_k長 + C_材料にする端の長さ;
            foreach (var l_橋 in p_橋群)
            {
                var l_左 = p_頂点[l_橋.A_u];
                var l_右 = p_頂点[l_橋.A_x];
                V_登録_kmer(l_集合, l_左[Math.Max(0, l_左.Length - l_端の長さ)..], p_k長);
                V_登録_kmer(l_集合, l_右[..Math.Min(l_右.Length, l_端の長さ)], p_k長);
            }

            return l_集合;
        }

        /// <summary>
        /// 配列の A/C/G/T だけからなる k-mer をすべて集合に入れる (集合は正準形で持つので片向きでよい)
        /// </summary>
        /// <param name="p_集合">k-mer 集合</param>
        /// <param name="p_配列">配列</param>
        /// <param name="p_k長">k-mer の長さ</param>
        private static void V_登録_kmer(LocalKmerSet p_集合, string p_配列, int p_k長)
        {
            if (p_配列.Length < p_k長)
            {
                return;
            }

            var l_塩基 = Util.V_変換_塩基列(p_配列);
            var l_連続 = 0;
            for (var i = 0; i < l_塩基.Length; i++)
            {
                l_連続 = l_塩基[i] is >= Consts.塩基ID.A and <= Consts.塩基ID.T ? l_連続 + 1 : 0;
                if (l_連続 >= p_k長)
                {
                    p_集合.V_登録(l_塩基.AsSpan(i - p_k長 + 1, p_k長));
                }
            }
        }

        /// <summary>
        /// 埋め方の TSV の印
        /// </summary>
        /// <param name="p_繋ぎ方">埋め方</param>
        /// <returns>印</returns>
        private static string Get_繋ぎ方の印(ロングリードの繋ぎ方 p_繋ぎ方)
        {
            return p_繋ぎ方.A_重なり is { } l_重なり ? C_重なりの印 + l_重なり : p_繋ぎ方.A_埋める配列 is { } l_配列 ? (p_繋ぎ方.A_Isロングリードで選んだ ? C_道選びの印 : C_経路の印) + l_配列.Length : C_繋がない印 + p_繋ぎ方.A_理由;
        }

        /// <summary>
        /// 当たりの点の並びから、点の数が最大の共線的な鎖を求める
        /// </summary>
        /// <param name="p_点群">(リード上の位置, 錨上の位置を向きに合わせて増える向きにしたもの)</param>
        /// <returns>鎖の点の数と、両端の点の位置。点が無ければ null</returns>
        private static (int A_数, int A_r始, int A_r終, int A_y始, int A_y終)? Get_最良の鎖(List<(int A_r, int A_y)> p_点群)
        {
            if (p_点群.Count == 0)
            {
                return null;
            }

            p_点群.Sort();
            var l_数 = new int[p_点群.Count];
            var l_前 = new int[p_点群.Count];
            var l_最良 = 0;
            for (var i = 0; i < p_点群.Count; i++)
            {
                l_数[i] = 1;
                l_前[i] = -1;
                for (var j = i - 1; j >= Math.Max(0, i - C_鎖を遡る数); j--)
                {
                    var l_dr = p_点群[i].A_r - p_点群[j].A_r;
                    var l_dy = p_点群[i].A_y - p_点群[j].A_y;
                    if (l_dr < 1 || l_dr > C_鎖の差の上限 || l_dy < 1 || l_dy > C_鎖の差の上限)
                    {
                        continue;
                    }

                    if (Math.Abs(l_dr - l_dy) > (C_鎖のずれの割合 * Math.Max(l_dr, l_dy)) + C_鎖のずれの余裕)
                    {
                        continue;
                    }

                    if (l_数[j] + 1 > l_数[i])
                    {
                        l_数[i] = l_数[j] + 1;
                        l_前[i] = j;
                    }
                }

                if (l_数[i] > l_数[l_最良])
                {
                    l_最良 = i;
                }
            }

            var l_始め = l_最良;
            while (l_前[l_始め] >= 0)
            {
                l_始め = l_前[l_始め];
            }

            return (l_数[l_最良], p_点群[l_始め].A_r, p_点群[l_最良].A_r, p_点群[l_始め].A_y, p_点群[l_最良].A_y);
        }

        /// <summary>
        /// 2 つの当たりが、リード上で短いほうの長さの C_重なりで捨てる割合 以上重なるか
        /// </summary>
        /// <param name="p_a">当たり</param>
        /// <param name="p_b">当たり</param>
        /// <returns>大きく重なるなら true</returns>
        private static bool Is大きく重なる(ロングリードの当たり p_a, ロングリードの当たり p_b)
        {
            var l_重なり = Math.Min(p_a.A_リード終, p_b.A_リード終) - Math.Max(p_a.A_リード始, p_b.A_リード始);
            var l_短いほう = Math.Min(p_a.A_リード終 - p_a.A_リード始, p_b.A_リード終 - p_b.A_リード始);
            return l_重なり >= C_重なりで捨てる割合 * l_短いほう;
        }

        /// <summary>
        /// 錨ごとの区画 (C_区画長) の深さ (区画を丸ごと覆う当たりの数) を求める
        /// </summary>
        /// <param name="p_当たり群">リードごとの当たり</param>
        /// <param name="p_長さ">配列ごとの長さ</param>
        /// <param name="p_Is錨">配列ごとに錨かどうか</param>
        /// <returns>配列ごとの区画の深さ (錨でなければ空)</returns>
        private static int[][] Get_区画の深さ(List<List<ロングリードの当たり>> p_当たり群, int[] p_長さ, bool[] p_Is錨)
        {
            var l_深さ = p_長さ.Select((x, i) => p_Is錨[i] ? new int[x / C_区画長] : []).ToArray();
            foreach (var l_当たり in p_当たり群.SelectMany(x => x))
            {
                var l_区画 = l_深さ[l_当たり.A_錨];
                var l_始め = (l_当たり.A_錨始 + C_区画長 - 1) / C_区画長;
                var l_終わり = Math.Min(l_区画.Length, l_当たり.A_錨終 / C_区画長);
                for (var b = l_始め; b < l_終わり; b++)
                {
                    l_区画[b]++;
                }
            }

            return l_深さ;
        }

        /// <summary>
        /// 全錨の、深さが 0 でない区画の深さの中央値を求める
        /// </summary>
        /// <param name="p_区画の深さ">配列ごとの区画の深さ</param>
        /// <returns>基準の深さ、区画が無ければ 0</returns>
        private static int Get_基準の深さ(int[][] p_区画の深さ)
        {
            var l_値 = p_区画の深さ.SelectMany(x => x).Where(x => x > 0).ToList();
            return l_値.Count == 0 ? 0 : StatsUtil.Get_中央値(l_値);
        }

        /// <summary>
        /// 錨の深さ (深さが 0 でない区画の深さの中央値) を求める
        /// </summary>
        /// <param name="p_区画">錨の区画の深さ</param>
        /// <returns>錨の深さ、区画が無ければ 0</returns>
        private static int Get_錨の深さ(int[] p_区画)
        {
            var l_値 = p_区画.Where(x => x > 0).ToList();
            return l_値.Count == 0 ? 0 : StatsUtil.Get_中央値(l_値);
        }

        /// <summary>
        /// 当たりが、丸ごと覆う反復でない区画を C_証拠に要る一意な区画の数 以上持つか
        /// </summary>
        /// <param name="p_当たり">当たり</param>
        /// <param name="p_区画">当たりの錨の区画の深さ</param>
        /// <param name="p_基準">反復の区画を決める基準の深さ (全体の基準と錨の深さの大きいほう)</param>
        /// <returns>橋の証拠に使えるなら true</returns>
        private static bool Is証拠に使える(ロングリードの当たり p_当たり, int[] p_区画, int p_基準)
        {
            var l_始め = (p_当たり.A_錨始 + C_区画長 - 1) / C_区画長;
            var l_終わり = Math.Min(p_区画.Length, p_当たり.A_錨終 / C_区画長);
            var l_数 = 0;
            for (var b = l_始め; b < l_終わり; b++)
            {
                if (p_区画[b] < C_多コピーとみなす深さの比 * p_基準)
                {
                    l_数++;
                }
            }

            return l_数 >= C_証拠に要る一意な区画の数;
        }

        /// <summary>
        /// 1 本のリードの証拠に使える当たりから、隣り合う当たりの組を橋の票として数える (同じ橋には 1 票だけ)
        /// </summary>
        /// <param name="p_当たり">証拠に使える当たり (リード上の始まりの順)</param>
        /// <param name="p_リード長">リードの長さ</param>
        /// <param name="p_長さ">配列ごとの長さ</param>
        /// <param name="p_票">向き付きの辺 (u→x と x^1→u^1 の両方) ごとの間の長さの一覧</param>
        private static void V_数える_証拠(List<ロングリードの当たり> p_当たり, int p_リード長, int[] p_長さ, Dictionary<(int A_u, int A_x), List<int>> p_票)
        {
            foreach (var (l_u, l_x, l_間, _, _) in Get_証拠群(p_当たり, p_リード長, p_長さ))
            {
                V_追加_票(p_票, (l_u, l_x), l_間);
                V_追加_票(p_票, (l_x ^ 1, l_u ^ 1), l_間);
            }
        }

        /// <summary>
        /// 1 本のリードの証拠に使える当たりから、隣り合う当たりの組 (橋の証拠) を求める<br/>
        /// 当たりの端から錨の端までが C_端までの残りの上限 を超える組は使わない<br/>
        /// 前の当たりは、リードの手前に C_端までの残りの上限 を超える配列があるのに錨の途中から入ったものなら使わない (後の当たりも向きを逆にして同じ)<br/>
        /// 別の場所から錨の端の反復の部分にだけ入ったリードを除くため。同じ橋は 1 回だけ
        /// </summary>
        /// <param name="p_当たり">証拠に使える当たり (リード上の始まりの順)</param>
        /// <param name="p_リード長">リードの長さ</param>
        /// <param name="p_長さ">配列ごとの長さ</param>
        /// <returns>(u, x, 間の長さ, 前の当たり, 後の当たり) の一覧</returns>
        private static List<(int A_u, int A_x, int A_間, ロングリードの当たり A_前, ロングリードの当たり A_後)> Get_証拠群(List<ロングリードの当たり> p_当たり, int p_リード長, int[] p_長さ)
        {
            var l_結果 = new List<(int A_u, int A_x, int A_間, ロングリードの当たり A_前, ロングリードの当たり A_後)>();
            var l_数えた = new HashSet<(int A_u, int A_x)>();
            for (var i = 1; i < p_当たり.Count; i++)
            {
                var l_前 = p_当たり[i - 1];
                var l_後 = p_当たり[i];
                if (l_前.A_錨 == l_後.A_錨)
                {
                    continue;
                }

                var l_u = (2 * l_前.A_錨) + (l_前.A_Is順向き ? 0 : 1);
                var l_x = (2 * l_後.A_錨) + (l_後.A_Is順向き ? 0 : 1);
                var l_前の残り = l_前.A_Is順向き ? p_長さ[l_前.A_錨] - l_前.A_錨終 : l_前.A_錨始;
                var l_後の手前 = l_後.A_Is順向き ? l_後.A_錨始 : p_長さ[l_後.A_錨] - l_後.A_錨終;
                if (l_前の残り > C_端までの残りの上限 || l_後の手前 > C_端までの残りの上限)
                {
                    continue;
                }

                var l_前の入り = l_前.A_Is順向き ? l_前.A_錨始 : p_長さ[l_前.A_錨] - l_前.A_錨終;
                var l_後の出 = l_後.A_Is順向き ? p_長さ[l_後.A_錨] - l_後.A_錨終 : l_後.A_錨始;
                if ((l_前の入り > C_端までの残りの上限 && l_前.A_リード始 > C_端までの残りの上限) || (l_後の出 > C_端までの残りの上限 && p_リード長 - l_後.A_リード終 > C_端までの残りの上限))
                {
                    continue;
                }

                if (l_数えた.Add(Get_正準辺(l_u, l_x)))
                {
                    l_結果.Add((l_u, l_x, l_後.A_リード始 - l_前.A_リード終 - l_前の残り - l_後の手前, l_前, l_後));
                }
            }

            return l_結果;
        }

        /// <summary>
        /// ロングリードをもう一度読み、決めた橋ごとに、それを跨いだリードの区間を取り出す (u→x の向き、橋ごとに C_区間の数の上限 本まで)<br/>
        /// 区間は、前の当たりの終わりの C_区間の外側の長さ 手前から、後の当たりの始まりの C_区間の外側の長さ 先まで
        /// </summary>
        /// <param name="p_ロングリード群">ロングリードの塩基列</param>
        /// <param name="p_索引">錨の minimizer の索引</param>
        /// <param name="p_長さ">配列ごとの長さ</param>
        /// <param name="p_Is使える">当たりが証拠に使えるか</param>
        /// <param name="p_橋群">決めた橋</param>
        /// <param name="p_並列数">並列数</param>
        /// <returns>橋ごとの区間 (橋と同じ順)</returns>
        private static List<List<string>> Get_区間群(IEnumerable<string> p_ロングリード群, Dictionary<ulong, List<(int A_錨, int A_位置, bool A_Is順向き)>> p_索引, int[] p_長さ, Func<ロングリードの当たり, bool> p_Is使える, List<ロングリードの橋> p_橋群, int p_並列数)
        {
            var l_番号 = new Dictionary<(int A_u, int A_x), int>();
            for (var b = 0; b < p_橋群.Count; b++)
            {
                l_番号[(p_橋群[b].A_u, p_橋群[b].A_x)] = b;
            }

            var l_区間群 = p_橋群.Select(_ => new List<string>()).ToList();
            if (p_橋群.Count == 0)
            {
                return l_区間群;
            }

            var l_取り出し = p_ロングリード群.Where(x => x.Length >= C_リードの最短長).AsParallel().AsOrdered().WithDegreeOfParallelism(p_並列数).Select(l_リード =>
            {
                var l_結果 = new List<(int A_橋, string A_区間)>();
                foreach (var (l_u, l_x, _, l_前, l_後) in Get_証拠群(Get_当たり群(l_リード, p_索引, p_長さ).Where(p_Is使える).ToList(), l_リード.Length, p_長さ))
                {
                    var l_正準 = Get_正準辺(l_u, l_x);
                    if (!l_番号.TryGetValue(l_正準, out var l_橋))
                    {
                        continue;
                    }

                    var l_始め = Math.Max(0, l_前.A_リード終 - C_区間の外側の長さ);
                    var l_終わり = Math.Min(l_リード.Length, l_後.A_リード始 + C_区間の外側の長さ);
                    if (l_終わり <= l_始め)
                    {
                        continue;
                    }

                    var l_区間 = l_リード[l_始め..l_終わり];
                    l_結果.Add((l_橋, l_正準 == (l_u, l_x) ? l_区間 : Util.V_逆相補_曖昧塩基あり(l_区間)));
                }

                return l_結果;
            });
            foreach (var (l_橋, l_区間) in l_取り出し.SelectMany(x => x))
            {
                if (l_区間群[l_橋].Count < C_区間の数の上限)
                {
                    l_区間群[l_橋].Add(l_区間);
                }
            }

            return l_区間群;
        }

        /// <summary>
        /// 向き付きの辺の票に間の長さを 1 つ足す
        /// </summary>
        /// <param name="p_票">向き付きの辺ごとの間の長さの一覧</param>
        /// <param name="p_辺">向き付きの辺</param>
        /// <param name="p_間">間の長さ</param>
        private static void V_追加_票(Dictionary<(int A_u, int A_x), List<int>> p_票, (int A_u, int A_x) p_辺, int p_間)
        {
            if (!p_票.TryGetValue(p_辺, out var l_一覧))
            {
                l_一覧 = [];
                p_票[p_辺] = l_一覧;
            }

            l_一覧.Add(p_間);
        }

        /// <summary>
        /// 頂点ごとの首位の相手を決め、互いに首位どうしの組を橋とする
        /// </summary>
        /// <param name="p_票">向き付きの辺ごとの間の長さの一覧</param>
        /// <returns>決めた橋 (正準の向きだけ)</returns>
        private static List<ロングリードの橋> Get_互いに一意な橋(Dictionary<(int A_u, int A_x), List<int>> p_票)
        {
            var l_首位 = new Dictionary<int, (int A_x, int A_票, int A_次点の票)>();
            foreach (var l_群 in p_票.GroupBy(x => x.Key.A_u))
            {
                var l_並び = l_群.OrderByDescending(x => x.Value.Count).ThenBy(x => x.Key.A_x).ToList();
                var l_首位の票 = l_並び[0].Value.Count;
                var l_次点の票 = l_並び.Count > 1 ? l_並び[1].Value.Count : 0;
                if (l_首位の票 >= C_要る票 && l_首位の票 >= C_首位の優勢比 * l_次点の票)
                {
                    l_首位[l_群.Key] = (l_並び[0].Key.A_x, l_首位の票, l_次点の票);
                }
            }

            var l_橋群 = new List<ロングリードの橋>();
            foreach (var (l_u, (l_x, l_票, l_次点の票)) in l_首位.OrderBy(x => x.Key))
            {
                if (Get_正準辺(l_u, l_x) != (l_u, l_x) || !l_首位.TryGetValue(l_x ^ 1, out var l_逆) || l_逆.A_x != (l_u ^ 1))
                {
                    continue;
                }

                l_橋群.Add(new ロングリードの橋(l_u, l_x, l_票, l_次点の票, StatsUtil.Get_中央値(p_票[(l_u, l_x)])));
            }

            return l_橋群;
        }

        /// <summary>
        /// 向き付きの辺 (u→x) を正準化する。u→x と x^1→u^1 は同じ繋ぎ
        /// </summary>
        /// <param name="p_u">左の頂点</param>
        /// <param name="p_x">右の頂点</param>
        /// <returns>正準化した (左, 右)</returns>
        private static (int A_u, int A_x) Get_正準辺(int p_u, int p_x)
        {
            return (p_u, p_x).CompareTo((p_x ^ 1, p_u ^ 1)) <= 0 ? (p_u, p_x) : (p_x ^ 1, p_u ^ 1);
        }

        /// <summary>
        /// 頂点の向きの印 (偶数は順向き、奇数は逆向き)
        /// </summary>
        /// <param name="p_頂点">頂点</param>
        /// <returns>向きの印</returns>
        private static string Get_向きの印(int p_頂点)
        {
            return p_頂点 % 2 == 0 ? C_順向きの印 : C_逆向きの印;
        }

        /// <summary>
        /// k-mer の正準形を混ぜた値 (minimizer の順位に使う)
        /// </summary>
        /// <param name="p_値">k-mer の正準形</param>
        /// <returns>混ぜた値</returns>
        private static ulong Get_混ぜた値(ulong p_値)
        {
            var l_値 = p_値 + 0x9E3779B97F4A7C15UL;
            l_値 = (l_値 ^ (l_値 >> 30)) * 0xBF58476D1CE4E5B9UL;
            l_値 = (l_値 ^ (l_値 >> 27)) * 0x94D049BB133111EBUL;
            return l_値 ^ (l_値 >> 31);
        }

        #endregion
    }
}
