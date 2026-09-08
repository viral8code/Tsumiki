using System.Text;
using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// 複数の k のアセンブリを1つに統合する。
    ///
    /// k を変えると壊れ方が逆になる。低い k は反復配列を潰して短く切れ、
    /// 高い k はカバレッジが薄くなってグラフが千切れる。したがって片方が
    /// 途切れている場所を、もう片方が繋いでいることがある。
    ///
    /// 最も良い1本を骨格に選び、他の k の配列がその骨格の2本を跨いでいる箇所を
    /// 隣接の証拠として集める。採否は相互一意性で判定し、跨いだ配列を
    /// そのまま繋ぎ目に使う。
    /// </summary>
    internal static class AssemblyMerger
    {
        /// <summary>
        /// 骨格側の末端でアンカーを探す長さ。ここより内側でしか一致しない配列は、
        /// 末端を跨いでいる証拠にならない。
        /// </summary>
        private const int 末端とみなす長さ = 2000;

        /// <summary>
        /// 橋渡しとして認める最大の挟み込み長。これを超える隙間を1本の配列で
        /// 埋めるのは、骨格側が丸ごと取りこぼした領域を持ち込むことになり、
        /// 骨格を選んだ判断と矛盾する。
        /// </summary>
        private const int 橋渡し長の上限 = 50_000;

        /// <summary>
        /// 連結を認めるために必要な、独立に同じ隣接を主張した k の数。
        ///
        /// 骨格が途切れているのは、そこで繋ぐ根拠が足りないと判断した結果である
        /// ことが多い。1つの k の1本の配列だけでその判断を覆すと、その配列自身が
        /// 誤アセンブリだった場合にそのまま持ち込むことになる。
        /// </summary>
        private const int 必要な独立支持数の既定値 = 2;

        /// <summary>
        /// 骨格に対して他の候補を統合し、結果を p_出力パス へ書き出す。
        /// 繋げた箇所が1つも無ければ false を返す(その場合、出力は行わない)。
        /// </summary>
        public static bool V_統合(
            アセンブリ実行結果 p_骨格,
            IReadOnlyList<アセンブリ実行結果> p_全候補,
            int p_アンカーk長,
            string p_出力パス,
            int p_必要な独立支持数 = 必要な独立支持数の既定値)
        {
            var (l_骨格名一覧, l_骨格配列) = Get_配列一覧(p_骨格.A_最終パス);
            if (l_骨格配列.Count == 0)
            {
                return false;
            }

            var l_索引 = Get_骨格索引(l_骨格配列, p_アンカーk長);
            var l_候補 = new List<橋渡し候補>();
            foreach (var l_他 in p_全候補)
            {
                if (l_他.A_k長 == p_骨格.A_k長)
                {
                    continue;
                }
                l_候補.AddRange(Get_橋渡し候補(l_他, l_索引, l_骨格配列, p_アンカーk長));
            }

            var l_確定 = Get_相互一意な橋渡し(l_候補, l_骨格配列.Count, p_必要な独立支持数);
            if (l_確定.Count == 0)
            {
                Console.WriteLine("[Merge] No junction in the backbone was spanned by another k; keeping it as-is.");
                return false;
            }

            V_書き出し(p_出力パス, l_骨格名一覧, l_骨格配列, l_確定);
            Console.WriteLine($"[Merge] Joined {l_確定.Count} backbone junction(s) using sequence from other k values.");
            return true;
        }

        private static (List<string> A_名前, List<string> A_配列) Get_配列一覧(string p_パス)
        {
            List<string> l_名前 = [];
            List<string> l_配列 = [];
            using var l_読み込み = new FastaReader(p_パス);
            while (l_読み込み.Get_続きがあるか())
            {
                var l_項目 = l_読み込み.Get_次の配列();
                l_名前.Add(l_項目.A_ID);
                l_配列.Add(l_項目.A_配列);
            }
            return (l_名前, l_配列);
        }

        /// <summary>
        /// 骨格の各配列の両端について、アンカー k-mer から
        /// (配列番号, その k-mer の開始位置, 順鎖で一致したか) を引ける索引を作る。
        /// 複数の配列に現れる k-mer は行き先を一意に決められないため捨てる。
        /// </summary>
        private static Dictionary<UInt128, (int A_配列番号, int A_位置, bool A_順鎖か)> Get_骨格索引(
            List<string> p_骨格配列, int p_アンカーk長)
        {
            Dictionary<UInt128, (int, int, bool)> l_索引 = [];
            HashSet<UInt128> l_重複 = [];

            for (var l_番号 = 0; l_番号 < p_骨格配列.Count; l_番号++)
            {
                var l_配列 = p_骨格配列[l_番号];
                foreach (var l_位置 in Get_末端の位置範囲(l_配列.Length, p_アンカーk長))
                {
                    if (!KmerPacking.Get_正規化パック(l_配列, l_位置, p_アンカーk長, out var l_鍵))
                    {
                        continue;
                    }
                    if (l_重複.Contains(l_鍵))
                    {
                        continue;
                    }
                    if (l_索引.ContainsKey(l_鍵))
                    {
                        _ = l_索引.Remove(l_鍵);
                        _ = l_重複.Add(l_鍵);
                        continue;
                    }
                    l_索引[l_鍵] = (l_番号, l_位置, Get_順鎖か(l_配列, l_位置, p_アンカーk長));
                }
            }
            return l_索引;
        }

        private static IEnumerable<int> Get_末端の位置範囲(int p_配列長, int p_アンカーk長)
        {
            var l_最終位置 = p_配列長 - p_アンカーk長;
            if (l_最終位置 < 0)
            {
                yield break;
            }
            var l_先頭の終わり = Math.Min(l_最終位置, 末端とみなす長さ);
            for (var i = 0; i <= l_先頭の終わり; i++)
            {
                yield return i;
            }
            var l_末尾の始まり = Math.Max(l_先頭の終わり + 1, l_最終位置 - 末端とみなす長さ);
            for (var i = l_末尾の始まり; i <= l_最終位置; i++)
            {
                yield return i;
            }
        }

        private static bool Get_順鎖か(string p_配列, int p_位置, int p_k長)
        {
            _ = KmerPacking.Get_正規化パック(p_配列, p_位置, p_k長, out var l_正規形);
            UInt128 l_順鎖 = 0;
            for (var i = 0; i < p_k長; i++)
            {
                l_順鎖 = (l_順鎖 << 2) | (UInt128)(Util.Get_塩基ID(p_配列[p_位置 + i]) - 1);
            }
            return l_順鎖 == l_正規形;
        }

        /// <summary>
        /// 他の k の配列を走査し、骨格の異なる2本の末端を順に跨いでいるものを
        /// 隣接の証拠として拾う。
        /// </summary>
        private static List<橋渡し候補> Get_橋渡し候補(
            アセンブリ実行結果 p_他,
            Dictionary<UInt128, (int A_配列番号, int A_位置, bool A_順鎖か)> p_索引,
            List<string> p_骨格配列,
            int p_アンカーk長)
        {
            List<橋渡し候補> l_結果 = [];
            var (_, l_他の配列) = Get_配列一覧(p_他.A_最終パス);

            foreach (var l_配列 in l_他の配列)
            {
                // この配列が骨格のどこにどの順で当たったかを並べる。
                List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_同じ向きか)> l_当たり = [];
                for (var i = 0; i + p_アンカーk長 <= l_配列.Length; i++)
                {
                    if (!KmerPacking.Get_正規化パック(l_配列, i, p_アンカーk長, out var l_鍵)
                        || !p_索引.TryGetValue(l_鍵, out var l_骨格側))
                    {
                        continue;
                    }
                    l_当たり.Add((i, l_骨格側.A_配列番号, l_骨格側.A_位置,
                        Get_順鎖か(l_配列, i, p_アンカーk長) == l_骨格側.A_順鎖か));
                }

                l_結果.AddRange(Get_連続する2本の跨ぎ(l_当たり, l_配列, p_骨格配列, p_アンカーk長, p_他.A_k長));
            }
            return l_結果;
        }

        /// <summary>
        /// 当たりの列を走査し、骨格配列が切り替わる箇所ごとに橋渡し候補を作る。
        /// 切り替わりの直前・直後の当たりが、それぞれの骨格配列の「出口」と
        /// 「入口」に当たっているときだけ隣接の証拠になる。
        /// </summary>
        private static IEnumerable<橋渡し候補> Get_連続する2本の跨ぎ(
            List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_同じ向きか)> p_当たり,
            string p_跨いだ配列,
            List<string> p_骨格配列,
            int p_アンカーk長,
            int p_由来のk長)
        {
            for (var i = 1; i < p_当たり.Count; i++)
            {
                var l_前 = p_当たり[i - 1];
                var l_後 = p_当たり[i];
                if (l_前.A_配列番号 == l_後.A_配列番号)
                {
                    continue;
                }

                // 前側は出口(その向きで見た末端)、後側は入口でなければならない。
                var l_前の頂点 = Get_出口頂点(l_前, p_骨格配列, p_アンカーk長);
                var l_後の頂点 = Get_入口頂点(l_後, p_骨格配列, p_アンカーk長);
                if (l_前の頂点 is not { } l_始点 || l_後の頂点 is not { } l_終点)
                {
                    continue;
                }

                // 跨いだ配列のうち、2つのアンカーに挟まれた部分が繋ぎ目になる。
                var l_開始 = l_前.A_自分の位置 + p_アンカーk長;
                var l_長さ = l_後.A_自分の位置 - l_開始;
                if (l_長さ < 0 || l_長さ > 橋渡し長の上限)
                {
                    continue;
                }

                yield return new 橋渡し候補(
                    l_始点, l_終点, p_跨いだ配列.Substring(l_開始, l_長さ), p_由来のk長);
            }
        }

        /// <summary>
        /// その当たりが骨格配列の出口側末端かを判定し、対応する頂点を返す。
        /// 順鎖で当たって末尾側にあるなら 2i、逆鎖で当たって先頭側にあるなら 2i+1。
        /// </summary>
        private static int? Get_出口頂点(
            (int A_自分の位置, int A_配列番号, int A_位置, bool A_同じ向きか) p_当たり,
            List<string> p_骨格配列, int p_アンカーk長)
        {
            var l_配列長 = p_骨格配列[p_当たり.A_配列番号].Length;
            var l_末尾からの距離 = l_配列長 - p_アンカーk長 - p_当たり.A_位置;
            if (p_当たり.A_同じ向きか)
            {
                return l_末尾からの距離 <= 末端とみなす長さ ? p_当たり.A_配列番号 << 1 : null;
            }
            return p_当たり.A_位置 <= 末端とみなす長さ ? (p_当たり.A_配列番号 << 1) | 1 : null;
        }

        private static int? Get_入口頂点(
            (int A_自分の位置, int A_配列番号, int A_位置, bool A_同じ向きか) p_当たり,
            List<string> p_骨格配列, int p_アンカーk長)
        {
            var l_配列長 = p_骨格配列[p_当たり.A_配列番号].Length;
            var l_末尾からの距離 = l_配列長 - p_アンカーk長 - p_当たり.A_位置;
            if (p_当たり.A_同じ向きか)
            {
                return p_当たり.A_位置 <= 末端とみなす長さ ? p_当たり.A_配列番号 << 1 : null;
            }
            return l_末尾からの距離 <= 末端とみなす長さ ? (p_当たり.A_配列番号 << 1) | 1 : null;
        }

        /// <summary>
        /// 相互一意な橋渡しだけを残す。始点から見て行き先が1つに定まり、かつ
        /// 終点から見た来訪元も1つに定まるときだけ採用する。これを課さないと、
        /// 同じ行き先を指す複数の候補のうち先着だけが繋がれ、残りが根拠なく落ちる。
        /// </summary>
        private static Dictionary<int, 橋渡し候補> Get_相互一意な橋渡し(
            List<橋渡し候補> p_候補, int p_骨格の本数, int p_必要な独立支持数)
        {
            // 同じ隣接を何個の k が独立に主張しているかを数える。
            Dictionary<(int, int), HashSet<int>> l_支持したk = [];
            foreach (var l_候補 in p_候補)
            {
                V_数える(l_支持したk, (l_候補.A_始点, l_候補.A_終点), l_候補.A_由来のk長);
                V_数える(l_支持したk, (l_候補.A_終点 ^ 1, l_候補.A_始点 ^ 1), l_候補.A_由来のk長);
            }

            // 頂点ごとに、行き先の集合と代表となる候補を集める。
            Dictionary<int, HashSet<int>> l_行き先 = [];
            Dictionary<(int, int), 橋渡し候補> l_代表 = [];

            foreach (var l_候補 in p_候補)
            {
                V_登録(l_行き先, l_代表, l_候補);
                // 逆鎖側の双子も同じ隣接を表す。橋渡し配列も逆相補にする。
                V_登録(l_行き先, l_代表, new 橋渡し候補(
                    l_候補.A_終点 ^ 1,
                    l_候補.A_始点 ^ 1,
                    Util.V_逆相補_曖昧塩基あり(l_候補.A_橋渡し配列),
                    l_候補.A_由来のk長));
            }

            Dictionary<int, 橋渡し候補> l_確定 = [];
            foreach (var (l_始点, l_集合) in l_行き先)
            {
                if (l_集合.Count != 1)
                {
                    continue;
                }
                var l_終点 = l_集合.First();
                if (l_支持したk[(l_始点, l_終点)].Count < p_必要な独立支持数)
                {
                    continue;
                }
                if (!l_行き先.TryGetValue(l_終点 ^ 1, out var l_逆側) || l_逆側.Count != 1
                    || l_逆側.First() != (l_始点 ^ 1))
                {
                    continue;
                }
                if (l_始点 >> 1 >= p_骨格の本数 || l_終点 >> 1 >= p_骨格の本数)
                {
                    continue;
                }
                l_確定[l_始点] = l_代表[(l_始点, l_終点)];
            }
            return l_確定;
        }

        private static void V_数える(
            Dictionary<(int, int), HashSet<int>> p_支持したk, (int, int) p_辺, int p_k長)
        {
            if (!p_支持したk.TryGetValue(p_辺, out var l_集合))
            {
                l_集合 = [];
                p_支持したk[p_辺] = l_集合;
            }
            _ = l_集合.Add(p_k長);
        }

        private static void V_登録(
            Dictionary<int, HashSet<int>> p_行き先,
            Dictionary<(int, int), 橋渡し候補> p_代表,
            橋渡し候補 p_候補)
        {
            if (!p_行き先.TryGetValue(p_候補.A_始点, out var l_集合))
            {
                l_集合 = [];
                p_行き先[p_候補.A_始点] = l_集合;
            }
            _ = l_集合.Add(p_候補.A_終点);
            p_代表.TryAdd((p_候補.A_始点, p_候補.A_終点), p_候補);
        }

        /// <summary>
        /// 確定した橋渡しに沿って骨格配列を連結し、書き出す。
        /// 各配列はちょうど1回だけ使う。
        /// </summary>
        private static void V_書き出し(
            string p_出力パス,
            List<string> p_骨格名一覧,
            List<string> p_骨格配列,
            Dictionary<int, 橋渡し候補> p_確定)
        {
            var l_使用済み = new bool[p_骨格配列.Count];
            var l_ID = 1;
            using var l_書き込み = new FastaWriter(p_出力パス);

            for (var l_番号 = 0; l_番号 < p_骨格配列.Count; l_番号++)
            {
                if (l_使用済み[l_番号] || Get_来訪元があるか(p_確定, l_番号))
                {
                    continue;
                }
                l_書き込み.V_書き込み(
                    Get_名前(p_骨格名一覧, l_番号, l_ID++),
                    Get_連結配列(l_番号 << 1, p_骨格配列, p_確定, l_使用済み));
            }

            // 環状に閉じていて起点が見つからなかったぶんを拾う。
            for (var l_番号 = 0; l_番号 < p_骨格配列.Count; l_番号++)
            {
                if (!l_使用済み[l_番号])
                {
                    l_書き込み.V_書き込み(
                        Get_名前(p_骨格名一覧, l_番号, l_ID++),
                        Get_連結配列(l_番号 << 1, p_骨格配列, p_確定, l_使用済み));
                }
            }
        }

        private static string Get_名前(List<string> p_骨格名一覧, int p_番号, int p_ID)
        {
            return p_番号 < p_骨格名一覧.Count ? p_骨格名一覧[p_番号] : $"MERGED{p_ID}";
        }

        /// <summary>
        /// その配列へ入ってくる確定辺があるか。あれば連結の起点にはしない。
        /// </summary>
        private static bool Get_来訪元があるか(Dictionary<int, 橋渡し候補> p_確定, int p_番号)
        {
            // v へ入る辺は、双子 v^1 から出る辺と同値。
            return p_確定.ContainsKey((p_番号 << 1) | 1);
        }

        private static string Get_連結配列(
            int p_開始頂点,
            List<string> p_骨格配列,
            Dictionary<int, 橋渡し候補> p_確定,
            bool[] p_使用済み)
        {
            var l_結果 = new StringBuilder();
            var l_頂点 = p_開始頂点;

            while (true)
            {
                var l_番号 = l_頂点 >> 1;
                if (p_使用済み[l_番号])
                {
                    break;
                }
                p_使用済み[l_番号] = true;

                var l_配列 = p_骨格配列[l_番号];
                _ = l_結果.Append((l_頂点 & 1) == 0 ? l_配列 : Util.V_逆相補_曖昧塩基あり(l_配列));

                if (!p_確定.TryGetValue(l_頂点, out var l_橋渡し))
                {
                    break;
                }
                _ = l_結果.Append(l_橋渡し.A_橋渡し配列);
                l_頂点 = l_橋渡し.A_終点;
            }
            return l_結果.ToString();
        }
    }
}
