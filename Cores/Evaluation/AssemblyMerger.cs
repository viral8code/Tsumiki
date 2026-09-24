using System.Text;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 複数の k のアセンブリを 1 つに統合する
    /// </summary>
    internal static class AssemblyMerger
    {
        #region 定数

        /// <summary>
        /// 骨格側の末端でアンカーを探す長さ
        /// </summary>
        private const int 末端とみなす長さ = 2_000;

        /// <summary>
        /// 橋渡しとして認める最大の挟み込み長
        /// </summary>
        private const int 橋渡し長の上限 = 50_000;

        /// <summary>
        /// 連結を認めるために必要な、独立に同じ隣接を主張した k の数
        /// </summary>
        private const int 必要な独立支持数の既定値 = 2;

        /// <summary>
        /// アンカーから骨格の末端までの区間を跨いだ配列と照合するときに許す不一致率
        /// </summary>
        private const double 末端照合の許容不一致率 = 0.01D;

        /// <summary>
        /// 隣接の証拠に使う当たりの塊に要求する、連続した当たりの数
        /// </summary>
        private const int 塊とみなす当たり数 = 20;

        /// <summary>
        /// 接合の端点に使う骨格配列の最小長
        /// </summary>
        private const int 端点に使う最小長 = 500;

        /// <summary>
        /// 繋ぎ目の k-mer のうち、骨格に既にあってよい割合
        /// </summary>
        private const double 繋ぎ目に許す既知kmerの割合 = 0.05D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 骨格に対して他の候補を統合し、結果を p_出力パス へ書き出す
        /// </summary>
        /// <param name="p_骨格"></param>
        /// <param name="p_全候補"></param>
        /// <param name="p_アンカーk長"></param>
        /// <param name="p_出力パス"></param>
        /// <param name="p_必要な独立支持数"></param>
        /// <returns></returns>
        public static bool Try統合(アセンブリ実行結果 p_骨格, IReadOnlyList<アセンブリ実行結果> p_全候補, int p_アンカーk長, string p_出力パス, int p_必要な独立支持数 = 必要な独立支持数の既定値)
        {
            var (l_骨格名一覧, l_骨格配列) = Get_配列一覧(p_骨格.A_最終パス);
            if (l_骨格配列.Count == 0)
            {
                return false;
            }

            var l_索引 = Get_骨格索引(l_骨格配列, p_アンカーk長);
            var l_骨格kmer = Get_骨格kmer集合(l_骨格配列, p_アンカーk長);
            var l_候補 = new List<橋渡し候補>();
            foreach (var l_他 in p_全候補)
            {
                if (l_他.A_k長 == p_骨格.A_k長)
                {
                    continue;
                }
                l_候補.AddRange(Get_橋渡し候補(l_他, l_索引, l_骨格配列, p_アンカーk長).Where(x => !Has既知の配列(x.A_橋渡し配列, l_骨格kmer, p_アンカーk長)));
            }

            var l_確定 = Get_相互一意な橋渡し(l_候補, l_骨格配列.Count, p_必要な独立支持数);
            if (l_確定.Count == 0)
            {
                Logger.V_出力(メッセージID.統合できる接合点なし);
                return false;
            }

            V_書き出し(p_出力パス, l_骨格名一覧, l_骨格配列, l_確定);
            Logger.V_出力(メッセージID.統合した接合点数, l_確定.Count);
            return true;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// FASTA から名前と配列の一覧を読み込む
        /// </summary>
        /// <param name="p_パス"></param>
        /// <returns></returns>
        private static (List<string> A_名前, List<string> A_配列) Get_配列一覧(string p_パス)
        {
            List<string> l_名前 = [];
            List<string> l_配列 = [];
            using var l_読み込み = new FastaReader(p_パス);
            while (l_読み込み.Has続き())
            {
                var l_項目 = l_読み込み.Get_次の配列();
                l_名前.Add(l_項目.A_ID);
                l_配列.Add(l_項目.A_配列);
            }
            return (l_名前, l_配列);
        }

        /// <summary>
        /// 骨格の各配列の両端について、アンカー k-mer から (配列番号, その k-mer の開始位置, 順鎖で一致したか) を引ける索引を作る
        /// </summary>
        /// <param name="p_骨格配列"></param>
        /// <param name="p_アンカーk長"></param>
        /// <returns></returns>
        private static Dictionary<UInt128, (int A_配列番号, int A_位置, bool A_Is順鎖)> Get_骨格索引(List<string> p_骨格配列, int p_アンカーk長)
        {
            Dictionary<UInt128, (int, int, bool)> l_索引 = [];
            HashSet<UInt128> l_重複 = [];

            for (var l_番号 = 0; l_番号 < p_骨格配列.Count; l_番号++)
            {
                var l_配列 = p_骨格配列[l_番号];
                if (l_配列.Length < 端点に使う最小長)
                {
                    continue;
                }
                foreach (var l_位置 in Get_末端位置範囲(l_配列.Length, p_アンカーk長))
                {
                    if (!KmerPacking.TryGet_正規化パック(l_配列, l_位置, p_アンカーk長, out var l_鍵))
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
                    l_索引[l_鍵] = (l_番号, l_位置, Is順鎖(l_配列, l_位置, p_アンカーk長));
                }
            }
            return l_索引;
        }

        /// <summary>
        /// 骨格配列に現れる全 k-mer の正規形
        /// </summary>
        /// <param name="p_骨格配列"></param>
        /// <param name="p_アンカーk長"></param>
        /// <returns></returns>
        private static HashSet<UInt128> Get_骨格kmer集合(List<string> p_骨格配列, int p_アンカーk長)
        {
            HashSet<UInt128> l_集合 = [];
            foreach (var l_配列 in p_骨格配列)
            {
                for (var i = 0; i + p_アンカーk長 <= l_配列.Length; i++)
                {
                    if (KmerPacking.TryGet_正規化パック(l_配列, i, p_アンカーk長, out var l_鍵))
                    {
                        _ = l_集合.Add(l_鍵);
                    }
                }
            }
            return l_集合;
        }

        /// <summary>
        /// 繋ぎ目が、骨格に既にある配列を許容量を超えて含むか
        /// </summary>
        /// <param name="p_橋渡し配列"></param>
        /// <param name="p_骨格kmer"></param>
        /// <param name="p_アンカーk長"></param>
        /// <returns></returns>
        private static bool Has既知の配列(string p_橋渡し配列, HashSet<UInt128> p_骨格kmer, int p_アンカーk長)
        {
            var l_件数 = 0;
            var l_既知 = 0;
            for (var i = 0; i + p_アンカーk長 <= p_橋渡し配列.Length; i++)
            {
                if (!KmerPacking.TryGet_正規化パック(p_橋渡し配列, i, p_アンカーk長, out var l_鍵))
                {
                    continue;
                }
                l_件数++;
                if (p_骨格kmer.Contains(l_鍵))
                {
                    l_既知++;
                }
            }
            return l_既知 > p_アンカーk長 && l_既知 > l_件数 * 繋ぎ目に許す既知kmerの割合;
        }

        /// <summary>
        /// 配列の両端で、アンカーを取る位置を返す
        /// </summary>
        /// <param name="p_配列長">配列の長さ</param>
        /// <param name="p_アンカーk長">アンカーの k 長</param>
        /// <returns>アンカーを取る位置</returns>
        private static IEnumerable<int> Get_末端位置範囲(int p_配列長, int p_アンカーk長)
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

        /// <summary>
        /// その位置の k-mer が順鎖のまま正規形になるか
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_位置">k-mer の開始位置</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>順鎖なら true</returns>
        private static bool Is順鎖(string p_配列, int p_位置, int p_k長)
        {
            _ = KmerPacking.TryGet_正規化パック(p_配列, p_位置, p_k長, out var l_正規形);
            UInt128 l_順鎖 = 0;
            for (var i = 0; i < p_k長; i++)
            {
                l_順鎖 = (l_順鎖 << 2) | (UInt128)(Util.Get_塩基ID(p_配列[p_位置 + i]) - 1);
            }
            return l_順鎖 == l_正規形;
        }

        /// <summary>
        /// 他の k の配列を走査し、骨格の異なる 2 本の末端を順に跨いでいるものを隣接の証拠として拾う
        /// </summary>
        /// <param name="p_他"></param>
        /// <param name="p_索引"></param>
        /// <param name="p_骨格配列"></param>
        /// <param name="p_アンカーk長"></param>
        /// <returns></returns>
        private static List<橋渡し候補> Get_橋渡し候補(アセンブリ実行結果 p_他, Dictionary<UInt128, (int A_配列番号, int A_位置, bool A_Is順鎖)> p_索引, List<string> p_骨格配列, int p_アンカーk長)
        {
            List<橋渡し候補> l_結果 = [];
            var (_, l_他の配列) = Get_配列一覧(p_他.A_最終パス);

            foreach (var l_配列 in l_他の配列)
            {
                List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)> l_当たり = [];
                for (var i = 0; i + p_アンカーk長 <= l_配列.Length; i++)
                {
                    if (!KmerPacking.TryGet_正規化パック(l_配列, i, p_アンカーk長, out var l_鍵)
                        || !p_索引.TryGetValue(l_鍵, out var l_骨格側))
                    {
                        continue;
                    }
                    l_当たり.Add((i, l_骨格側.A_配列番号, l_骨格側.A_位置, Is順鎖(l_配列, i, p_アンカーk長) == l_骨格側.A_Is順鎖));
                }

                l_結果.AddRange(Get_連続2本跨ぎ(Get_塊として続く当たり(l_当たり), l_配列, p_骨格配列, p_アンカーk長, p_他.A_k長));
            }
            return l_結果;
        }

        /// <summary>
        /// 同じ骨格配列の同じ対角線上に一定数以上続く当たりだけを残す
        /// </summary>
        /// <param name="p_当たり">跨いだ配列の位置順に並んだ当たり</param>
        /// <returns></returns>
        private static List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)> Get_塊として続く当たり(List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)> p_当たり)
        {
            List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)> l_結果 = [];
            var l_塊の開始 = 0;
            for (var i = 1; i <= p_当たり.Count; i++)
            {
                if (i < p_当たり.Count && Is同じ塊(p_当たり[i - 1], p_当たり[i]))
                {
                    continue;
                }

                if (i - l_塊の開始 >= 塊とみなす当たり数)
                {
                    l_結果.AddRange(p_当たり.Skip(l_塊の開始).Take(i - l_塊の開始));
                }
                l_塊の開始 = i;
            }
            return l_結果;
        }

        /// <summary>
        /// 2 つの当たりが、同じ骨格配列の同じ向き・同じ対角線上にあるか
        /// </summary>
        /// <param name="p_前"></param>
        /// <param name="p_後"></param>
        /// <returns></returns>
        private static bool Is同じ塊((int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向) p_前, (int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向) p_後)
        {
            return p_前.A_配列番号 == p_後.A_配列番号
                && p_前.A_Is同方向 == p_後.A_Is同方向
                && (p_前.A_Is同方向 ? p_前.A_自分の位置 - p_前.A_位置 == p_後.A_自分の位置 - p_後.A_位置 : p_前.A_自分の位置 + p_前.A_位置 == p_後.A_自分の位置 + p_後.A_位置);
        }

        /// <summary>
        /// 当たりの列を走査し、骨格配列が切り替わる箇所ごとに橋渡し候補を作る
        /// </summary>
        /// <param name="p_当たり"></param>
        /// <param name="p_跨いだ配列"></param>
        /// <param name="p_骨格配列"></param>
        /// <param name="p_アンカーk長"></param>
        /// <param name="p_由来のk長"></param>
        /// <returns></returns>
        private static IEnumerable<橋渡し候補> Get_連続2本跨ぎ(List<(int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向)> p_当たり, string p_跨いだ配列, List<string> p_骨格配列, int p_アンカーk長, int p_由来のk長)
        {
            for (var i = 1; i < p_当たり.Count; i++)
            {
                var l_前 = p_当たり[i - 1];
                var l_後 = p_当たり[i];

                if (l_前.A_配列番号 == l_後.A_配列番号)
                {
                    continue;
                }

                var l_前の頂点 = Get_出口頂点(l_前, p_骨格配列, p_アンカーk長);
                var l_後の頂点 = Get_入口頂点(l_後, p_骨格配列, p_アンカーk長);

                if (l_前の頂点 is not { } l_始点 || l_後の頂点 is not { } l_終点)
                {
                    continue;
                }

                var l_前配列 = Get_向き付き配列(p_骨格配列, l_始点);
                var l_後配列 = Get_向き付き配列(p_骨格配列, l_終点);
                var l_前の残り = l_前.A_Is同方向 ? l_前配列.Length - p_アンカーk長 - l_前.A_位置 : l_前.A_位置;
                var l_後の手前 = l_後.A_Is同方向 ? l_後.A_位置 : l_後配列.Length - p_アンカーk長 - l_後.A_位置;

                var l_開始 = l_前.A_自分の位置 + p_アンカーk長 + l_前の残り;
                var l_終了 = l_後.A_自分の位置 - l_後の手前;
                if (l_開始 > p_跨いだ配列.Length || l_終了 < 0)
                {
                    continue;
                }

                if (!Is一致(p_跨いだ配列, l_前.A_自分の位置, l_前配列, l_前配列.Length - p_アンカーk長 - l_前の残り, p_アンカーk長 + l_前の残り)
                    || !Is一致(p_跨いだ配列, l_終了, l_後配列, 0, l_後の手前 + p_アンカーk長))
                {
                    continue;
                }

                var l_長さ = l_終了 - l_開始;
                if (l_長さ > 橋渡し長の上限)
                {
                    continue;
                }

                if (l_長さ >= 0)
                {
                    yield return new 橋渡し候補(l_始点, l_終点, p_跨いだ配列.Substring(l_開始, l_長さ), p_由来のk長);
                    continue;
                }

                var l_重なり長 = -l_長さ;
                if (l_重なり長 >= l_前配列.Length || l_重なり長 >= l_後配列.Length)
                {
                    continue;
                }
                yield return new 橋渡し候補(l_始点, l_終点, string.Empty, p_由来のk長, l_重なり長);
            }
        }

        /// <summary>
        /// 頂点の向きで見た骨格配列
        /// </summary>
        /// <param name="p_骨格配列"></param>
        /// <param name="p_頂点"></param>
        /// <returns></returns>
        private static string Get_向き付き配列(List<string> p_骨格配列, int p_頂点)
        {
            var l_配列 = p_骨格配列[p_頂点 >> 1];
            return (p_頂点 & 1) == 0 ? l_配列 : Util.V_逆相補_曖昧塩基あり(l_配列);
        }

        /// <summary>
        /// 跨いだ配列の区間が、骨格配列の区間と許容範囲内で一致するか
        /// </summary>
        /// <param name="p_跨いだ配列"></param>
        /// <param name="p_跨いだ側の開始"></param>
        /// <param name="p_骨格側配列"></param>
        /// <param name="p_骨格側の開始"></param>
        /// <param name="p_長さ"></param>
        /// <returns></returns>
        private static bool Is一致(string p_跨いだ配列, int p_跨いだ側の開始, string p_骨格側配列, int p_骨格側の開始, int p_長さ)
        {
            if (p_跨いだ側の開始 < 0 || p_骨格側の開始 < 0 || p_跨いだ側の開始 + p_長さ > p_跨いだ配列.Length || p_骨格側の開始 + p_長さ > p_骨格側配列.Length)
            {
                return false;
            }

            var l_許容不一致数 = (int)(p_長さ * 末端照合の許容不一致率);
            var l_不一致数 = 0;
            for (var i = 0; i < p_長さ; i++)
            {
                if (p_跨いだ配列[p_跨いだ側の開始 + i] != p_骨格側配列[p_骨格側の開始 + i] && ++l_不一致数 > l_許容不一致数)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// その当たりが骨格配列の出口側末端かを判定し、対応する頂点を返す
        /// </summary>
        /// <param name="p_当たり"></param>
        /// <param name="p_骨格配列"></param>
        /// <param name="p_アンカーk長"></param>
        /// <returns></returns>
        private static int? Get_出口頂点((int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向) p_当たり, List<string> p_骨格配列, int p_アンカーk長)
        {
            var l_配列長 = p_骨格配列[p_当たり.A_配列番号].Length;
            var l_末尾からの距離 = l_配列長 - p_アンカーk長 - p_当たり.A_位置;
            return p_当たり.A_Is同方向 ? l_末尾からの距離 <= 末端とみなす長さ ? p_当たり.A_配列番号 << 1 : null : p_当たり.A_位置 <= 末端とみなす長さ ? (p_当たり.A_配列番号 << 1) | 1 : null;
        }

        /// <summary>
        /// アンカーの当たりから、骨格側で入っていく頂点を返す
        /// </summary>
        /// <param name="p_当たり">アンカーが当たった位置と向き</param>
        /// <param name="p_骨格配列">骨格の配列</param>
        /// <param name="p_アンカーk長">アンカーの k 長</param>
        /// <returns>入口の頂点、決められなければ null</returns>
        private static int? Get_入口頂点((int A_自分の位置, int A_配列番号, int A_位置, bool A_Is同方向) p_当たり, List<string> p_骨格配列, int p_アンカーk長)
        {
            var l_配列長 = p_骨格配列[p_当たり.A_配列番号].Length;
            var l_末尾からの距離 = l_配列長 - p_アンカーk長 - p_当たり.A_位置;
            return p_当たり.A_Is同方向 ? p_当たり.A_位置 <= 末端とみなす長さ ? p_当たり.A_配列番号 << 1 : null : l_末尾からの距離 <= 末端とみなす長さ ? (p_当たり.A_配列番号 << 1) | 1 : null;
        }

        /// <summary>
        /// 相互一意な橋渡しだけを残す
        /// </summary>
        /// <param name="p_候補"></param>
        /// <param name="p_骨格の本数"></param>
        /// <param name="p_必要な独立支持数"></param>
        /// <returns></returns>
        private static Dictionary<int, 橋渡し候補> Get_相互一意な橋渡し(List<橋渡し候補> p_候補, int p_骨格の本数, int p_必要な独立支持数)
        {
            Dictionary<(int, int), HashSet<int>> l_支持したk = [];
            foreach (var l_候補 in p_候補)
            {
                V_登録_k支持(l_支持したk, (l_候補.A_始点, l_候補.A_終点), l_候補.A_由来のk長);
                V_登録_k支持(l_支持したk, (l_候補.A_終点 ^ 1, l_候補.A_始点 ^ 1), l_候補.A_由来のk長);
            }

            Dictionary<int, HashSet<int>> l_行き先 = [];
            Dictionary<(int, int), 橋渡し候補> l_代表 = [];

            foreach (var l_候補 in p_候補)
            {
                V_登録(l_行き先, l_代表, l_候補);

                V_登録(l_行き先, l_代表, new 橋渡し候補(l_候補.A_終点 ^ 1, l_候補.A_始点 ^ 1, Util.V_逆相補_曖昧塩基あり(l_候補.A_橋渡し配列), l_候補.A_由来のk長, l_候補.A_重なり長));
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

        /// <summary>
        /// 辺を支持した k を書き留める
        /// </summary>
        /// <param name="p_支持したk">辺ごとに支持した k</param>
        /// <param name="p_辺">支持された辺</param>
        /// <param name="p_k長">支持した k</param>
        private static void V_登録_k支持(Dictionary<(int, int), HashSet<int>> p_支持したk, (int, int) p_辺, int p_k長)
        {
            if (!p_支持したk.TryGetValue(p_辺, out var l_集合))
            {
                l_集合 = [];
                p_支持したk[p_辺] = l_集合;
            }
            _ = l_集合.Add(p_k長);
        }

        /// <summary>
        /// 橋渡しの候補を、行き先と代表の両方へ書き留める
        /// </summary>
        /// <param name="p_行き先">頂点ごとの行き先</param>
        /// <param name="p_代表">辺ごとの代表となる候補</param>
        /// <param name="p_候補">書き留める候補</param>
        private static void V_登録(Dictionary<int, HashSet<int>> p_行き先, Dictionary<(int, int), 橋渡し候補> p_代表, 橋渡し候補 p_候補)
        {
            if (!p_行き先.TryGetValue(p_候補.A_始点, out var l_集合))
            {
                l_集合 = [];
                p_行き先[p_候補.A_始点] = l_集合;
            }
            _ = l_集合.Add(p_候補.A_終点);
            _ = p_代表.TryAdd((p_候補.A_始点, p_候補.A_終点), p_候補);
        }

        /// <summary>
        /// 確定した橋渡しに沿って骨格配列を連結し、書き出す
        /// </summary>
        /// <param name="p_出力パス"></param>
        /// <param name="p_骨格名一覧"></param>
        /// <param name="p_骨格配列"></param>
        /// <param name="p_確定"></param>
        private static void V_書き出し(string p_出力パス, List<string> p_骨格名一覧, List<string> p_骨格配列, Dictionary<int, 橋渡し候補> p_確定)
        {
            var l_使用済み = new bool[p_骨格配列.Count];
            var l_ID = 1;
            using var l_書き込み = new FastaWriter(p_出力パス);

            for (var l_番号 = 0; l_番号 < p_骨格配列.Count; l_番号++)
            {
                if (l_使用済み[l_番号] || Has来訪元(p_確定, l_番号))
                {
                    continue;
                }
                l_書き込み.V_書き込み(Get_名前(p_骨格名一覧, l_番号, l_ID++), Get_連結配列(l_番号 << 1, p_骨格配列, p_確定, l_使用済み));
            }

            for (var l_番号 = 0; l_番号 < p_骨格配列.Count; l_番号++)
            {
                if (!l_使用済み[l_番号])
                {
                    l_書き込み.V_書き込み(Get_名前(p_骨格名一覧, l_番号, l_ID++), Get_連結配列(l_番号 << 1, p_骨格配列, p_確定, l_使用済み));
                }
            }
        }

        /// <summary>
        /// 統合後の配列に付ける名前を返す
        /// </summary>
        /// <param name="p_骨格名一覧">骨格側の名前</param>
        /// <param name="p_番号">骨格の番号</param>
        /// <param name="p_ID">通し番号</param>
        /// <returns>付ける名前</returns>
        private static string Get_名前(List<string> p_骨格名一覧, int p_番号, int p_ID)
        {
            return p_番号 < p_骨格名一覧.Count ? p_骨格名一覧[p_番号] : $"MERGED{p_ID}";
        }

        /// <summary>
        /// その配列へ入ってくる確定辺があるか
        /// </summary>
        /// <param name="p_確定"></param>
        /// <param name="p_番号"></param>
        /// <returns></returns>
        private static bool Has来訪元(Dictionary<int, 橋渡し候補> p_確定, int p_番号)
        {
            return p_確定.ContainsKey((p_番号 << 1) | 1);
        }

        /// <summary>
        /// 確定した橋渡しを辿って、繋がった 1 本の配列を返す
        /// </summary>
        /// <param name="p_開始頂点">辿り始める頂点</param>
        /// <param name="p_骨格配列">骨格の配列</param>
        /// <param name="p_確定">頂点ごとに確定した橋渡し</param>
        /// <param name="p_使用済み">既に使った配列</param>
        /// <returns>繋がった配列</returns>
        private static string Get_連結配列(int p_開始頂点, List<string> p_骨格配列, Dictionary<int, 橋渡し候補> p_確定, bool[] p_使用済み)
        {
            var l_結果 = new StringBuilder();
            var l_頂点 = p_開始頂点;
            var l_削る長さ = 0;

            while (true)
            {
                var l_番号 = l_頂点 >> 1;
                if (p_使用済み[l_番号])
                {
                    break;
                }
                p_使用済み[l_番号] = true;

                var l_配列 = Get_向き付き配列(p_骨格配列, l_頂点);
                _ = l_結果.Append(l_配列, Math.Min(l_削る長さ, l_配列.Length), l_配列.Length - Math.Min(l_削る長さ, l_配列.Length));

                if (!p_確定.TryGetValue(l_頂点, out var l_橋渡し))
                {
                    break;
                }
                _ = l_結果.Append(l_橋渡し.A_橋渡し配列);
                l_削る長さ = l_橋渡し.A_重なり長;
                l_頂点 = l_橋渡し.A_終点;
            }
            return l_結果.ToString();
        }

        #endregion
    }
}
