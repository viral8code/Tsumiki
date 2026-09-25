using System.Text;
using Tsumiki.Cores.UnitigBuilding;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker のうち、確定した結合を辿って contig 配列を組み立てる部分
    /// </summary>
    internal partial class ContigMaker
    {
        #region プロパティ

        /// <summary>
        /// 書き出した contig が分岐のある継ぎ目で通った辺 ((k+1)-mer、両向き)
        /// </summary>
        public HashSet<string> A_分岐の継ぎ目 { get; } = new(StringComparer.Ordinal);

        #endregion

        #region 公開メソッド

        /// <summary>
        /// walk が分岐のある継ぎ目で通った辺を (k+1)-mer で返す
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_walk順"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        internal static IEnumerable<string> Get_分岐の継ぎ目(UnitigGraph p_グラフ, IReadOnlyList<string> p_unitig配列, IReadOnlyList<int> p_walk順, int p_k長)
        {
            for (var w = 1; w < p_walk順.Count; w++)
            {
                var l_前 = p_walk順[w - 1];
                var l_次 = p_walk順[w];
                if (p_グラフ.A_出辺[l_前].Count < 2 && p_グラフ.Get_入次数(l_次) < 2)
                {
                    continue;
                }

                var l_前の配列 = p_unitig配列[l_前];
                var l_次の配列 = p_unitig配列[l_次];
                if (l_前の配列.Length < p_k長 || l_次の配列.Length < p_k長)
                {
                    continue;
                }

                yield return string.Concat(l_前の配列.AsSpan(l_前の配列.Length - p_k長), l_次の配列.AsSpan(p_k長 - 1, 1));
            }
        }

        /// <summary>
        /// 符号付き unitig ID (正=順鎖、負=逆鎖) をグラフの頂点番号に変換する
        /// </summary>
        /// <param name="p_符号付きunitigID"></param>
        /// <returns></returns>
        public static int Get_頂点番号(int p_符号付きunitigID)
        {
            return (Math.Abs(p_符号付きunitigID) << 1) | (p_符号付きunitigID > 0 ? 0 : 1);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 始点から結合を辿って 1 本の contig を組み立て、結果を各一覧へ追加する
        /// </summary>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_結合"></param>
        /// <param name="p_訪問済み"></param>
        /// <param name="p_重なり長"></param>
        /// <param name="p_始点"></param>
        /// <param name="p_contig群"></param>
        /// <param name="p_walk順群"></param>
        /// <param name="p_環状フラグ群"></param>
        private static void V_実行_walk(List<string> p_unitig配列, int[] p_結合, bool[] p_訪問済み, int p_重なり長, int p_始点, List<string> p_contig群, List<List<int>> p_walk順群, List<bool> p_環状フラグ群)
        {
            List<int> l_walk順 = [];
            var (l_配列, l_Is環状) = Get_walk結果(p_unitig配列, p_結合, p_訪問済み, p_重なり長, p_始点, l_walk順);
            p_contig群.Add(l_配列);
            p_walk順群.Add(l_walk順);
            p_環状フラグ群.Add(l_Is環状);
        }

        /// <summary>
        /// 始点から結合を辿って配列を組み立てる
        /// </summary>
        /// <param name="p_unitig配列"></param>
        /// <param name="p_結合"></param>
        /// <param name="p_訪問済み"></param>
        /// <param name="p_重なり長"></param>
        /// <param name="p_始点"></param>
        /// <param name="p_walk順"></param>
        /// <returns></returns>
        private static (string A_配列, bool A_Is環状) Get_walk結果(List<string> p_unitig配列, int[] p_結合, bool[] p_訪問済み, int p_重なり長, int p_始点, List<int> p_walk順)
        {
            var l_出力 = new StringBuilder(p_unitig配列[p_始点]);
            p_walk順.Add(p_始点);
            p_訪問済み[p_始点 >> 1] = true;
            var l_現在 = p_始点;
            var l_Is環状 = false;
            while (true)
            {
                var l_次 = p_結合[l_現在];

                if (l_次 < 0)
                {
                    break;
                }

                if (p_訪問済み[l_次 >> 1])
                {
                    l_Is環状 = l_次 == p_始点;
                    break;
                }

                var l_配列 = p_unitig配列[l_次];

                if (l_配列.Length < p_重なり長 || l_出力.Length < p_重なり長)
                {
                    break;
                }

                if (!Is重なり一致(l_出力, l_配列, p_重なり長))
                {
                    break;
                }

                _ = l_出力.Append(l_配列[p_重なり長..]);
                p_訪問済み[l_次 >> 1] = true;
                p_walk順.Add(l_次);
                l_現在 = l_次;
            }

            if (l_Is環状 && l_出力.Length > p_重なり長)
            {
                _ = l_出力.Remove(l_出力.Length - p_重なり長, p_重なり長);
            }

            return (l_出力.ToString(), l_Is環状);
        }

        /// <summary>
        /// 出力文字列末尾と unitig の文字列先頭が指定された長さだけ一致するかどうかを判定する
        /// </summary>
        /// <param name="p_出力"></param>
        /// <param name="p_unitig"></param>
        /// <param name="p_重なり長"></param>
        /// <returns></returns>
        private static bool Is重なり一致(StringBuilder p_出力, string p_unitig, int p_重なり長)
        {
            var l_開始位置 = p_出力.Length - p_重なり長;
            for (var j = 0; j < p_重なり長; j++)
            {
                if (p_出力[l_開始位置 + j] != p_unitig[j])
                {
                    return false;
                }
            }

            return true;
        }

        #endregion
    }
}
