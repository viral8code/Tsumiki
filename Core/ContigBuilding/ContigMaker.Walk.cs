using System.Text;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker のうち、確定した結合を辿って contig 配列を組み立てる部分。
    /// (unitig へのマッピングは ContigMaker.Mapping.cs、辺の選択・結合の確定は
    /// ContigMaker.cs の V_結合_コンティグ を参照)
    /// </summary>
    internal partial class ContigMaker
    {
        /// <summary>
        /// 符号付き unitig ID(正=順鎖、負=逆鎖)をグラフの頂点番号に変換する。
        /// </summary>
        internal static int Get_頂点番号(int p_符号付きユニティグID)
        {
            return (Math.Abs(p_符号付きユニティグID) << 1) | (p_符号付きユニティグID > 0 ? 0 : 1);
        }

        /// <summary>
        /// 始点から結合を辿って1本の contig を組み立て、結果を各一覧へ追加する。
        /// </summary>
        private static void V_実行_walk(
            List<string> p_ユニティグ配列, int[] p_結合, bool[] p_訪問済み, int p_重なり長, int p_始点,
            List<string> p_コンティグ群, List<List<int>> p_walk順群, List<bool> p_環状フラグ群)
        {
            List<int> l_walk順 = [];
            var (l_配列, l_環状か) = Get_walk結果(p_ユニティグ配列, p_結合, p_訪問済み, p_重なり長, p_始点, l_walk順);
            p_コンティグ群.Add(l_配列);
            p_walk順群.Add(l_walk順);
            p_環状フラグ群.Add(l_環状か);
        }

        /// <summary>
        /// 始点から結合を辿って配列を組み立てる。経路が始点へ戻ってきた場合は
        /// 環状として報告する。
        /// </summary>
        private static (string A_配列, bool A_環状か) Get_walk結果(
            List<string> p_ユニティグ配列, int[] p_結合, bool[] p_訪問済み, int p_重なり長, int p_始点,
            List<int> p_walk順)
        {
            var l_出力 = new StringBuilder(p_ユニティグ配列[p_始点]);
            p_walk順.Add(p_始点);
            p_訪問済み[p_始点 >> 1] = true;
            var l_現在 = p_始点;
            var l_環状か = false;
            while (true)
            {
                var l_次 = p_結合[l_現在];
                if (l_次 < 0)
                {
                    break;
                }
                if (p_訪問済み[l_次 >> 1])
                {
                    // 始点へ戻ってきた = 経路が閉じている。細菌の染色体と
                    // プラスミドは環状なので、これは「その複製単位を
                    // 完全に1周組み上げられた」ことを意味する。
                    l_環状か = l_次 == p_始点;
                    break;
                }
                var l_配列 = p_ユニティグ配列[l_次];
                if (l_配列.Length < p_重なり長 || l_出力.Length < p_重なり長)
                {
                    break;
                }
                // 構築方法より k-1 のオーバーラップは保証されているが、
                // 万一崩れていた場合に誤った配列を作らないよう検証する。
                if (!Get_重なりが一致するか(l_出力, l_配列, p_重なり長))
                {
                    break;
                }
                _ = l_出力.Append(l_配列[p_重なり長..]);
                p_訪問済み[l_次 >> 1] = true;
                p_walk順.Add(l_次);
                l_現在 = l_次;
            }

            // 環状では末尾 unitig が始点との重なり k-1 塩基を含んでおり、
            // それは配列の先頭にも現れる。線状の連結では次の unitig 側から
            // 取り除くが、環状では「次」が出力済みの始点なので末尾から取り除く。
            if (l_環状か && l_出力.Length > p_重なり長)
            {
                _ = l_出力.Remove(l_出力.Length - p_重なり長, p_重なり長);
            }

            return (l_出力.ToString(), l_環状か);
        }

        /// <summary>
        /// 出力の末尾 p_重なり長 文字と unitig の先頭 p_重なり長 文字が
        /// 一致するかどうかを判定する。
        /// </summary>
        private static bool Get_重なりが一致するか(StringBuilder p_出力, string p_ユニティグ, int p_重なり長)
        {
            var l_開始位置 = p_出力.Length - p_重なり長;
            for (var j = 0; j < p_重なり長; j++)
            {
                if (p_出力[l_開始位置 + j] != p_ユニティグ[j])
                {
                    return false;
                }
            }
            return true;
        }
    }
}
