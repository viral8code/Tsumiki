using System.Text;

namespace Tsumiki.Core
{
    /// <summary>
    /// ContigMaker のうち、確定した結合を辿って contig 配列を組み立てる部分
    /// </summary>
    /// <remarks>
    /// (unitig へのマッピングは ContigMaker.Mapping.cs、辺の選択・結合の確定は ContigMaker.cs の <see cref="V_結合_Contig"/> を参照)
    /// </remarks>
    internal partial class ContigMaker
    {
        #region 公開メソッド

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
        /// <remarks>
        /// 経路が始点へ戻ってきた場合は環状として報告する
        /// </remarks>
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
                    // 始点へ戻ってきた = 経路が閉じている
                    // 細菌の染色体と
                    // プラスミドは環状なので、これは「その複製単位を
                    // 完全に 1 周組み上げられた」ことを意味する
                    l_Is環状 = l_次 == p_始点;
                    break;
                }
                var l_配列 = p_unitig配列[l_次];

                if (l_配列.Length < p_重なり長 || l_出力.Length < p_重なり長)
                {
                    break;
                }

                // 構築方法より k-1 のオーバーラップは保証されているが、
                // 万一崩れていた場合に誤った配列を作らないよう検証する
                if (!Is重なり一致(l_出力, l_配列, p_重なり長))
                {
                    break;
                }
                _ = l_出力.Append(l_配列[p_重なり長..]);
                p_訪問済み[l_次 >> 1] = true;
                p_walk順.Add(l_次);
                l_現在 = l_次;
            }

            // 環状では末尾 unitig が始点との重なり k-1 塩基を含んでおり、
            // それは配列の先頭にも現れる
            // 線状の連結では次の unitig 側から
            // 取り除くが、環状では「次」が出力済みの始点なので末尾から取り除く
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
