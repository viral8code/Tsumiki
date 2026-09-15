using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 反復解決の検証で使う unitig グラフを配列から組み立てる
    /// </summary>
    internal static class RepeatGraphFixture
    {
        #region 定数

        /// <summary>
        /// 曖昧 k-mer の番兵
        /// </summary>
        private const int 曖昧kmer番号 = int.MinValue;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ"></param>
        /// <param name="p_乱数種"></param>
        /// <returns></returns>
        public static string Get_乱数配列(int p_長さ, int p_乱数種)
        {
            var l_乱数生成器 = new Random(p_乱数種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数生成器.Next(4)]));
        }

        /// <summary>
        /// ContigMaker と同じ規則で k-mer 辞書を作り、unitig グラフを組み立てる
        /// </summary>
        /// <param name="p_k長"></param>
        /// <param name="p_unitig群">ID 1 から順に並べる unitig 配列</param>
        /// <returns></returns>
        public static (UnitigGraph A_グラフ, List<string> A_unitig配列) Get_グラフ(int p_k長, params string[] p_unitig群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = p_k長, A_スレッド数 = 1 };
            List<string> l_unitig配列 = [string.Empty, string.Empty];
            Dictionary<KmerKey, (int A_unitigID, int A_位置)> l_kmer辞書 = [];

            var l_ID = 1;
            foreach (var l_配列 in p_unitig群)
            {
                l_unitig配列.Add(l_配列);
                l_unitig配列.Add(Util.V_逆相補(l_配列));
                for (var i = p_k長; i <= l_配列.Length; i++)
                {
                    var l_キー = new KmerKey(l_配列.AsSpan(i - p_k長, p_k長));
                    V_登録(l_kmer辞書, l_キー, l_ID, i - p_k長);
                    V_登録(l_kmer辞書, l_キー.Get_逆相補(), -l_ID, l_配列.Length - i);
                }
                l_ID++;
            }
            return (UnitigGraph.Get_グラフ(l_unitig配列, l_kmer辞書, p_k長, 曖昧kmer番号), l_unitig配列);
        }

        /// <summary>
        /// 頂点の唯一の出辺の行き先
        /// </summary>
        /// <param name="p_グラフ"></param>
        /// <param name="p_頂点"></param>
        /// <returns></returns>
        public static int Get_唯一の次(UnitigGraph p_グラフ, int p_頂点)
        {
            return Assert.Single(p_グラフ.A_出辺[p_頂点]);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// k-mer を、それが載る unitig と開始位置の辞書へ登録する
        /// </summary>
        /// <param name="p_辞書"></param>
        /// <param name="p_キー"></param>
        /// <param name="p_ID"></param>
        /// <param name="p_位置"></param>
        private static void V_登録(Dictionary<KmerKey, (int A_unitigID, int A_位置)> p_辞書, KmerKey p_キー, int p_ID, int p_位置)
        {
            if (p_辞書.TryGetValue(p_キー, out var l_既存))
            {
                if (l_既存.A_unitigID is not 曖昧kmer番号 && l_既存.A_unitigID != p_ID)
                {
                    p_辞書[p_キー] = (曖昧kmer番号, 0);
                }
                return;
            }
            p_辞書[p_キー] = (p_ID, p_位置);
        }

        #endregion
    }
}
