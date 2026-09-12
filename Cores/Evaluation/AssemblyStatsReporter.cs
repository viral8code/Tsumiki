using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// アセンブリの長さの統計を求めて出力する
    /// </summary>
    internal static class AssemblyStatsReporter
    {
        #region 定数

        /// <summary>
        /// 他アセンブラとの比較で慣習的に使われる最小長 (abyss-fac の既定)
        /// </summary>
        /// <remarks>
        /// 全件の統計だけでは、短い断片を含むぶん公表値と比較にならない
        /// </remarks>
        public const int 比較用の最小長 = 500;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列群から統計を求めて返す
        /// </summary>
        /// <param name="p_配列群">対象の配列</param>
        /// <returns>アセンブリ統計</returns>
        public static アセンブリ統計 Get_統計(IEnumerable<string> p_配列群)
        {
            var l_長さ一覧 = new List<int>();
            var l_総延長 = 0L;
            var l_GC数 = 0L;
            var l_塩基数 = 0L;

            foreach (var l_配列 in p_配列群)
            {
                l_長さ一覧.Add(l_配列.Length);
                l_総延長 += l_配列.Length;
                foreach (var l_文字 in l_配列)
                {
                    if (l_文字 is 'N' or 'n')
                    {
                        continue;
                    }

                    l_塩基数++;
                    if (l_文字 is 'G' or 'g' or 'C' or 'c')
                    {
                        l_GC数++;
                    }
                }
            }

            if (l_長さ一覧.Count == 0)
            {
                return new アセンブリ統計(0, 0L, 0, 0, 0, 0, 0D);
            }

            var (l_N50, l_L50) = StatsUtil.Get_N50([.. l_長さ一覧.Select(x => (long)x)]);
            var l_GC率 = l_塩基数 == 0D ? 0D : (100D * l_GC数 / l_塩基数);

            return new アセンブリ統計(A_配列数: l_長さ一覧.Count, A_総延長: l_総延長, A_最大長: l_長さ一覧.Max(), A_最小長: l_長さ一覧.Min(), A_N50: (int)l_N50, A_L50: l_L50, A_GC率: l_GC率);
        }

        /// <summary>
        /// FASTA から統計を求めて返す
        /// </summary>
        /// <param name="p_FASTAパス">対象の FASTA のパス</param>
        /// <returns>アセンブリ統計</returns>
        public static アセンブリ統計 Get_統計_FASTA(string p_FASTAパス)
        {
            return Get_統計(Get_配列群(p_FASTAパス));
        }

        /// <summary>
        /// FASTA の統計量を計算し、"[Stats] ラベル: ..." の形式でコンソールへ出力する
        /// </summary>
        /// <param name="p_ラベル">出力の見出しに使うラベル</param>
        /// <param name="p_FASTAパス">対象の FASTA のパス</param>
        /// <remarks>
        /// 全配列を対象とした統計に加えて、他アセンブラの公表値と直接比較できるよう比較用の最小長以上の配列だけに絞った統計も併記する
        /// </remarks>
        public static void V_出力_統計(string p_ラベル, string p_FASTAパス)
        {
            if (!File.Exists(p_FASTAパス))
            {
                Logger.V_出力(メッセージID.統計_ファイルなし, p_ラベル, p_FASTAパス);
                return;
            }

            var l_統計 = Get_統計_FASTA(p_FASTAパス);
            Logger.V_出力(メッセージID.統計, p_ラベル, l_統計);

            var l_絞り込み統計 = Get_統計(Get_配列群(p_FASTAパス).Where(x => x.Length >= 比較用の最小長));
            Logger.V_出力(メッセージID.統計_長さで絞り込み, p_ラベル, 比較用の最小長, l_絞り込み統計);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// FASTA から配列だけを取り出して返す
        /// </summary>
        /// <param name="p_FASTAパス">対象の FASTA のパス</param>
        /// <returns>配列</returns>
        private static IEnumerable<string> Get_配列群(string p_FASTAパス)
        {
            using var l_読み込み = new FastaReader(p_FASTAパス);
            while (l_読み込み.Has続き())
            {
                yield return l_読み込み.Get_次の配列().A_配列;
            }
        }

        #endregion
    }
}
