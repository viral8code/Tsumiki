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
        public const int C_比較用の最小長 = 500;

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
        /// N 連続区間で分割した contig の統計を求めて返す
        /// </summary>
        /// <param name="p_配列群"></param>
        /// <param name="p_最小長"></param>
        /// <returns></returns>
        public static アセンブリ統計 Get_N分割統計(IEnumerable<string> p_配列群, int p_最小長 = 0)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(p_最小長);

            List<string> l_断片群 = [];
            foreach (var l_配列 in p_配列群)
            {
                l_断片群.AddRange(l_配列.Split(['N', 'n'], StringSplitOptions.RemoveEmptyEntries).Where(x => x.Length >= p_最小長));
            }

            return Get_統計(l_断片群);
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
        public static void V_出力_統計(string p_ラベル, string p_FASTAパス)
        {
            if (!File.Exists(p_FASTAパス))
            {
                Logger.V_出力(メッセージID.統計_ファイルなし, p_ラベル, p_FASTAパス);
                return;
            }

            var l_統計 = Get_統計_FASTA(p_FASTAパス);
            Logger.V_出力(メッセージID.統計, p_ラベル, l_統計);

            var l_絞り込み統計 = Get_統計(Get_配列群(p_FASTAパス).Where(x => x.Length >= C_比較用の最小長));
            Logger.V_出力(メッセージID.統計_長さで絞り込み, p_ラベル, C_比較用の最小長, l_絞り込み統計);

            var l_N分割統計 = Get_N分割統計(Get_配列群(p_FASTAパス), C_比較用の最小長);
            Logger.V_出力_そのまま($"[Stats] {p_ラベル} (N-split, >= {C_比較用の最小長}bp): {l_N分割統計}");
        }

        /// <summary>
        /// 複数の FASTA の統計を 1 つの Markdown の表の行にする
        /// </summary>
        /// <param name="p_対象群">ラベルと FASTA のパス、存在しないパスは飛ばす</param>
        /// <returns>見出し行を含む表の行</returns>
        public static List<string> Get_統計表(IReadOnlyList<(string A_ラベル, string A_FASTAパス)> p_対象群)
        {
            List<string> l_行群 =
            [
                "| file | filter | count | total length | max | min | N50 | L50 | GC% |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|",
            ];
            foreach (var (l_ラベル, l_パス) in p_対象群)
            {
                if (!File.Exists(l_パス))
                {
                    continue;
                }

                var l_配列群 = Get_配列群(l_パス).ToList();
                (string A_条件, アセンブリ統計 A_統計)[] l_絞り込み群 =
                [
                    ("all", Get_統計(l_配列群)),
                    ($">= {C_比較用の最小長}bp", Get_統計(l_配列群.Where(x => x.Length >= C_比較用の最小長))),
                    ($"N-split, >= {C_比較用の最小長}bp", Get_N分割統計(l_配列群, C_比較用の最小長)),
                ];
                foreach (var (l_条件, l_統計) in l_絞り込み群)
                {
                    l_行群.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"| {l_ラベル} | {l_条件} | {l_統計.A_配列数:N0} | {l_統計.A_総延長:N0} | {l_統計.A_最大長:N0} | {l_統計.A_最小長:N0} | {l_統計.A_N50:N0} | {l_統計.A_L50:N0} | {l_統計.A_GC率:F2} |"));
                }
            }

            return l_行群;
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
