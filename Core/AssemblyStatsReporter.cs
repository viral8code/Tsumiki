using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    internal static class AssemblyStatsReporter
    {
        /// <summary>
        /// 他アセンブラとの比較で慣習的に使われる最小長。ABySS 付属の
        /// abyss-fac が既定でこの長さ以上の配列だけを集計するため、
        /// 公表されている N50 等はほぼこの条件で計算されている。
        /// 全件の統計だけを出していると、数十bpの断片が大量に混じった
        /// こちらの数字と比較して不当に悪く見える(あるいはその逆になる)。
        /// </summary>
        public const int 比較用の最小長 = 500;

        public static アセンブリ統計 Get_統計(IEnumerable<string> p_配列群)
        {
            var l_長さ一覧 = new List<int>();
            long l_総延長 = 0;
            long l_GC数 = 0;
            long l_塩基数 = 0;

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
                return new アセンブリ統計(0, 0, 0, 0, 0, 0, 0);
            }

            l_長さ一覧.Sort();
            l_長さ一覧.Reverse();

            var l_半分 = l_総延長 / 2.0;
            long l_累積 = 0;
            var l_N50 = l_長さ一覧[^1];
            var l_L50 = l_長さ一覧.Count;
            for (var i = 0; i < l_長さ一覧.Count; i++)
            {
                l_累積 += l_長さ一覧[i];
                if (l_累積 >= l_半分)
                {
                    l_N50 = l_長さ一覧[i];
                    l_L50 = i + 1;
                    break;
                }
            }

            var l_GC率 = l_塩基数 == 0 ? 0.0 : (100.0 * l_GC数 / l_塩基数);

            return new アセンブリ統計(
                p_配列数: l_長さ一覧.Count,
                p_総延長: l_総延長,
                p_最大長: l_長さ一覧[0],
                p_最小長: l_長さ一覧[^1],
                p_N50: l_N50,
                p_L50: l_L50,
                p_GC率: l_GC率);
        }

        public static アセンブリ統計 Get_統計_FASTA(string p_FASTAパス)
        {
            return Get_統計(Get_配列群(p_FASTAパス));
        }

        private static IEnumerable<string> Get_配列群(string p_FASTAパス)
        {
            using var l_読み込み = new FastaReader(p_FASTAパス);
            while (l_読み込み.Get_続きがあるか())
            {
                yield return l_読み込み.Get_次の配列().A_配列;
            }
        }

        /// <summary>
        /// FASTA の統計量を計算し、"[Stats] ラベル: ..." の形式でコンソールへ出力する。
        /// 全配列を対象とした統計に加えて、他アセンブラの公表値と直接比較できるよう
        /// 比較用の最小長 以上の配列だけに絞った統計も併記する。
        /// </summary>
        public static void V_出力_統計(string p_ラベル, string p_FASTAパス)
        {
            if (!File.Exists(p_FASTAパス))
            {
                Console.WriteLine($"[Stats] {p_ラベル}: (file not found: {p_FASTAパス})");
                return;
            }
            var l_統計 = Get_統計_FASTA(p_FASTAパス);
            Console.WriteLine($"[Stats] {p_ラベル}: {l_統計}");

            var l_絞り込み統計 = Get_統計(Get_配列群(p_FASTAパス).Where(x => x.Length >= 比較用の最小長));
            Console.WriteLine($"[Stats] {p_ラベル} (>={比較用の最小長}bp, comparable to abyss-fac): {l_絞り込み統計}");
        }
    }
}
