using Tsumiki.Common;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// 各 unitig がゲノム中に何回現れるか(コピー数)を、カバレッジから推定する。
    ///
    /// 考え方: ゲノム中に1回しか現れない領域のカバレッジを基準値とすると、
    /// n 回現れる反復配列にはリードが n 倍集まるので、カバレッジも約 n 倍になる。
    /// したがって「その unitig の平均カバレッジ / 基準値」を丸めればコピー数になる。
    ///
    /// なぜ必要か:
    /// - 反復配列かどうかを、グラフの形(入次数・出次数)ではなく量的な根拠で判定できる。
    ///   入次数2・出次数2でも実は単一コピー(バブルの残骸)ということがありうるし、
    ///   逆に次数が1でも高カバレッジならタンデムリピートを1本に潰している疑いがある。
    /// - ゲノム全体の経路探索で「この unitig は2回使ってよい」という予算になる。
    ///   予算が無いと、探索は同じ反復を何度でも通れてしまい破綻する。
    ///
    /// 基準値は全 unitig の平均カバレッジの「長さ加重中央値」を使う。
    /// 単純平均・単純中央値だと、本数として多い短い断片(エラー由来の残骸など)に
    /// 引きずられる。塩基数で重み付けすれば、ゲノムの大部分を占める単一コピー領域の
    /// 水準が選ばれる。
    /// </summary>
    internal static class CopyNumberEstimator
    {
        /// <summary>
        /// これを下回るカバレッジ比の unitig は、コピー数を推定できるだけの
        /// 根拠が無いとみなして 1 として扱う(0 コピーにはしない)。
        /// </summary>
        private const double 多コピーとみなす比の下限 = 1.5;

        /// <summary>
        /// コピー数の上限。これを超える比が出た場合、rRNA オペロンのような
        /// 高コピー反復か、あるいはカバレッジ異常のどちらかで区別がつかない。
        /// 経路探索の予算としては大きすぎると探索が発散するため頭打ちにする。
        /// </summary>
        private const int コピー数の上限 = 12;

        /// <summary>
        /// unitig ID(1始まり) -> その unitig を構成する k-mer の平均カバレッジ、を計算する。
        /// </summary>
        public static Dictionary<int, double> Get_カバレッジ(
            TrustedKmerIndex p_kmerインデックス,
            IReadOnlyDictionary<int, string> p_ユニティグ配列,
            int p_k長)
        {
            Dictionary<int, double> l_カバレッジ = [];
            foreach (var (l_ID, l_配列) in p_ユニティグ配列)
            {
                if (l_配列.Length < p_k長)
                {
                    l_カバレッジ[l_ID] = 0;
                    continue;
                }

                var l_塩基列 = new byte[l_配列.Length];
                for (var i = 0; i < l_配列.Length; i++)
                {
                    l_塩基列[i] = Util.Get_塩基ID(l_配列[i]);
                }

                ulong l_合計 = 0;
                var l_件数 = 0;
                for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
                {
                    l_合計 += p_kmerインデックス.Get_カバレッジ(l_塩基列.AsSpan(i, p_k長));
                    l_件数++;
                }
                l_カバレッジ[l_ID] = l_件数 == 0 ? 0 : (double)l_合計 / l_件数;
            }
            return l_カバレッジ;
        }

        /// <summary>
        /// カバレッジからコピー数を推定する。
        /// </summary>
        public static コピー数推定結果 Get_推定結果(
            IReadOnlyDictionary<int, double> p_カバレッジ,
            IReadOnlyDictionary<int, int> p_ユニティグ長)
        {
            var l_基準値 = Get_長さ加重中央値(p_カバレッジ, p_ユニティグ長);

            Dictionary<int, int> l_コピー数 = [];
            foreach (var (l_ID, l_カバレッジ値) in p_カバレッジ)
            {
                if (l_基準値 <= 0)
                {
                    l_コピー数[l_ID] = 1;
                    continue;
                }

                var l_比 = l_カバレッジ値 / l_基準値;
                if (l_比 < 多コピーとみなす比の下限)
                {
                    // 単一コピー(あるいは低カバレッジで判断できない)。
                    // 0 にはしない: 実際に配列は存在しており、経路から
                    // 締め出してしまうと組み立てられなくなる。
                    l_コピー数[l_ID] = 1;
                    continue;
                }

                l_コピー数[l_ID] = Math.Clamp((int)Math.Round(l_比), 1, コピー数の上限);
            }

            return new コピー数推定結果(l_基準値, p_カバレッジ, l_コピー数);
        }

        /// <summary>
        /// 長さで重み付けしたカバレッジの中央値。ゲノムの大部分を占める
        /// 単一コピー領域の水準を推定するために使う。
        /// </summary>
        private static double Get_長さ加重中央値(
            IReadOnlyDictionary<int, double> p_カバレッジ,
            IReadOnlyDictionary<int, int> p_ユニティグ長)
        {
            var l_組 = p_カバレッジ
                .Where(x => p_ユニティグ長.ContainsKey(x.Key) && x.Value > 0)
                .Select(x => (A_長さ: (long)p_ユニティグ長[x.Key], A_カバレッジ: x.Value))
                .OrderBy(x => x.A_カバレッジ)
                .ToList();
            if (l_組.Count == 0)
            {
                return 0;
            }

            var l_総延長 = l_組.Sum(x => x.A_長さ);
            if (l_総延長 == 0)
            {
                return 0;
            }

            var l_半分 = l_総延長 / 2.0;
            long l_累積 = 0;
            foreach (var (l_長さ, l_カバレッジ値) in l_組)
            {
                l_累積 += l_長さ;
                if (l_累積 >= l_半分)
                {
                    return l_カバレッジ値;
                }
            }
            return l_組[^1].A_カバレッジ;
        }

        /// <summary>
        /// 推定結果の要約をコンソールへ出力する。
        /// 「単一コピーが何本・何bp、2コピー以上が何本・何bp」が分かると、
        /// 反復配列がアセンブリのどれだけを占めているかが把握できる。
        /// </summary>
        public static void V_出力_推定結果(コピー数推定結果 p_推定結果, IReadOnlyDictionary<int, int> p_ユニティグ長)
        {
            Console.WriteLine($"[Info] Single-copy coverage baseline estimated as {p_推定結果.A_単一コピー基準値:0.#} (length-weighted median).");

            var l_コピー数別 = p_推定結果.A_コピー数
                .GroupBy(x => x.Value)
                .OrderBy(x => x.Key)
                .Select(x => (A_コピー数: x.Key,
                              A_本数: x.Count(),
                              A_塩基数: x.Sum(y => (long)p_ユニティグ長.GetValueOrDefault(y.Key, 0))))
                .ToList();

            var l_要約 = string.Join(", ", l_コピー数別.Select(x => $"x{x.A_コピー数}: {x.A_本数} unitig(s)/{x.A_塩基数:N0}bp"));
            Console.WriteLine($"[Info] Estimated copy numbers -- {l_要約}");

            var l_反復塩基数 = l_コピー数別.Where(x => x.A_コピー数 >= 2).Sum(x => x.A_塩基数);
            var l_総塩基数 = l_コピー数別.Sum(x => x.A_塩基数);
            if (l_総塩基数 > 0)
            {
                Console.WriteLine(
                    $"[Info] Multi-copy (repeat) content: {l_反復塩基数:N0}bp of {l_総塩基数:N0}bp " +
                    $"({100.0 * l_反復塩基数 / l_総塩基数:0.0}% of the assembly is sequence that occurs more than once).");
            }
        }
    }
}
