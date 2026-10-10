using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 選んだアセンブリに無い配列を、他の k のアセンブリから補う
    /// </summary>
    internal static class AssemblyComplementer
    {
        #region 定数

        /// <summary>
        /// 補う配列の最小長 (評価に含める配列の最小長と揃える)
        /// </summary>
        public const int C_補う最小長 = 500;

        /// <summary>
        /// 補う配列のアンカー k-mer のうち、リードで信頼できるものに求める割合
        /// </summary>
        public const double C_信頼kmerの最小割合 = 0.9D;

        /// <summary>
        /// 欠けた区間どうしをつなぐときに、間に挟んでよい既知の区間の最長 (短い反復で区間が割れるのを防ぐ)
        /// </summary>
        public const int C_挟める既知区間の最長 = 500;

        /// <summary>
        /// つないだ区間のうち、欠けた区間が占めるべき割合の下限
        /// </summary>
        public const double C_欠けの最小割合 = 0.5D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 骨格に無いアンカー k-mer が続く区間 (短い既知の区間を挟んでつないだもの) を他の候補から拾い、骨格の後ろに足して p_出力パス へ書き出す
        /// </summary>
        /// <param name="p_骨格">選んだアセンブリ</param>
        /// <param name="p_全候補">全 k の実行結果 (骨格と同じ k は使わない)</param>
        /// <param name="p_Is信頼">アンカー k-mer の正規形がリードで信頼できるか</param>
        /// <param name="p_アンカーk長">アンカー k 長</param>
        /// <param name="p_出力パス">補った結果の書き出し先 (補う配列が無ければ書き出さない)</param>
        /// <returns>補った配列の本数と延長</returns>
        public static (int A_本数, long A_延長) Get_補完(アセンブリ実行結果 p_骨格, IReadOnlyList<アセンブリ実行結果> p_全候補, Func<UInt128, bool> p_Is信頼, int p_アンカーk長, string p_出力パス)
        {
            var l_骨格 = FastaReader.Get_全エントリ(p_骨格.A_最終パス);
            var l_既知 = new HashSet<UInt128>();
            foreach (var (_, l_配列) in l_骨格)
            {
                V_登録_kmer(l_既知, l_配列, p_アンカーk長);
            }

            List<(string A_ID, string A_配列)> l_補った = [];
            foreach (var l_他 in p_全候補.Where(x => x.A_k長 != p_骨格.A_k長).OrderByDescending(x => x.A_k長))
            {
                var l_連番 = 1;
                foreach (var (_, l_配列) in FastaReader.Get_全エントリ(l_他.A_最終パス))
                {
                    foreach (var (l_開始, l_長さ) in Get_つないだ区間(l_配列, [.. Get_欠けた区間(l_配列, l_既知, p_アンカーk長)]))
                    {
                        var l_区間 = l_配列.Substring(l_開始, l_長さ);
                        if (l_長さ < C_補う最小長 || Get_信頼kmerの割合(l_区間, p_Is信頼, p_アンカーk長) < C_信頼kmerの最小割合)
                        {
                            continue;
                        }

                        l_補った.Add(($"{Consts.補完配列の接頭辞}{l_他.A_k長}_{l_連番++}", l_区間));
                        V_登録_kmer(l_既知, l_区間, p_アンカーk長);
                    }
                }
            }

            if (l_補った.Count == 0)
            {
                return (0, 0L);
            }

            using (var l_書き込み = new FastaWriter(p_出力パス))
            {
                foreach (var (l_ID, l_配列) in l_骨格.Concat(l_補った))
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }

            return (l_補った.Count, l_補った.Sum(x => (long)x.A_配列.Length));
        }

        /// <summary>
        /// 既知の集合に無い k-mer が続く区間を、配列上の開始位置と長さで列挙する
        /// </summary>
        /// <param name="p_配列">調べる配列</param>
        /// <param name="p_既知">既知の k-mer の正規形</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>区間の開始位置と長さ</returns>
        public static IEnumerable<(int A_開始, int A_長さ)> Get_欠けた区間(string p_配列, IReadOnlySet<UInt128> p_既知, int p_k長)
        {
            var l_インデックス1 = 0;
            while (l_インデックス1 + p_k長 <= p_配列.Length)
            {
                if (!Is未知(p_配列, l_インデックス1, p_既知, p_k長))
                {
                    l_インデックス1++;
                    continue;
                }

                var l_インデックス2 = l_インデックス1;
                while (l_インデックス2 + p_k長 <= p_配列.Length && Is未知(p_配列, l_インデックス2, p_既知, p_k長))
                {
                    l_インデックス2++;
                }

                yield return (l_インデックス1, l_インデックス2 - 1 + p_k長 - l_インデックス1);
                l_インデックス1 = l_インデックス2;
            }
        }

        /// <summary>
        /// 間の既知の区間が短く塩基だけから成る欠けた区間どうしをつなぎ、欠けが十分に占めるものだけを返す
        /// </summary>
        /// <param name="p_配列">区間を取った配列</param>
        /// <param name="p_欠けた区間">配列上の順に並んだ欠けた区間</param>
        /// <returns>つないだ区間の開始位置と長さ</returns>
        public static IEnumerable<(int A_開始, int A_長さ)> Get_つないだ区間(string p_配列, IReadOnlyList<(int A_開始, int A_長さ)> p_欠けた区間)
        {
            var l_インデックス1 = 0;
            while (l_インデックス1 < p_欠けた区間.Count)
            {
                var l_開始 = p_欠けた区間[l_インデックス1].A_開始;
                var l_終了 = l_開始 + p_欠けた区間[l_インデックス1].A_長さ;
                var l_欠け = p_欠けた区間[l_インデックス1].A_長さ;
                var l_インデックス2 = l_インデックス1 + 1;
                while (l_インデックス2 < p_欠けた区間.Count && p_欠けた区間[l_インデックス2].A_開始 - l_終了 <= C_挟める既知区間の最長 && Is塩基のみ(p_配列, l_終了, p_欠けた区間[l_インデックス2].A_開始))
                {
                    l_終了 = Math.Max(l_終了, p_欠けた区間[l_インデックス2].A_開始 + p_欠けた区間[l_インデックス2].A_長さ);
                    l_欠け += p_欠けた区間[l_インデックス2].A_長さ;
                    l_インデックス2++;
                }

                if (l_欠け >= (l_終了 - l_開始) * C_欠けの最小割合)
                {
                    yield return (l_開始, l_終了 - l_開始);
                }

                l_インデックス1 = l_インデックス2;
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 位置 p_位置 の k-mer が塩基だけから成り、既知の集合に無いか
        /// </summary>
        /// <param name="p_配列">調べる配列</param>
        /// <param name="p_位置">k-mer の開始位置</param>
        /// <param name="p_既知">既知の k-mer の正規形</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>未知なら true</returns>
        private static bool Is未知(string p_配列, int p_位置, IReadOnlySet<UInt128> p_既知, int p_k長)
        {
            return KmerPacking.Is成功_正規化パック(p_配列, p_位置, p_k長, out var l_鍵) && !p_既知.Contains(l_鍵);
        }

        /// <summary>
        /// 配列の [p_開始, p_終了) が塩基だけから成るか (重なっていれば true)
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_開始">開始位置</param>
        /// <param name="p_終了">終了位置 (含まない)</param>
        /// <returns>塩基だけなら true</returns>
        private static bool Is塩基のみ(string p_配列, int p_開始, int p_終了)
        {
            return p_終了 <= p_開始 || p_配列.AsSpan(p_開始, p_終了 - p_開始).IndexOfAnyExcept(Consts.塩基文字) < 0;
        }

        /// <summary>
        /// 配列の k-mer の正規形を集合へ足す
        /// </summary>
        /// <param name="p_集合">足す先</param>
        /// <param name="p_配列">配列</param>
        /// <param name="p_k長">k 長</param>
        private static void V_登録_kmer(HashSet<UInt128> p_集合, string p_配列, int p_k長)
        {
            for (var i = 0; i + p_k長 <= p_配列.Length; i++)
            {
                if (KmerPacking.Is成功_正規化パック(p_配列, i, p_k長, out var l_鍵))
                {
                    _ = p_集合.Add(l_鍵);
                }
            }
        }

        /// <summary>
        /// 配列のアンカー k-mer のうち、リードで信頼できるものの割合
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_Is信頼">アンカー k-mer の正規形がリードで信頼できるか</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>割合 (k-mer が無ければ 0)</returns>
        private static double Get_信頼kmerの割合(string p_配列, Func<UInt128, bool> p_Is信頼, int p_k長)
        {
            var l_件数 = 0;
            var l_信頼 = 0;
            for (var i = 0; i + p_k長 <= p_配列.Length; i++)
            {
                if (KmerPacking.Is成功_正規化パック(p_配列, i, p_k長, out var l_鍵))
                {
                    l_件数++;
                    l_信頼 += p_Is信頼(l_鍵) ? 1 : 0;
                }
            }

            return l_件数 == 0 ? 0D : (double)l_信頼 / l_件数;
        }

        #endregion
    }
}
