using Tsumiki.Commons;
using Tsumiki.Cores.Output;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 最終配列の継ぎ目 (反復の区間を挟んで固有配列を繋いだ箇所と N のギャップ) ごとに、繋ぐ相手を誤っている確率を出す
    /// </summary>
    internal static class JunctionRiskEvaluator
    {
        #region 定数

        /// <summary>
        /// 項目 junction risk
        /// </summary>
        private const string C_項目_junction_risk = "junction-risk";

        /// <summary>
        /// 項目
        /// </summary>
        private const string C_項目 = "_";

        /// <summary>
        /// 配列の中の反復を見る k-mer の長さ
        /// </summary>
        private const int C_反復のk長 = 31;

        /// <summary>
        /// 候補の両側に要る固有配列の長さ
        /// </summary>
        private const int C_固有の最小長 = 100;

        /// <summary>
        /// 候補とする反復の区間の最小長
        /// </summary>
        private const int C_反復の最小長 = 20;

        /// <summary>
        /// 反復の区間を 1 つにつなぐ隙間の上限
        /// </summary>
        private const int C_つなぐ隙間 = 50;

        /// <summary>
        /// 跨いだとみなすのに区間の外へ踏み込む長さ
        /// </summary>
        private const int C_端の許容 = 10;

        /// <summary>
        /// 畳まれた反復とみなす、深さの中央値に対する比
        /// </summary>
        private const double C_畳まれた反復の深さ比 = 1.8D;

        /// <summary>
        /// 畳まれた反復とみなす、一意に当たったリードの割合の上限
        /// </summary>
        private const double C_一意な割合の下限 = 0.5D;

        /// <summary>
        /// 深さの中央値を取る間隔
        /// </summary>
        private const int C_深さを取る間隔 = 50;

        /// <summary>
        /// 両側の深さを取る長さ
        /// </summary>
        private const int C_両側の深さの窓 = 100;

        /// <summary>
        /// 切れ端を数える、区間の端からの距離
        /// </summary>
        private const int C_切れ端の窓 = 5;

        /// <summary>
        /// 跨ぐ読みを探しに戻る距離
        /// </summary>
        private const int C_読みを探す戻り幅 = 1_000;

        /// <summary>
        /// 断片長の分布から取る標本の、配列あたりの目安
        /// </summary>
        private const int C_配列あたりの断片の標本数 = 2_000;

        /// <summary>
        /// 期待の組に掛ける、両側の深さ比の上限
        /// </summary>
        private const double C_深さ比の上限 = 2D;

        /// <summary>
        /// 評価に要る、組めたペアの数
        /// </summary>
        private const int C_必要な断片数 = 100;

        /// <summary>
        /// 安全版で切る継ぎ目の、誤りの確率の下限 (誤りの方が起きやすい継ぎ目だけを切る)
        /// </summary>
        internal const double C_切る確率 = 0.5D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_FASTAパス の継ぎ目を評価して p_出力パス へ書き出す
        /// </summary>
        /// <param name="p_FASTAパス">最終配列</param>
        /// <param name="p_ライブラリ群">ライブラリごとのリードの組</param>
        /// <param name="p_出力パス">書き出し先の TSV</param>
        /// <returns>評価、ペアのリードが足りず評価できなければ null</returns>
        public static IReadOnlyList<継ぎ目の評価>? Get_評価結果(string p_FASTAパス, IReadOnlyList<(string A_順リード, string A_逆リード)> p_ライブラリ群, string p_出力パス)
        {
            using var l_計測 = new StageTimer(C_項目_junction_risk);
            Logger.V_出力(メッセージID.継ぎ目の評価開始);
            if (Get_評価(p_FASTAパス, p_ライブラリ群) is not { } l_評価群)
            {
                Logger.V_出力(メッセージID.継ぎ目の評価を行えず);
                return null;
            }

            ReportWriter.V_書き出し_継ぎ目(p_出力パス, l_評価群);
            Logger.V_出力(メッセージID.継ぎ目の評価結果, l_評価群.Count, l_評価群.Sum(x => x.A_誤りの確率), l_評価群.Count(x => x.A_誤りの確率 >= C_切る確率), C_切る確率, p_出力パス);
            return l_評価群;
        }

        /// <summary>
        /// p_FASTAパス の継ぎ目ごとに、繋ぐ相手を誤っている確率を出す (書き出さない)
        /// </summary>
        /// <param name="p_FASTAパス">評価する配列</param>
        /// <param name="p_ライブラリ群">ライブラリごとのリードの組</param>
        /// <returns>評価、ペアのリードが足りず評価できなければ null</returns>
        public static List<継ぎ目の評価>? Get_評価(string p_FASTAパス, IReadOnlyList<(string A_順リード, string A_逆リード)> p_ライブラリ群)
        {
            var l_エントリ群 = FastaReader.Get_全エントリ(p_FASTAパス);
            if (l_エントリ群.Count == 0 || p_ライブラリ群.All(x => string.IsNullOrWhiteSpace(x.A_逆リード)))
            {
                return null;
            }

            var l_配列群 = l_エントリ群.Select(x => x.A_配列.ToUpperInvariant()).ToList();
            var l_証拠 = JunctionEvidenceCollector.Get_証拠(l_配列群, p_ライブラリ群);
            if (l_証拠.A_断片.Sum(x => x.Length) < C_必要な断片数)
            {
                return null;
            }

            var l_名前群 = l_エントリ群.Select(x => x.A_ID.TrimStart('>').Split(' ', 2)[0]).ToList();
            return [.. Get_候補群(l_証拠, l_名前群).Select(JunctionRiskModel.Get_評価)];
        }

        /// <summary>
        /// 誤りの確率が 切る確率 以上の継ぎ目で配列を切り、p_出力パス へ書き出す
        /// </summary>
        /// <param name="p_FASTAパス">最終配列</param>
        /// <param name="p_評価群">継ぎ目の評価</param>
        /// <param name="p_出力パス">書き出し先の FASTA</param>
        public static void V_書き出し_切った配列(string p_FASTAパス, IReadOnlyList<継ぎ目の評価> p_評価群, string p_出力パス)
        {
            var l_切る所 = p_評価群.Where(x => x.A_誤りの確率 >= C_切る確率).GroupBy(x => x.A_候補.A_配列名).ToDictionary(x => x.Key, x => x.Select(y => (y.A_候補.A_開始, y.A_候補.A_終了)).ToList());
            var l_本数 = 0;
            using (var l_書き込み = new FastaWriter(p_出力パス))
            {
                foreach (var (l_ID, l_配列) in FastaReader.Get_全エントリ(p_FASTAパス))
                {
                    var l_名前 = l_ID.TrimStart('>').Split(' ', 2)[0];
                    foreach (var (l_片の名前, l_片) in Get_切った配列群(l_名前, l_配列, l_切る所.GetValueOrDefault(l_名前) ?? []))
                    {
                        l_書き込み.V_書き込み(l_片の名前, l_片);
                        l_本数++;
                    }
                }
            }

            Logger.V_出力(メッセージID.継ぎ目で切った配列を書き出した, l_切る所.Sum(x => x.Value.Count), C_切る確率, l_本数, p_出力パス);
        }

        /// <summary>
        /// 配列を切る所 (反復の区間) で切る<br/>
        /// 反復の区間は両側の片に残し (どちらのコピーかは分からないため) 、区間に N のギャップがあればその手前と後ろで切る<br/>
        /// 切ったら環状の目印は外す
        /// </summary>
        /// <param name="p_名前">配列名</param>
        /// <param name="p_配列"></param>
        /// <param name="p_切る所">反復の区間 (始まり, 終わり)</param>
        /// <returns>(片の名前, 片)</returns>
        public static List<(string A_名前, string A_配列)> Get_切った配列群(string p_名前, string p_配列, IReadOnlyList<(int A_開始, int A_終了)> p_切る所)
        {
            if (p_切る所.Count == 0)
            {
                return [(p_名前, p_配列)];
            }

            var l_名前 = p_名前.Replace(C_項目 + Consts.環状の目印, string.Empty, StringComparison.OrdinalIgnoreCase);
            List<(string A_名前, string A_配列)> l_片群 = [];
            var l_始まり = 0;
            foreach (var (l_開始, l_終了) in p_切る所.Order())
            {
                var l_区間 = p_配列.AsSpan(l_開始, l_終了 - l_開始);
                var l_最初のN = l_区間.IndexOfAny('N', 'n');
                if (l_最初のN < 0)
                {
                    V_追加_片(l_片群, l_名前, p_配列[l_始まり..l_終了]);
                    l_始まり = l_開始;
                    continue;
                }

                V_追加_片(l_片群, l_名前, p_配列[l_始まり..(l_開始 + l_最初のN)]);
                l_始まり = l_開始 + l_区間.LastIndexOfAny('N', 'n') + 1;
            }

            V_追加_片(l_片群, l_名前, p_配列[l_始まり..]);
            return l_片群;
        }

        /// <summary>
        /// 証拠から継ぎ目の候補を拾い、特徴量を付ける
        /// </summary>
        /// <param name="p_証拠"></param>
        /// <param name="p_名前群">配列名</param>
        /// <returns>候補</returns>
        public static List<継ぎ目候補> Get_候補群(継ぎ目の証拠 p_証拠, IReadOnlyList<string> p_名前群)
        {
            var l_配列群 = p_証拠.A_配列群;
            var l_反復の印 = Get_反復の印(l_配列群);
            var l_中央 = Get_固有の深さの中央値(p_証拠, l_反復の印);
            var l_断片長 = p_証拠.A_断片.SelectMany(x => x).Select(x => x.A_終了 - x.A_開始).Order().ToArray();
            var l_p1 = l_断片長[(int)(0.01D * (l_断片長.Length - 1))];
            var l_p50 = l_断片長[(int)(0.5D * (l_断片長.Length - 1))];
            var l_p99 = l_断片長[(int)(0.99D * (l_断片長.Length - 1))];
            var l_断片密度 = (double)l_断片長.Count(x => x >= l_p1 && x <= l_p99) / l_配列群.Sum(x => (long)x.Length);
            var l_断片長の標本 = p_証拠.A_断片.SelectMany(x => x.Where((_, i) => i % Math.Max(1, x.Length / C_配列あたりの断片の標本数) == 0)).Select(x => x.A_終了 - x.A_開始).Where(x => x >= l_p1 && x <= l_p99).ToArray();
            List<継ぎ目候補> l_候補群 = [];
            for (var c = 0; c < l_配列群.Count; c++)
            {
                var l_配列 = l_配列群[c];
                var l_長さ = l_配列.Length;
                foreach (var (l_a, l_b) in Get_反復の区間(l_反復の印[c], p_証拠.A_深さ[c], p_証拠.A_一意な深さ[c], l_中央))
                {
                    if (l_b - l_a < C_反復の最小長 || l_a < C_固有の最小長 || l_長さ - l_b < C_固有の最小長)
                    {
                        continue;
                    }

                    var l_N数 = l_配列.AsSpan(l_a, l_b - l_a).Count('N');
                    var l_左深 = Get_平均(p_証拠.A_一意な深さ[c].AsSpan(l_a - C_両側の深さの窓, C_両側の深さの窓)) / l_中央;
                    var l_右深 = Get_平均(p_証拠.A_一意な深さ[c].AsSpan(l_b, C_両側の深さの窓)) / l_中央;
                    var l_期待 = Get_期待の組(l_a, l_b - l_N数, l_長さ - l_N数, l_断片長の標本, l_断片密度) * Math.Clamp(Math.Min(l_左深, l_右深), 0D, C_深さ比の上限);
                    var l_外れ錨 = Get_外れ錨の数(p_証拠.A_外れ錨[c], Math.Max(0, l_a - l_p99), l_a, 1, l_p50, l_長さ) + Get_外れ錨の数(p_証拠.A_外れ錨[c], l_b, Math.Min(l_長さ, l_b + l_p99), -1, l_p50, l_長さ);
                    l_候補群.Add(new 継ぎ目候補(p_名前群[c], l_a, l_b, l_N数 > 0, Get_中央値(p_証拠.A_深さ[c].AsSpan(l_a, l_b - l_a)) / l_中央, Get_跨ぐ読みの数(p_証拠.A_読み[c], l_a - C_端の許容, l_b + C_端の許容), Get_跨ぐ組の数(p_証拠.A_断片[c], l_a, l_b, l_N数, l_p1, l_p99), l_期待, l_外れ錨, Get_切れ端の数(p_証拠, c, l_a, l_b) / l_中央, l_左深, l_右深));
                }
            }

            return l_候補群;
        }

        /// <summary>
        /// 配列の中で 2 回以上現れる k-mer と N に覆われる位置に印を付ける
        /// </summary>
        /// <param name="p_配列群"></param>
        /// <returns>配列ごとの印</returns>
        public static bool[][] Get_反復の印(IReadOnlyList<string> p_配列群)
        {
            var l_全kmer = new List<ulong>(p_配列群.Sum(x => Math.Max(0, x.Length - C_反復のk長 + 1)));
            foreach (var l_配列 in p_配列群)
            {
                foreach (var (_, l_正規形) in Get_正規kmer列(l_配列))
                {
                    l_全kmer.Add(l_正規形);
                }
            }

            l_全kmer.Sort();
            HashSet<ulong> l_反復kmer = [];
            for (var i = 1; i < l_全kmer.Count; i++)
            {
                if (l_全kmer[i] == l_全kmer[i - 1])
                {
                    _ = l_反復kmer.Add(l_全kmer[i]);
                }
            }

            l_全kmer.Clear();
            var l_印 = new bool[p_配列群.Count][];
            for (var c = 0; c < p_配列群.Count; c++)
            {
                var l_配列 = p_配列群[c];
                var l_配列の印 = new bool[l_配列.Length];
                var l_有効 = new bool[Math.Max(0, l_配列.Length - C_反復のk長 + 1)];
                foreach (var (l_位置, l_正規形) in Get_正規kmer列(l_配列))
                {
                    l_有効[l_位置] = true;
                    if (l_反復kmer.Contains(l_正規形))
                    {
                        l_配列の印.AsSpan(l_位置, C_反復のk長).Fill(true);
                    }
                }

                for (var p = 0; p < l_有効.Length; p++)
                {
                    if (!l_有効[p])
                    {
                        l_配列の印.AsSpan(p, C_反復のk長).Fill(true);
                    }
                }

                for (var p = 0; p < l_配列.Length; p++)
                {
                    if (l_配列[p] == 'N')
                    {
                        l_配列の印[p] = true;
                    }
                }

                l_印[c] = l_配列の印;
            }

            return l_印;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 両端の N を落とした片を、通し番号の名前で足す (N だけの片は捨てる)
        /// </summary>
        /// <param name="p_片群"></param>
        /// <param name="p_名前"></param>
        /// <param name="p_片"></param>
        private static void V_追加_片(List<(string A_名前, string A_配列)> p_片群, string p_名前, string p_片)
        {
            var l_片 = p_片.Trim('N', 'n');
            if (l_片.Length > 0)
            {
                p_片群.Add(($"{p_名前}_{p_片群.Count + 1}", l_片));
            }
        }

        /// <summary>
        /// 配列の k-mer を、N を含まないものだけ (位置, 正規形) で順に返す
        /// </summary>
        /// <param name="p_配列"></param>
        /// <returns></returns>
        private static IEnumerable<(int A_位置, ulong A_正規形)> Get_正規kmer列(string p_配列)
        {
            const ulong l_マスク = (1UL << (C_反復のk長 * 2)) - 1UL;
            const int l_最上位へ = (C_反復のk長 - 1) * 2;
            var l_順鎖 = 0UL;
            var l_逆鎖 = 0UL;
            var l_直近のN = -1;
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_配列[i]);
                var l_Is有効 = l_塩基ID is >= Consts.塩基ID.A and <= Consts.塩基ID.T;
                var l_値 = l_Is有効 ? (ulong)(l_塩基ID - 1) : 0UL;
                l_順鎖 = ((l_順鎖 << 2) | l_値) & l_マスク;
                l_逆鎖 = (l_逆鎖 >> 2) | ((l_値 ^ 3UL) << l_最上位へ);
                if (!l_Is有効)
                {
                    l_直近のN = i;
                }

                var l_開始 = i - C_反復のk長 + 1;
                if (l_開始 >= 0 && l_直近のN < l_開始)
                {
                    yield return (l_開始, Math.Min(l_順鎖, l_逆鎖));
                }
            }
        }

        /// <summary>
        /// 反復の印のない位置での深さの中央値
        /// </summary>
        /// <param name="p_証拠"></param>
        /// <param name="p_反復の印"></param>
        /// <returns>中央値、取れなければ 1</returns>
        private static double Get_固有の深さの中央値(継ぎ目の証拠 p_証拠, bool[][] p_反復の印)
        {
            List<int> l_深さ群 = [];
            for (var c = 0; c < p_証拠.A_深さ.Length; c++)
            {
                for (var p = 0; p < p_証拠.A_深さ[c].Length; p += C_深さを取る間隔)
                {
                    if (!p_反復の印[c][p])
                    {
                        l_深さ群.Add(p_証拠.A_深さ[c][p]);
                    }
                }
            }

            var l_中央 = l_深さ群.Count > 0 ? Get_中央値(l_深さ群.ToArray()) : 0D;
            return l_中央 > 0D ? l_中央 : 1D;
        }

        /// <summary>
        /// 反復の印に、深さの厚い位置 (畳まれた反復) と一意に当たるリードの少ない位置を加え、近い区間をつないで返す
        /// </summary>
        /// <param name="p_反復の印"></param>
        /// <param name="p_深さ"></param>
        /// <param name="p_一意な深さ"></param>
        /// <param name="p_中央">固有の位置での深さの中央値</param>
        /// <returns>区間 (始まり, 終わり)</returns>
        private static List<(int A_開始, int A_終了)> Get_反復の区間(bool[] p_反復の印, int[] p_深さ, int[] p_一意な深さ, double p_中央)
        {
            var l_長さ = p_反復の印.Length;
            var l_印 = new bool[l_長さ];
            for (var p = 0; p < l_長さ; p++)
            {
                l_印[p] = p_反復の印[p] || p_深さ[p] >= C_畳まれた反復の深さ比 * p_中央 || (p_深さ[p] > 0 && p_一意な深さ[p] < C_一意な割合の下限 * p_深さ[p]);
            }

            List<(int A_開始, int A_終了)> l_区間群 = [];
            var l_位置 = 0;
            while (l_位置 < l_長さ)
            {
                if (!l_印[l_位置])
                {
                    l_位置++;
                    continue;
                }

                var l_終わり = l_位置;
                while (l_終わり < l_長さ && l_印[l_終わり])
                {
                    l_終わり++;
                }

                if (l_区間群.Count > 0 && l_位置 - l_区間群[^1].A_終了 < C_つなぐ隙間)
                {
                    l_区間群[^1] = (l_区間群[^1].A_開始, l_終わり);
                }
                else
                {
                    l_区間群.Add((l_位置, l_終わり));
                }

                l_位置 = l_終わり;
            }

            return l_区間群;
        }

        /// <summary>
        /// [p_左, p_右] を 1 本で読んだリードの数
        /// </summary>
        /// <param name="p_読み">始まりの順の整列範囲</param>
        /// <param name="p_左"></param>
        /// <param name="p_右"></param>
        /// <returns></returns>
        private static int Get_跨ぐ読みの数((int A_開始, int A_終了)[] p_読み, int p_左, int p_右)
        {
            var l_数 = 0;
            for (var i = Get_最初の位置(p_読み, p_左 - C_読みを探す戻り幅); i < p_読み.Length && p_読み[i].A_開始 <= p_左; i++)
            {
                if (p_読み[i].A_終了 >= p_右)
                {
                    l_数++;
                }
            }

            return l_数;
        }

        /// <summary>
        /// 断片が [p_左 - 端の許容, p_右 + 端の許容] を覆い、断片長から間の N の数を引いた長さが p1〜p99 に入るペアの数
        /// </summary>
        /// <param name="p_断片">始まりの順の断片</param>
        /// <param name="p_左"></param>
        /// <param name="p_右"></param>
        /// <param name="p_N数">区間の N の数 (リードは N に当たらず、見かけの断片長がその分伸びる)</param>
        /// <param name="p_p1"></param>
        /// <param name="p_p99"></param>
        /// <returns></returns>
        private static int Get_跨ぐ組の数((int A_開始, int A_終了)[] p_断片, int p_左, int p_右, int p_N数, int p_p1, int p_p99)
        {
            var l_数 = 0;
            for (var i = Get_最初の位置(p_断片, p_左 - p_p99 - p_N数); i < p_断片.Length && p_断片[i].A_開始 <= p_左 - C_端の許容; i++)
            {
                var l_長さ = p_断片[i].A_終了 - p_断片[i].A_開始 - p_N数;
                if (p_断片[i].A_終了 >= p_右 + C_端の許容 && l_長さ >= p_p1 && l_長さ <= p_p99)
                {
                    l_数++;
                }
            }

            return l_数;
        }

        /// <summary>
        /// 長さの標本の断片が [p_左 - 端の許容, p_右 + 端の許容] を覆い、かつ配列 [0, p_配列長) に収まる置き方の数の平均に、断片の密度を掛けた値
        /// </summary>
        /// <param name="p_左"></param>
        /// <param name="p_右"></param>
        /// <param name="p_配列長"></param>
        /// <param name="p_断片長の標本"></param>
        /// <param name="p_断片密度">1 塩基あたりの断片の数</param>
        /// <returns>期待の組</returns>
        private static double Get_期待の組(int p_左, int p_右, int p_配列長, int[] p_断片長の標本, double p_断片密度)
        {
            if (p_断片長の標本.Length == 0)
            {
                return 0D;
            }

            var l_左 = p_左 - C_端の許容;
            var l_右 = p_右 + C_端の許容;
            var l_合計 = 0L;
            foreach (var l_断片長 in p_断片長の標本)
            {
                l_合計 += Math.Max(0, Math.Min(l_左, p_配列長 - l_断片長) - Math.Max(0, l_右 - l_断片長) + 1);
            }

            return p_断片密度 * l_合計 / p_断片長の標本.Length;
        }

        /// <summary>
        /// [p_左, p_右) にある p_向き の外れ錨のうち、相方が来るはずの位置 (断片長の中央値だけ先) が配列の中にあるものの数
        /// </summary>
        /// <param name="p_外れ錨">位置の順の外れ錨</param>
        /// <param name="p_左"></param>
        /// <param name="p_右"></param>
        /// <param name="p_向き"></param>
        /// <param name="p_p50">断片長の中央値</param>
        /// <param name="p_配列長"></param>
        /// <returns></returns>
        private static int Get_外れ錨の数((int A_位置, int A_向き)[] p_外れ錨, int p_左, int p_右, int p_向き, int p_p50, int p_配列長)
        {
            var l_数 = 0;
            for (var i = Get_最初の位置(p_外れ錨, p_左); i < p_外れ錨.Length && p_外れ錨[i].A_位置 < p_右; i++)
            {
                var l_相方 = p_外れ錨[i].A_位置 + (p_外れ錨[i].A_向き * p_p50);
                if (p_外れ錨[i].A_向き == p_向き && l_相方 >= 0 && l_相方 < p_配列長)
                {
                    l_数++;
                }
            }

            return l_数;
        }

        /// <summary>
        /// 区間の両端の前後 切れ端の窓 で、整列が切れていたリードの数
        /// </summary>
        /// <param name="p_証拠"></param>
        /// <param name="p_配列番号"></param>
        /// <param name="p_開始"></param>
        /// <param name="p_終了"></param>
        /// <returns></returns>
        private static int Get_切れ端の数(継ぎ目の証拠 p_証拠, int p_配列番号, int p_開始, int p_終了)
        {
            var l_左 = p_証拠.A_左の切れ端[p_配列番号];
            var l_右 = p_証拠.A_右の切れ端[p_配列番号];
            var l_数 = 0;
            foreach (var l_端 in new[] { p_開始, p_終了, })
            {
                for (var x = Math.Max(0, l_端 - C_切れ端の窓); x <= Math.Min(l_左.Length - 1, l_端 + C_切れ端の窓); x++)
                {
                    l_数 += l_左[x] + l_右[x];
                }
            }

            return l_数;
        }

        /// <summary>
        /// 始まり (1 つ目の値) の順に並んだ配列で、1 つ目の値が p_下限 以上になる最初の添字
        /// </summary>
        /// <param name="p_並び"></param>
        /// <param name="p_下限"></param>
        /// <returns></returns>
        private static int Get_最初の位置((int, int)[] p_並び, int p_下限)
        {
            var l_下 = 0;
            var l_上 = p_並び.Length;
            while (l_下 < l_上)
            {
                var l_中 = (l_下 + l_上) >> 1;
                if (p_並び[l_中].Item1 < p_下限)
                {
                    l_下 = l_中 + 1;
                }
                else
                {
                    l_上 = l_中;
                }
            }

            return l_下;
        }

        /// <summary>
        /// 中央値 (偶数個なら真ん中 2 つの平均)
        /// </summary>
        /// <param name="p_値群"></param>
        /// <returns></returns>
        private static double Get_中央値(ReadOnlySpan<int> p_値群)
        {
            var l_並べ替え = p_値群.ToArray();
            Array.Sort(l_並べ替え);
            var l_中 = l_並べ替え.Length >> 1;
            return (l_並べ替え.Length & 1) == 1 ? l_並べ替え[l_中] : (l_並べ替え[l_中 - 1] + l_並べ替え[l_中]) / 2D;
        }

        /// <summary>
        /// 平均
        /// </summary>
        /// <param name="p_値群"></param>
        /// <returns></returns>
        private static double Get_平均(ReadOnlySpan<int> p_値群)
        {
            var l_合計 = 0L;
            foreach (var l_値 in p_値群)
            {
                l_合計 += l_値;
            }

            return p_値群.Length > 0 ? (double)l_合計 / p_値群.Length : 0D;
        }

        #endregion
    }
}
