using Tsumiki.Model.Foundation;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer 出現回数ヒストグラムを「誤り成分(幾何分布、裾が重い) ＋
    /// 真のk-mer成分(単一コピー平均 λ の整数倍に山を持つ、コピー数上限までの
    /// ポアソン混合)」の 2 成分混合モデルとして EM で推定する<br/>
    /// KmerHistogram の谷検出は、低カバレッジなど二峰性が視認できないデータで
    /// 壊れる(「the spectrum may not be bimodal」に落ちる)<br/>
    /// 混合モデルは谷の
    /// 目視判定に頼らず、各出現回数がどちらの成分に属するかを尤度で判定する<br/>
    /// EM は初期値に収束先が左右されるため、谷検出には頼らず複数の初期値
    /// (頻度の高い局所極大 + 対数間隔グリッド)から独立に実行し、対数尤度が
    /// 最良のものを採用する<br/>
    /// k-mer カットオフ・単一コピー基準値・tip 判定の「無条件に信頼する下限」を、
    /// 別々のヒューリスティックではなくこの 1 つのモデルから導出できるようにする
    /// (ConfigurationManager.A_スペクトルモデル として公開する)
    /// </summary>
    internal static class KmerSpectrumMixtureModel
    {
        /// <summary>
        /// 混合する真のk-mer成分のコピー数上限<br/>
        /// これを超える倍率は稀な高コピー反復と見分けがつかない
        /// </summary>
        private const int コピー数の上限 = 10;

        private const int 最大反復数 = 500;

        /// <summary>
        /// 対数尤度の変化がこれを下回ったら収束とみなす(相対値)
        /// </summary>
        private const double 収束判定 = 1e-8;

        /// <summary>
        /// 事後誤り確率がこれを下回ったら「誤りではない」と判定する有意水準
        /// </summary>
        private const double 有意水準 = 0.5;

        /// <summary>
        /// これに満たない走査範囲では、単一コピーの山とその整数倍の山を
        /// 区別する材料が無いためモデルを当てない
        /// </summary>
        private const int 最小走査範囲 = コピー数の上限 * 2;

        /// <summary>
        /// 単一コピー平均がこれを下回る場合、誤り成分(平均カバレッジが低い)との
        /// 分離が実質的にできていないとみなし、その初期値・結果を採用しない
        /// </summary>
        private const double 単一コピー平均の下限 = 3.0;

        /// <summary>
        /// 局所極大から拾う初期値候補の上限数<br/>
        /// 頻度上位のものだけを試す
        /// </summary>
        private const int 局所極大候補の上限数 = 20;

        /// <summary>
        /// コピー数別の混合比 π_k は自由な10 パラメータにせず、
        /// π_k ∝ r^(k-1) という単一の減衰率 r で表す(k=1 が最大、以降単調減少)<br/>
        /// 自由な 10 パラメータのままだと、λ' = λ/d (d は 2 以上の約数) にして
        /// k'=d, 2d, 3d... にだけ重みを乗せれば、真のピーク位置(dλ', 2dλ', ...)を
        /// そっくりそのまま再現できてしまう「調波エイリアシング」の別解にEMが
        /// 収束しうる<br/>
        /// 単調減少という 1 パラメータの制約は、真のゲノムで高コピー
        /// 配列ほど少ないという実態にも合致し、この別解(低いkを飛ばして高いkにだけ
        /// 重みが乗る形)を作れなくする
        /// </summary>
        private static double[] Get_コピー数別混合比(double p_r)
        {
            var l_重み = new double[コピー数の上限];
            var l_合計 = 0.0;
            for (var k = 0; k < コピー数の上限; k++)
            {
                l_重み[k] = Math.Pow(p_r, k);
                l_合計 += l_重み[k];
            }
            for (var k = 0; k < コピー数の上限; k++)
            {
                l_重み[k] /= l_合計;
            }
            return l_重み;
        }

        public static 混合スペクトル解析結果? Get_解析結果(
            IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限 = 10_000)
        {
            if (p_ヒストグラム.Count == 0)
            {
                return null;
            }

            var l_最大キー = p_ヒストグラム.Keys.Max();
            var l_走査上限 = (int)Math.Min(l_最大キー, p_走査上限);
            if (l_走査上限 < 最小走査範囲)
            {
                return null;
            }

            var l_出現回数 = new double[l_走査上限];
            var l_頻度 = new double[l_走査上限];
            for (var i = 0; i < l_走査上限; i++)
            {
                l_出現回数[i] = i + 1;
                l_頻度[i] = p_ヒストグラム.GetValueOrDefault((ulong)(i + 1), 0L);
            }
            var l_総数 = l_頻度.Sum();
            if (l_総数 <= 0)
            {
                return null;
            }

            var l_log階乗 = Get_log階乗テーブル(l_走査上限);

            試行結果? l_最良 = null;
            foreach (var l_初期λ in Get_初期λ候補(l_出現回数, l_頻度, l_走査上限))
            {
                var l_試行 = Get_単一試行(l_出現回数, l_頻度, l_log階乗, l_総数, l_初期λ);
                if (l_試行 is not { } l_結果)
                {
                    continue;
                }
                if (l_最良 is null || l_結果.A_対数尤度 > l_最良.Value.A_対数尤度)
                {
                    l_最良 = l_結果;
                }
            }

            if (l_最良 is not { } l_採用)
            {
                return null;
            }

            var l_事後誤り確率 = Get_事後誤り確率(
                l_出現回数, l_log階乗, l_採用.A_誤り平均, l_採用.A_誤り混合比, l_採用.A_λ,
                Get_コピー数別混合比(l_採用.A_r));

            ulong? l_カットオフ = null;
            for (var i = 0; i < l_走査上限; i++)
            {
                if (l_事後誤り確率[i] < 有意水準)
                {
                    l_カットオフ = (ulong)l_出現回数[i];
                    break;
                }
            }
            if (l_カットオフ is not { } l_確定カットオフ)
            {
                // 全域が誤り成分寄りと判定された
                // モデルが当てはまっていない
                return null;
            }

            var l_信頼下限 = 1UL;
            for (var i = l_走査上限 - 1; i >= 0; i--)
            {
                if (l_事後誤り確率[i] >= 有意水準)
                {
                    l_信頼下限 = (ulong)l_出現回数[i] + 1;
                    break;
                }
            }

            return new 混合スペクトル解析結果(
                A_単一コピー平均: l_採用.A_λ,
                A_カットオフ: l_確定カットオフ,
                A_信頼下限: l_信頼下限,
                A_誤り成分の混合比: l_採用.A_誤り混合比,
                A_誤り成分の平均: l_採用.A_誤り平均,
                A_反復回数: l_採用.A_反復回数);
        }

        private readonly record struct 試行結果(
            double A_λ, double A_誤り平均, double A_誤り混合比, double A_r,
            double A_対数尤度, int A_反復回数);

        /// <summary>
        /// 1 つの初期値から EM を収束(または反復上限)まで回す<br/>
        /// どちらかの成分が完全に空になる、または最終的な λ が下限未満になる
        /// (誤り成分と分離できていない)場合は null を返す
        /// </summary>
        private static 試行結果? Get_単一試行(
            double[] p_出現回数, double[] p_頻度, double[] p_log階乗, double p_総数, double p_初期λ)
        {
            var l_λ = p_初期λ;
            var l_誤り平均 = Get_初期誤り平均(p_出現回数, p_頻度);
            var l_誤り混合比 = 0.5;
            var l_r = 0.3;

            var l_r誤り = new double[p_出現回数.Length];
            var l_rコピー = new double[コピー数の上限][];
            for (var k = 0; k < コピー数の上限; k++)
            {
                l_rコピー[k] = new double[p_出現回数.Length];
            }

            var l_前回対数尤度 = double.NegativeInfinity;
            var l_対数尤度 = double.NegativeInfinity;
            int l_反復数;
            for (l_反復数 = 1; l_反復数 <= 最大反復数; l_反復数++)
            {
                var l_コピー数別混合比 = Get_コピー数別混合比(l_r);
                l_対数尤度 = Get_Estep(
                    p_出現回数, p_頻度, p_log階乗, l_誤り平均, l_誤り混合比, l_λ, l_コピー数別混合比,
                    l_r誤り, l_rコピー);

                if (double.IsNaN(l_対数尤度) || double.IsInfinity(l_対数尤度))
                {
                    return null;
                }

                if (!Get_Mstep(
                    p_出現回数, p_頻度, l_r誤り, l_rコピー, p_総数,
                    ref l_誤り平均, ref l_誤り混合比, ref l_λ, ref l_r))
                {
                    return null;
                }

                var l_基準 = Math.Abs(l_前回対数尤度) < 1 ? 1 : Math.Abs(l_前回対数尤度);
                if (Math.Abs(l_対数尤度 - l_前回対数尤度) < 収束判定 * l_基準)
                {
                    break;
                }
                l_前回対数尤度 = l_対数尤度;
            }

            return double.IsNaN(l_λ) || double.IsInfinity(l_λ) || l_λ < 単一コピー平均の下限 ? null : new 試行結果(l_λ, l_誤り平均, l_誤り混合比, l_r, l_対数尤度, l_反復数);
        }

        /// <summary>
        /// 誤り成分の初期平均<br/>
        /// 出現回数 1..20 の頻度加重平均を使う
        /// (低頻度域の大まかな水準を見るだけの粗い初期値で、EMが精密化する)<br/>
        /// 定数 1 で初期化すると、幾何分布が c=1 にのみ質量を持つ退化点になり、
        /// そこでは c&gt;=2 の責任度が最初から実質ゼロになって動けなくなる
        /// (EMが抜け出せない自明な吸収点に落ちる)
        /// </summary>
        private static double Get_初期誤り平均(double[] p_出現回数, double[] p_頻度)
        {
            var l_上限 = Math.Min(20, p_出現回数.Length);
            var l_重み合計 = 0.0;
            var l_加重和 = 0.0;
            for (var i = 0; i < l_上限; i++)
            {
                l_重み合計 += p_頻度[i];
                l_加重和 += p_頻度[i] * p_出現回数[i];
            }
            return l_重み合計 > 0 ? Math.Max(1.5, l_加重和 / l_重み合計) : 2.0;
        }

        /// <summary>
        /// EM の初期値候補<br/>
        /// 頻度の高い局所極大(データが実際に示す山)と、
        /// 対数間隔のグリッド(局所極大が谷に埋もれて見えない場合の保険)を併用する<br/>
        /// 谷の目視判定には依存しない
        /// </summary>
        private static List<double> Get_初期λ候補(double[] p_出現回数, double[] p_頻度, int p_走査上限)
        {
            var l_極大インデックス = new List<int>();
            for (var i = 1; i < p_走査上限 - 1; i++)
            {
                if (p_頻度[i] > p_頻度[i - 1] && p_頻度[i] > p_頻度[i + 1] && p_頻度[i] > 0)
                {
                    l_極大インデックス.Add(i);
                }
            }

            var l_候補 = l_極大インデックス
                .OrderByDescending(i => p_頻度[i])
                .Take(局所極大候補の上限数)
                .Select(i => p_出現回数[i])
                .ToList();

            for (var l_値 = 単一コピー平均の下限; l_値 <= p_走査上限; l_値 *= 2)
            {
                l_候補.Add(l_値);
            }

            return [.. l_候補.Where(x => x >= 単一コピー平均の下限).Distinct()];
        }

        /// <summary>
        /// E-step<br/>
        /// 各出現回数について、誤り成分・コピー数 1..上限の各成分への
        /// 事後責任(責任度)を計算し、対数尤度を返す
        /// </summary>
        private static double Get_Estep(
            double[] p_出現回数, double[] p_頻度, double[] p_log階乗,
            double p_誤り平均, double p_誤り混合比, double p_λ, double[] p_コピー数別混合比,
            double[] p_r誤り, double[][] p_rコピー)
        {
            var l_p = 1.0 / p_誤り平均;
            var l_log1マイナスp = Math.Log(Math.Max(1e-300, 1 - l_p));
            var l_logP = Math.Log(Math.Max(1e-300, l_p));
            var l_logw誤り = Math.Log(Math.Max(1e-300, p_誤り混合比));
            var l_logw真 = Math.Log(Math.Max(1e-300, 1 - p_誤り混合比));

            var l_対数尤度 = 0.0;
            var l_項 = new double[1 + コピー数の上限];

            for (var i = 0; i < p_出現回数.Length; i++)
            {
                var c = p_出現回数[i];

                // 幾何分布(誤り成分): P(c) = (1-p)^(c-1) * p
                l_項[0] = l_logw誤り + (c - 1) * l_log1マイナスp + l_logP;

                for (var k = 1; k <= コピー数の上限; k++)
                {
                    var l_μ = k * p_λ;
                    // ポアソン分布(コピー数kのゲノム成分): P(c) = exp(-μ) μ^c / c!
                    l_項[k] = l_logw真 + Math.Log(Math.Max(1e-300, p_コピー数別混合比[k - 1]))
                        + (-l_μ + c * Math.Log(l_μ) - p_log階乗[(int)c]);
                }

                var l_logP_c = Get_LogSumExp(l_項);
                if (p_頻度[i] > 0)
                {
                    l_対数尤度 += p_頻度[i] * l_logP_c;
                }

                p_r誤り[i] = Math.Exp(l_項[0] - l_logP_c);
                for (var k = 1; k <= コピー数の上限; k++)
                {
                    p_rコピー[k - 1][i] = Math.Exp(l_項[k] - l_logP_c);
                }
            }

            return l_対数尤度;
        }

        /// <summary>
        /// M-step<br/>
        /// 責任度で重み付けした最尤推定でパラメータを更新する<br/>
        /// どちらかの成分が完全に空(重みの総和が 0)になった場合は false を返す<br/>
        /// r(コピー数の減衰率)は、コピー数の重み付き平均 m から
        /// (幾何分布の平均 = 1/(1-r) の関係を使って) r = 1 - 1/m として更新する<br/>
        /// 打ち切り(コピー数上限 10)の影響を無視した近似だが、r が極端に 1 に
        /// 近くない限り無視できる誤差であり、閉形式で軽量に更新できる
        /// </summary>
        private static bool Get_Mstep(
            double[] p_出現回数, double[] p_頻度, double[] p_r誤り, double[][] p_rコピー, double p_総数,
            ref double p_誤り平均, ref double p_誤り混合比, ref double p_λ, ref double p_r)
        {
            var l_誤りの重み合計 = 0.0;
            var l_誤りの加重カウント合計 = 0.0;
            for (var i = 0; i < p_出現回数.Length; i++)
            {
                var w = p_頻度[i] * p_r誤り[i];
                l_誤りの重み合計 += w;
                l_誤りの加重カウント合計 += w * p_出現回数[i];
            }

            var l_コピー別重み合計 = new double[コピー数の上限];
            var l_コピー別加重カウント合計 = new double[コピー数の上限];
            var l_ゲノム成分の重み合計 = 0.0;
            for (var k = 0; k < コピー数の上限; k++)
            {
                for (var i = 0; i < p_出現回数.Length; i++)
                {
                    var w = p_頻度[i] * p_rコピー[k][i];
                    l_コピー別重み合計[k] += w;
                    l_コピー別加重カウント合計[k] += w * p_出現回数[i];
                }
                l_ゲノム成分の重み合計 += l_コピー別重み合計[k];
            }

            if (l_誤りの重み合計 <= 0 || l_ゲノム成分の重み合計 <= 0)
            {
                return false;
            }

            p_誤り混合比 = l_誤りの重み合計 / p_総数;
            p_誤り平均 = Math.Max(1.0, l_誤りの加重カウント合計 / l_誤りの重み合計);

            var l_分子 = 0.0;
            var l_コピー数の加重合計 = 0.0;
            for (var k = 0; k < コピー数の上限; k++)
            {
                l_分子 += l_コピー別加重カウント合計[k];
                l_コピー数の加重合計 += (k + 1) * l_コピー別重み合計[k];
            }
            if (l_コピー数の加重合計 <= 0)
            {
                return false;
            }
            p_λ = l_分子 / l_コピー数の加重合計;

            var l_平均コピー数 = l_コピー数の加重合計 / l_ゲノム成分の重み合計;
            p_r = Math.Clamp(1 - 1 / Math.Max(1.0, l_平均コピー数), 1e-6, 1 - 1e-6);

            return true;
        }

        /// <summary>
        /// 最終パラメータでの、各出現回数における事後誤り確率 P(誤り成分 | c)
        /// </summary>
        private static double[] Get_事後誤り確率(
            double[] p_出現回数, double[] p_log階乗,
            double p_誤り平均, double p_誤り混合比, double p_λ, double[] p_コピー数別混合比)
        {
            var l_p = 1.0 / p_誤り平均;
            var l_log1マイナスp = Math.Log(Math.Max(1e-300, 1 - l_p));
            var l_logP = Math.Log(Math.Max(1e-300, l_p));
            var l_logw誤り = Math.Log(Math.Max(1e-300, p_誤り混合比));
            var l_logw真 = Math.Log(Math.Max(1e-300, 1 - p_誤り混合比));

            var l_結果 = new double[p_出現回数.Length];
            var l_項 = new double[1 + コピー数の上限];
            for (var i = 0; i < p_出現回数.Length; i++)
            {
                var c = p_出現回数[i];
                l_項[0] = l_logw誤り + (c - 1) * l_log1マイナスp + l_logP;
                for (var k = 1; k <= コピー数の上限; k++)
                {
                    var l_μ = k * p_λ;
                    l_項[k] = l_logw真 + Math.Log(Math.Max(1e-300, p_コピー数別混合比[k - 1]))
                        + (-l_μ + c * Math.Log(l_μ) - p_log階乗[(int)c]);
                }
                var l_logP_c = Get_LogSumExp(l_項);
                l_結果[i] = Math.Exp(l_項[0] - l_logP_c);
            }
            return l_結果;
        }

        private static double Get_LogSumExp(double[] p_対数値)
        {
            var l_最大 = p_対数値.Max();
            if (double.IsNegativeInfinity(l_最大))
            {
                return double.NegativeInfinity;
            }
            var l_合計 = 0.0;
            foreach (var l_値 in p_対数値)
            {
                l_合計 += Math.Exp(l_値 - l_最大);
            }
            return l_最大 + Math.Log(l_合計);
        }

        /// <summary>
        /// c=0..p_上限 の log(c!) の表<br/>
        /// ポアソン対数尤度の計算に使う
        /// </summary>
        private static double[] Get_log階乗テーブル(int p_上限)
        {
            var l_表 = new double[p_上限 + 1];
            l_表[0] = 0;
            for (var i = 1; i <= p_上限; i++)
            {
                l_表[i] = l_表[i - 1] + Math.Log(i);
            }
            return l_表;
        }
    }
}
