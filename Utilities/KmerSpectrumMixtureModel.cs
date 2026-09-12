using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// k-mer 出現回数ヒストグラムを「誤り成分 (幾何分布、裾が重い) ＋真の k-mer 成分 (単一コピー平均 λ の整数倍に山を持つ、コピー数上限までのポアソン混合) 」の 2 成分混合モデルとして EM で推定する
    /// </summary>
    internal static class KmerSpectrumMixtureModel
    {
        #region 定数

        /// <summary>
        /// 混合する真の k-mer 成分のコピー数上限
        /// </summary>
        /// <remarks>
        /// これを超える倍率は稀な高コピー反復と見分けがつかない
        /// </remarks>
        private const int コピー数の上限 = 10;

        /// <summary>
        /// 最大反復数
        /// </summary>
        private const int 最大反復数 = 500;

        /// <summary>
        /// 対数尤度の変化がこれを下回ったら収束とみなす (相対値)
        /// </summary>
        private const double 収束判定 = 1e-8D;

        /// <summary>
        /// 事後誤り確率がこれを下回ったら「誤りではない」と判定する有意水準
        /// </summary>
        private const double 有意水準 = 0.5D;

        /// <summary>
        /// これに満たない走査範囲では、単一コピーの山とその整数倍の山を区別する材料が無いためモデルを当てない
        /// </summary>
        private const int 最小走査範囲 = コピー数の上限 << 1;

        /// <summary>
        /// 単一コピー平均がこれを下回る場合、誤り成分 (平均カバレッジが低い) との分離が実質的にできていないとみなし、その初期値・結果を採用しない
        /// </summary>
        private const double 単一コピー平均の下限 = 3D;

        /// <summary>
        /// 局所極大から拾う初期値候補の上限数
        /// </summary>
        /// <remarks>
        /// 頻度上位のものだけを試す
        /// </remarks>
        private const int 局所極大候補の上限数 = 20;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// k-mer スペクトルへ 2 成分の混合モデルを当てはめて返す
        /// </summary>
        /// <param name="p_ヒストグラム">出現回数ごとの k-mer 種類数</param>
        /// <param name="p_走査上限">当てはめに使う出現回数の上限</param>
        /// <returns>当てはめた結果、収束しなければ null</returns>
        public static 混合スペクトル解析結果? Get_解析結果(IReadOnlyDictionary<ulong, long> p_ヒストグラム, ulong p_走査上限 = 10_000UL)
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
            if (l_総数 <= 0D)
            {
                return null;
            }

            var l_log階乗 = Get_log階乗テーブル(l_走査上限);

            ヒストグラム試行結果? l_最良 = null;
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

            var l_事後誤り確率 = Get_事後誤り確率(l_出現回数, l_log階乗, l_採用.A_誤り平均, l_採用.A_誤り混合比, l_採用.A_λ, Get_コピー数別混合比(l_採用.A_コピー数減衰率));

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
                    l_信頼下限 = (ulong)l_出現回数[i] + 1UL;
                    break;
                }
            }

            return new 混合スペクトル解析結果(A_単一コピー平均: l_採用.A_λ, A_カットオフ: l_確定カットオフ, A_信頼下限: l_信頼下限, A_誤り成分の混合比: l_採用.A_誤り混合比, A_誤り成分の平均: l_採用.A_誤り平均, A_反復回数: l_採用.A_反復回数);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// コピー数別の混合比 π_k は自由な 10 パラメータにせず、π_k ∝ r^ (k-1) という単一の減衰率 r で表す (k=1 が最大、以降単調減少)
        /// </summary>
        /// <param name="p_コピー数減衰率">コピー数別混合比の減衰率</param>
        /// <remarks>
        /// 自由な 10 パラメータのままだと、λ' = λ/d (d は 2 以上の約数) にして k'=d, 2 d, 3 d... にだけ重みを乗せれば、真のピーク位置 (dλ', 2 dλ', ...) をそっくりそのまま再現できてしまう「調波エイリアシング」の別解に EM が収束しうる<br/>
        /// 単調減少という 1 パラメータの制約は、真のゲノムで高コピー配列ほど少ないという実態にも合致し、この別解 (低い k を飛ばして高い k にだけ重みが乗る形) を作れなくする
        /// </remarks>
        /// <returns></returns>
        private static double[] Get_コピー数別混合比(double p_コピー数減衰率)
        {
            var l_重み = new double[コピー数の上限];
            var l_合計 = 0D;
            for (var k = 0; k < コピー数の上限; k++)
            {
                l_重み[k] = Math.Pow(p_コピー数減衰率, k);
                l_合計 += l_重み[k];
            }
            for (var k = 0; k < コピー数の上限; k++)
            {
                l_重み[k] /= l_合計;
            }
            return l_重み;
        }

        /// <summary>
        /// 1 つの初期値から EM を収束 (または反復上限) まで回す
        /// </summary>
        /// <param name="p_出現回数"></param>
        /// <param name="p_頻度"></param>
        /// <param name="p_log階乗"></param>
        /// <param name="p_総数"></param>
        /// <param name="p_初期λ"></param>
        /// <remarks>
        /// どちらかの成分が完全に空になる、または最終的な λ が下限未満になる (誤り成分と分離できていない) 場合は null を返す
        /// </remarks>
        /// <returns></returns>
        private static ヒストグラム試行結果? Get_単一試行(double[] p_出現回数, double[] p_頻度, double[] p_log階乗, double p_総数, double p_初期λ)
        {
            var l_λ = p_初期λ;
            var l_誤り平均 = Get_初期誤り平均(p_出現回数, p_頻度);
            var l_誤り混合比 = 0.5D;
            var l_コピー数減衰率 = 0.3D;

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
                var l_コピー数別混合比 = Get_コピー数別混合比(l_コピー数減衰率);
                l_対数尤度 = Get_Estep(p_出現回数, p_頻度, p_log階乗, l_誤り平均, l_誤り混合比, l_λ, l_コピー数別混合比, l_r誤り, l_rコピー);

                if (double.IsNaN(l_対数尤度) || double.IsInfinity(l_対数尤度))
                {
                    return null;
                }

                if (!Try更新_Mstep(p_出現回数, p_頻度, l_r誤り, l_rコピー, p_総数, ref l_誤り平均, ref l_誤り混合比, ref l_λ, ref l_コピー数減衰率))
                {
                    return null;
                }

                var l_基準 = Math.Abs(l_前回対数尤度) < 1D ? 1D : Math.Abs(l_前回対数尤度);
                if (Math.Abs(l_対数尤度 - l_前回対数尤度) < 収束判定 * l_基準)
                {
                    break;
                }
                l_前回対数尤度 = l_対数尤度;
            }

            return double.IsNaN(l_λ) || double.IsInfinity(l_λ) || l_λ < 単一コピー平均の下限 ? null : new ヒストグラム試行結果(l_λ, l_誤り平均, l_誤り混合比, l_コピー数減衰率, l_対数尤度, l_反復数);
        }

        /// <summary>
        /// 誤り成分の初期平均
        /// </summary>
        /// <param name="p_出現回数"></param>
        /// <param name="p_頻度"></param>
        /// <remarks>
        /// 出現回数 1..20 の頻度加重平均を使う (低頻度域の大まかな水準を見るだけの粗い初期値で、EM が精密化する) <br/>
        /// 定数 1 で初期化すると、幾何分布が c=1 にのみ質量を持つ退化点になり、そこでは c&gt;=2 の責任度が最初から実質ゼロになって動けなくなる (EM が抜け出せない自明な吸収点に落ちる)
        /// </remarks>
        /// <returns></returns>
        private static double Get_初期誤り平均(double[] p_出現回数, double[] p_頻度)
        {
            var l_上限 = Math.Min(20, p_出現回数.Length);
            var l_重み合計 = 0D;
            var l_加重和 = 0D;
            for (var i = 0; i < l_上限; i++)
            {
                l_重み合計 += p_頻度[i];
                l_加重和 += p_頻度[i] * p_出現回数[i];
            }
            return l_重み合計 > 0D ? Math.Max(1.5D, l_加重和 / l_重み合計) : 2.0D;
        }

        /// <summary>
        /// EM の初期値候補
        /// </summary>
        /// <param name="p_出現回数"></param>
        /// <param name="p_頻度"></param>
        /// <param name="p_走査上限"></param>
        /// <remarks>
        /// 頻度の高い局所極大 (データが実際に示す山) と、対数間隔のグリッド (局所極大が谷に埋もれて見えない場合の保険) を併用する<br/>
        /// 谷の目視判定には依存しない
        /// </remarks>
        /// <returns></returns>
        private static List<double> Get_初期λ候補(double[] p_出現回数, double[] p_頻度, int p_走査上限)
        {
            var l_極大インデックス = new List<int>();
            for (var i = 1; i < p_走査上限 - 1; i++)
            {
                if (p_頻度[i] > p_頻度[i - 1] && p_頻度[i] > p_頻度[i + 1] && p_頻度[i] > 0D)
                {
                    l_極大インデックス.Add(i);
                }
            }

            var l_候補 = l_極大インデックス
                .OrderByDescending(i => p_頻度[i])
                .Take(局所極大候補の上限数)
                .Select(i => p_出現回数[i])
                .ToList();

            for (var l_値 = 単一コピー平均の下限; l_値 <= p_走査上限; l_値 *= 2D)
            {
                l_候補.Add(l_値);
            }

            return [.. l_候補.Where(x => x >= 単一コピー平均の下限).Distinct()];
        }

        /// <summary>
        /// E-step
        /// </summary>
        /// <param name="p_出現回数"></param>
        /// <param name="p_頻度"></param>
        /// <param name="p_log階乗"></param>
        /// <param name="p_誤り平均"></param>
        /// <param name="p_誤り混合比"></param>
        /// <param name="p_λ"></param>
        /// <param name="p_コピー数別混合比"></param>
        /// <param name="p_r誤り">出現回数ごとの誤り成分への責任度、この呼び出しで書き込む</param>
        /// <param name="p_rコピー">出現回数・コピー数ごとの責任度、この呼び出しで書き込む</param>
        /// <remarks>
        /// 各出現回数について、誤り成分・コピー数 1..上限の各成分への事後責任 (責任度) を計算し、対数尤度を返す
        /// </remarks>
        /// <returns></returns>
        private static double Get_Estep(double[] p_出現回数, double[] p_頻度, double[] p_log階乗, double p_誤り平均, double p_誤り混合比, double p_λ, double[] p_コピー数別混合比, double[] p_r誤り, double[][] p_rコピー)
        {
            var l_p = 1D / p_誤り平均;
            var l_log1マイナスp = Math.Log(Math.Max(1e-300D, 1D - l_p));
            var l_logP = Math.Log(Math.Max(1e-300D, l_p));
            var l_logw誤り = Math.Log(Math.Max(1e-300D, p_誤り混合比));
            var l_logw真 = Math.Log(Math.Max(1e-300D, 1D - p_誤り混合比));

            var l_対数尤度 = 0D;
            var l_項 = new double[1 + コピー数の上限];

            for (var i = 0; i < p_出現回数.Length; i++)
            {
                var l_末尾配列 = p_出現回数[i];

                // 幾何分布 (誤り成分) : P (c) = (1-p) ^ (c-1) * p
                l_項[0] = l_logw誤り + (l_末尾配列 - 1D) * l_log1マイナスp + l_logP;

                for (var k = 1; k <= コピー数の上限; k++)
                {
                    var l_μ = k * p_λ;
                    // ポアソン分布 (コピー数 k のゲノム成分) : P (c) = exp (-μ) μ^c / c!
                    l_項[k] = l_logw真 + Math.Log(Math.Max(1e-300D, p_コピー数別混合比[k - 1]))
                        + (-l_μ + l_末尾配列 * Math.Log(l_μ) - p_log階乗[(int)l_末尾配列]);
                }

                var l_logP_c = Get_LogSumExp(l_項);
                if (p_頻度[i] > 0D)
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
        /// M-step
        /// </summary>
        /// <param name="p_出現回数"></param>
        /// <param name="p_頻度"></param>
        /// <param name="p_r誤り"></param>
        /// <param name="p_rコピー"></param>
        /// <param name="p_総数"></param>
        /// <param name="p_誤り平均">更新後の誤り成分の平均</param>
        /// <param name="p_誤り混合比">更新後の誤り成分の混合比</param>
        /// <param name="p_λ">更新後の単一コピー平均</param>
        /// <param name="p_コピー数減衰率">更新後のコピー数別混合比の減衰率</param>
        /// <returns></returns>
        private static bool Try更新_Mstep(double[] p_出現回数, double[] p_頻度, double[] p_r誤り, double[][] p_rコピー, double p_総数, ref double p_誤り平均, ref double p_誤り混合比, ref double p_λ, ref double p_コピー数減衰率)
        {
            var l_誤りの重み合計 = 0D;
            var l_誤りの加重カウント合計 = 0D;
            for (var i = 0; i < p_出現回数.Length; i++)
            {
                var l_終点 = p_頻度[i] * p_r誤り[i];
                l_誤りの重み合計 += l_終点;
                l_誤りの加重カウント合計 += l_終点 * p_出現回数[i];
            }

            var l_コピー別重み合計 = new double[コピー数の上限];
            var l_コピー別加重カウント合計 = new double[コピー数の上限];
            var l_ゲノム成分の重み合計 = 0D;
            for (var k = 0; k < コピー数の上限; k++)
            {
                for (var i = 0; i < p_出現回数.Length; i++)
                {
                    var l_終点 = p_頻度[i] * p_rコピー[k][i];
                    l_コピー別重み合計[k] += l_終点;
                    l_コピー別加重カウント合計[k] += l_終点 * p_出現回数[i];
                }
                l_ゲノム成分の重み合計 += l_コピー別重み合計[k];
            }

            if (l_誤りの重み合計 <= 0D || l_ゲノム成分の重み合計 <= 0D)
            {
                return false;
            }

            p_誤り混合比 = l_誤りの重み合計 / p_総数;
            p_誤り平均 = Math.Max(1D, l_誤りの加重カウント合計 / l_誤りの重み合計);

            var l_分子 = 0D;
            var l_コピー数の加重合計 = 0D;
            for (var k = 0; k < コピー数の上限; k++)
            {
                l_分子 += l_コピー別加重カウント合計[k];
                l_コピー数の加重合計 += (k + 1) * l_コピー別重み合計[k];
            }
            if (l_コピー数の加重合計 <= 0D)
            {
                return false;
            }
            p_λ = l_分子 / l_コピー数の加重合計;

            var l_平均コピー数 = l_コピー数の加重合計 / l_ゲノム成分の重み合計;
            p_コピー数減衰率 = Math.Clamp(1D - 1D / Math.Max(1D, l_平均コピー数), 1e-6D, 1D - 1e-6D);

            return true;
        }

        /// <summary>
        /// 最終パラメータでの、各出現回数における事後誤り確率 P (誤り成分 | c)
        /// </summary>
        /// <param name="p_出現回数"></param>
        /// <param name="p_log階乗"></param>
        /// <param name="p_誤り平均"></param>
        /// <param name="p_誤り混合比"></param>
        /// <param name="p_λ"></param>
        /// <param name="p_コピー数別混合比"></param>
        /// <returns></returns>
        private static double[] Get_事後誤り確率(double[] p_出現回数, double[] p_log階乗, double p_誤り平均, double p_誤り混合比, double p_λ, double[] p_コピー数別混合比)
        {
            var l_p = 1.0D / p_誤り平均;
            var l_log1マイナスp = Math.Log(Math.Max(1e-300D, 1D - l_p));
            var l_logP = Math.Log(Math.Max(1e-300D, l_p));
            var l_logw誤り = Math.Log(Math.Max(1e-300D, p_誤り混合比));
            var l_logw真 = Math.Log(Math.Max(1e-300D, 1D - p_誤り混合比));

            var l_結果 = new double[p_出現回数.Length];
            var l_項 = new double[1 + コピー数の上限];
            for (var i = 0; i < p_出現回数.Length; i++)
            {
                var l_末尾配列 = p_出現回数[i];
                l_項[0] = l_logw誤り + (l_末尾配列 - 1D) * l_log1マイナスp + l_logP;
                for (var k = 1; k <= コピー数の上限; k++)
                {
                    var l_μ = k * p_λ;
                    l_項[k] = l_logw真 + Math.Log(Math.Max(1e-300D, p_コピー数別混合比[k - 1]))
                        + (-l_μ + l_末尾配列 * Math.Log(l_μ) - p_log階乗[(int)l_末尾配列]);
                }
                var l_logP_c = Get_LogSumExp(l_項);
                l_結果[i] = Math.Exp(l_項[0] - l_logP_c);
            }
            return l_結果;
        }

        /// <summary>
        /// 対数のまま和を取る
        /// </summary>
        /// <remarks>
        /// そのまま指数へ戻すと桁が溢れるため、最大値を括り出してから足す
        /// </remarks>
        /// <param name="p_対数値">足し合わせる対数値</param>
        /// <returns>和の対数</returns>
        private static double Get_LogSumExp(double[] p_対数値)
        {
            var l_最大 = p_対数値.Max();
            if (double.IsNegativeInfinity(l_最大))
            {
                return double.NegativeInfinity;
            }
            var l_合計 = 0D;
            foreach (var l_値 in p_対数値)
            {
                l_合計 += Math.Exp(l_値 - l_最大);
            }
            return l_最大 + Math.Log(l_合計);
        }

        /// <summary>
        /// c=0..p_上限 の log (c!) の表
        /// </summary>
        /// <param name="p_上限"></param>
        /// <remarks>
        /// ポアソン対数尤度の計算に使う
        /// </remarks>
        /// <returns></returns>
        private static double[] Get_log階乗テーブル(int p_上限)
        {
            var l_表 = new double[p_上限 + 1];
            l_表[0] = 0D;
            for (var i = 1; i <= p_上限; i++)
            {
                l_表[i] = l_表[i - 1] + Math.Log(i);
            }
            return l_表;
        }

        #endregion
    }
}
