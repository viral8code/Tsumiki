using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer スペクトルの 2 成分混合モデル (誤り=幾何分布、真の k-mer=単一コピー平均の
    /// 整数倍に山を持つポアソン混合) の EM 推定を、理論分布そのものから作った
    /// ヒストグラムで固定する<br/>
    /// 理論分布を使うことで、サンプリング由来のノイズを
    /// 排して「モデルが正しいパラメータへ収束するか」だけを検証できる
    /// </summary>
    public class KmerSpectrumMixtureModelTests
    {
        private static double Get_幾何分布確率(double p_出現回数, double p_平均)
        {
            var l_p = 1.0 / p_平均;
            return Math.Pow(1 - l_p, p_出現回数 - 1) * l_p;
        }

        private static double Get_ポアソン確率(double p_出現回数, double p_μ)
        {
            var l_log階乗 = 0.0;
            for (var i = 2; i <= p_出現回数; i++)
            {
                l_log階乗 += Math.Log(i);
            }
            return Math.Exp(-p_μ + p_出現回数 * Math.Log(p_μ) - l_log階乗);
        }

        /// <summary>
        /// 指定したパラメータ通りの理論混合分布から、その通りのヒストグラムを作る
        /// (サンプリングはしない<br/>
        /// EM が真のパラメータへ収束するかだけを見るため)
        /// </summary>
        private static Dictionary<ulong, long> Get_理論ヒストグラム(
            double p_λ, double p_誤り混合比, double p_誤り平均, double[] p_コピー数別混合比,
            int p_上限, long p_総数)
        {
            Dictionary<ulong, long> l_ヒストグラム = [];
            for (var c = 1UL; c <= (ulong)p_上限; c++)
            {
                var l_確率 = p_誤り混合比 * Get_幾何分布確率(c, p_誤り平均);
                for (var k = 1; k <= p_コピー数別混合比.Length; k++)
                {
                    l_確率 += (1 - p_誤り混合比) * p_コピー数別混合比[k - 1] * Get_ポアソン確率(c, k * p_λ);
                }
                l_ヒストグラム[c] = (long)Math.Round(l_確率 * p_総数);
            }
            return l_ヒストグラム;
        }

        private static double[] Get_コピー数別混合比_単一コピー優勢()
        {
            // モデル (π_k ∝ r^(k-1))と同じ形の生成分布
            const double r = 0.15;
            var l_比 = new double[10];
            var l_合計 = 0.0;
            for (var k = 0; k < 10; k++)
            {
                l_比[k] = Math.Pow(r, k);
                l_合計 += l_比[k];
            }
            for (var k = 0; k < 10; k++)
            {
                l_比[k] /= l_合計;
            }
            return l_比;
        }

        [Fact]
        public void Get_解析結果_EmptyHistogram_ReturnsNull()
        {
            Assert.Null(KmerSpectrumMixtureModel.Get_解析結果(new Dictionary<ulong, long>()));
        }

        [Fact]
        public void Get_解析結果_TooNarrowAScanRange_ReturnsNull()
        {
            // コピー数上限 (10) の 2 倍に満たない走査範囲では、単一コピーの山と
            // その倍数の山を区別する材料が無い
            Dictionary<ulong, long> l_ヒストグラム = new() { [1] = 100, [2] = 50, [5] = 200 };
            Assert.Null(KmerSpectrumMixtureModel.Get_解析結果(l_ヒストグラム));
        }

        [Fact]
        public void Get_解析結果_RecoversSingleCopyMeanFromATheoreticalMixture()
        {
            const double 真のλ = 30.0;
            const double 真の誤り平均 = 3.0;
            const double 真の誤り混合比 = 0.35;
            var l_コピー数別混合比 = Get_コピー数別混合比_単一コピー優勢();

            var l_ヒストグラム = Get_理論ヒストグラム(
                真のλ, 真の誤り混合比, 真の誤り平均, l_コピー数別混合比, p_上限: 300, p_総数: 1_000_000);

            var l_結果 = KmerSpectrumMixtureModel.Get_解析結果(l_ヒストグラム);

            Assert.NotNull(l_結果);
            Assert.InRange(l_結果!.A_単一コピー平均, 真のλ - 2, 真のλ + 2);
            // カットオフは誤り成分側、信頼下限はカットオフ以上、どちらも単一コピー峰 (30) 未満のはず
            Assert.InRange((double)l_結果.A_カットオフ, 1, 真のλ);
            Assert.True(l_結果.A_信頼下限 >= l_結果.A_カットオフ);
        }

        /// <summary>
        /// 低カバレッジ気味 (単一コピー平均が誤り成分の平均に近い) でも、
        /// 谷が視認できるかどうかに関わらずモデルが分離できることを確かめる
        /// </summary>
        [Fact]
        public void Get_解析結果_SeparatesComponentsEvenAtLowCoverage()
        {
            const double 真のλ = 12.0;
            const double 真の誤り平均 = 2.0;
            const double 真の誤り混合比 = 0.5;
            var l_コピー数別混合比 = Get_コピー数別混合比_単一コピー優勢();

            var l_ヒストグラム = Get_理論ヒストグラム(
                真のλ, 真の誤り混合比, 真の誤り平均, l_コピー数別混合比, p_上限: 150, p_総数: 500_000);

            var l_結果 = KmerSpectrumMixtureModel.Get_解析結果(l_ヒストグラム);

            Assert.NotNull(l_結果);
            Assert.InRange(l_結果!.A_単一コピー平均, 真のλ - 2, 真のλ + 2);
        }
    }
}
