using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// k-mer スペクトルの 2 成分混合モデル (誤り=幾何分布、真の k-mer=単一コピー平均の
    /// 整数倍に山を持つポアソン混合) の EM 推定を、理論分布そのものから作った
    /// ヒストグラムで固定する
    /// </summary>
    /// <remarks>
    /// 理論分布を使うことで、サンプリング由来のノイズを
    /// 排して「モデルが正しいパラメータへ収束するか」だけを検証できる
    /// </remarks>
    public class KmerSpectrumMixtureModelTests
    {
        /// <summary>
        /// 幾何分布の確率を返す
        /// </summary>
        /// <param name="p_出現回数">出現回数</param>
        /// <param name="p_平均">分布の平均</param>
        /// <returns>確率</returns>
        private static double Get_幾何分布確率(double p_出現回数, double p_平均)
        {
            var l_p = 1.0 / p_平均;
            return Math.Pow(1 - l_p, p_出現回数 - 1) * l_p;
        }

        /// <summary>
        /// ポアソン分布の確率を返す
        /// </summary>
        /// <param name="p_出現回数">出現回数</param>
        /// <param name="p_μ">分布の平均</param>
        /// <returns>確率</returns>
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
        /// (サンプリングはしない
        /// </summary>
        /// <remarks>
        /// EM が真のパラメータへ収束するかだけを見るため)
        /// </remarks>
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

        /// <summary>
        /// 単一コピーが優勢なコピー数ごとの混合比を返す
        /// </summary>
        /// <returns>コピー数ごとの混合比</returns>
        private static double[] Get_コピー数別混合比_単一コピー優勢()
        {
            // モデル (π_k ∝ r^(k-1))と同じ形の生成分布
            const double r = 0.15D;
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

        /// <summary>
        /// 空のヒストグラムでは null を返す
        /// </summary>
        [Fact]
        public void 空のヒストグラムではnullを返す()
        {
            Assert.Null(KmerSpectrumMixtureModel.Get_解析結果(new Dictionary<ulong, long>()));
        }

        /// <summary>
        /// 走査範囲が狭すぎる場合は null を返す
        /// </summary>
        [Fact]
        public void 走査範囲が狭すぎる場合はnullを返す()
        {
            // コピー数上限 (10) の 2 倍に満たない走査範囲では、単一コピーの山と
            // その倍数の山を区別する材料が無い
            Dictionary<ulong, long> l_ヒストグラム = new() { [1] = 100, [2] = 50, [5] = 200 };
            Assert.Null(KmerSpectrumMixtureModel.Get_解析結果(l_ヒストグラム));
        }

        /// <summary>
        /// 理論混合分布から、単一コピー平均を復元できる
        /// </summary>
        [Fact]
        public void 理論混合分布から単一コピー平均を復元できる()
        {
            const double l_真のλ = 30.0D;
            const double l_真の誤り平均 = 3.0D;
            const double l_真の誤り混合比 = 0.35D;
            var l_コピー数別混合比 = Get_コピー数別混合比_単一コピー優勢();

            var l_ヒストグラム = Get_理論ヒストグラム(
                l_真のλ, l_真の誤り混合比, l_真の誤り平均, l_コピー数別混合比, p_上限: 300, p_総数: 1_000_000L);

            var l_結果 = KmerSpectrumMixtureModel.Get_解析結果(l_ヒストグラム);

            Assert.NotNull(l_結果);
            Assert.InRange(l_結果!.A_単一コピー平均, l_真のλ - 2, l_真のλ + 2);
            // カットオフは誤り成分側、信頼下限はカットオフ以上、どちらも単一コピー峰 (30) 未満のはず
            Assert.InRange((double)l_結果.A_カットオフ, 1, l_真のλ);
            Assert.True(l_結果.A_信頼下限 >= l_結果.A_カットオフ);
        }

        /// <summary>
        /// 低カバレッジ気味 (単一コピー平均が誤り成分の平均に近い) でも、
        /// 谷が視認できるかどうかに関わらずモデルが分離できることを確かめる
        /// </summary>
        [Fact]
        public void 低カバレッジでも成分を分離できる()
        {
            const double l_真のλ = 12.0D;
            const double l_真の誤り平均 = 2.0D;
            const double l_真の誤り混合比 = 0.5D;
            var l_コピー数別混合比 = Get_コピー数別混合比_単一コピー優勢();

            var l_ヒストグラム = Get_理論ヒストグラム(
                l_真のλ, l_真の誤り混合比, l_真の誤り平均, l_コピー数別混合比, p_上限: 150, p_総数: 500_000L);

            var l_結果 = KmerSpectrumMixtureModel.Get_解析結果(l_ヒストグラム);

            Assert.NotNull(l_結果);
            Assert.InRange(l_結果!.A_単一コピー平均, l_真のλ - 2, l_真のλ + 2);
        }
    }
}
