using Tsumiki.Commons;
using Tsumiki.Cores.Pipeline;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 間引いて余ったリードで端に触れるものを選ぶ処理 (FinalAssemblyPipeline) の検証
    /// </summary>
    public class FinalAssemblyPipelineEndTests
    {
        #region 定数

        /// <summary>
        /// 配列の長さ
        /// </summary>
        private const int C_配列長 = 1000;

        /// <summary>
        /// 切り出すリードの長さ
        /// </summary>
        private const int C_リードの長さ = 40;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 先頭 300 bp と末尾 300 bp に触れるリードだけが残り、中央だけのリードは残らない
        /// </summary>
        [Fact]
        public void Get_端に触れるリード_先頭と末尾に触れるものだけ残る()
        {
            var l_配列 = Get_ランダムな塩基(new Random(7), C_配列長);
            var l_先頭リード = l_配列.Substring(100, C_リードの長さ);
            var l_末尾リード = Util.V_逆相補_曖昧塩基あり(l_配列.Substring(900, C_リードの長さ));
            var l_中央リード = l_配列.Substring(450, C_リードの長さ);

            var l_種集合 = FinalAssemblyPipeline.Get_端の種集合(new[] { l_配列 });
            var l_リード群 = new List<string> { l_先頭リード, l_中央リード, l_末尾リード };
            var l_結果 = FinalAssemblyPipeline.Get_端に触れるリード(l_種集合, l_リード群);

            var l_期待 = new List<string> { l_先頭リード, l_末尾リード };
            Assert.Equal(l_期待, l_結果);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 指定した乱数で、A/C/G/T からなる塩基列を作る
        /// </summary>
        /// <param name="p_乱数">乱数</param>
        /// <param name="p_長さ">長さ</param>
        /// <returns>塩基列</returns>
        private static string Get_ランダムな塩基(Random p_乱数, int p_長さ)
        {
            const string l_塩基 = Consts.塩基文字;
            var l_配列 = new char[p_長さ];
            for (var i = 0; i < p_長さ; i++)
            {
                l_配列[i] = l_塩基[p_乱数.Next(4)];
            }

            return new string(l_配列);
        }

        #endregion
    }
}
