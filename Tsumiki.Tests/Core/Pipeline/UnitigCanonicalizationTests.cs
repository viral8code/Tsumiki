using Tsumiki.Commons;
using Tsumiki.Cores.Pipeline;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// unitig の向きと閉路の回転を、起点に左右されない形にそろえる処理の検証
    /// </summary>
    public class UnitigCanonicalizationTests
    {
        #region 定数

        /// <summary>
        /// 塩基配列 TTGACCAGTA
        /// </summary>
        private const string C_塩基配列_TTGACCAGTA = "TTGACCAGTA";

        /// <summary>
        /// 塩基配列 ACGTTGCAAGCTTAGGCATC
        /// </summary>
        private const string C_塩基配列_ACGTTGCAAGCTTAGGCATC = "ACGTTGCAAGCTTAGGCATC";

        /// <summary>
        /// 塩基配列 ACG
        /// </summary>
        private const string C_塩基配列_ACG = "ACG";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 配列とその逆相補は同じ向きにそろうことを確かめる
        /// </summary>
        [Fact]
        public void Get_正規向き_逆相補どうしは同じになる()
        {
            const string l_配列 = C_塩基配列_TTGACCAGTA;
            Assert.Equal(AssemblyPipeline.Get_正規向き(l_配列), AssemblyPipeline.Get_正規向き(Util.V_逆相補(l_配列)));
        }

        /// <summary>
        /// 同じ閉路をどの k-mer から、どちらの向きで回っても同じ配列になることを確かめる
        /// </summary>
        [Fact]
        public void Get_正規回転_起点と向きに左右されない()
        {
            const int l_k長 = 5;
            const string l_本体 = C_塩基配列_ACGTTGCAAGCTTAGGCATC;
            var l_期待 = AssemblyPipeline.Get_正規回転(Get_一周配列(l_本体, 0, l_k長), l_k長);
            var l_逆本体 = Util.V_逆相補(l_本体);
            for (var i = 0; i < l_本体.Length; i++)
            {
                Assert.Equal(l_期待, AssemblyPipeline.Get_正規回転(Get_一周配列(l_本体, i, l_k長), l_k長));
                Assert.Equal(l_期待, AssemblyPipeline.Get_正規回転(Get_一周配列(l_逆本体, i, l_k長), l_k長));
            }

            Assert.Equal(l_本体.Length + l_k長 - 1, l_期待.Length);
        }

        /// <summary>
        /// k より短い閉路でも範囲外を読まずに回転をそろえることを確かめる
        /// </summary>
        [Fact]
        public void Get_正規回転_kより短い閉路でも扱える()
        {
            const int l_k長 = 7;
            const string l_本体 = C_塩基配列_ACG;
            var l_1周 = l_本体 + string.Concat(Enumerable.Repeat(l_本体, 3))[..(l_k長 - 1)];
            var l_結果 = AssemblyPipeline.Get_正規回転(l_1周, l_k長);
            Assert.Equal(l_本体.Length + l_k長 - 1, l_結果.Length);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 指定した開始位置から一周した閉路の配列
        /// </summary>
        /// <param name="p_本体"></param>
        /// <param name="p_開始"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        private static string Get_一周配列(string p_本体, int p_開始, int p_k長)
        {
            var l_回転 = p_本体[p_開始..] + p_本体[..p_開始];
            return l_回転 + l_回転[..(p_k長 - 1)];
        }

        #endregion
    }
}
