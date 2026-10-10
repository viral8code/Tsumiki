using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 向きをそのままにした k-mer のキーを転がして作る処理の検証
    /// </summary>
    public class 順鎖Kmer転がしTests
    {
        #region 公開メソッド

        /// <summary>
        /// 転がして作ったキーは、窓ごとに詰めたキーと同じで、A・C・G・T 以外を含む窓は飛ばす
        /// </summary>
        /// <param name="p_k長">k-mer の長さ</param>
        [Theory]
        [InlineData(21)]
        [InlineData(32)]
        [InlineData(33)]
        [InlineData(64)]
        [InlineData(65)]
        [InlineData(135)]
        public void Get_追加_窓ごとに詰めたキーと同じ(int p_k長)
        {
            var l_乱数 = new Random(p_k長);
            var l_配列 = new string([..Enumerable.Range(0, 600).Select(i => i % 97 == 50 ? 'N' : Consts.塩基文字[l_乱数.Next(4)])]);
            var l_窓 = new 順鎖Kmer転がし(p_k長);
            List<(int A_終端, KmerKey A_キー)> l_転がし = [];
            for (var i = 0; i < l_配列.Length; i++)
            {
                if (l_窓.Is成功_追加(l_配列[i], out var l_キー))
                {
                    l_転がし.Add((i + 1, l_キー.Get_複製()));
                }
            }

            List<(int A_終端, KmerKey A_キー)> l_詰め直し = [];
            for (var l_終端 = p_k長; l_終端 <= l_配列.Length; l_終端++)
            {
                var l_窓の配列 = l_配列.AsSpan(l_終端 - p_k長, p_k長);
                if (!l_窓の配列.Contains('N'))
                {
                    l_詰め直し.Add((l_終端, new KmerKey(l_窓の配列)));
                }
            }

            Assert.Equal(l_詰め直し, l_転がし);
        }

        #endregion
    }
}
