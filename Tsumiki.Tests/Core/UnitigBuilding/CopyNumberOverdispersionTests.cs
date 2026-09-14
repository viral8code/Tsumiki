using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 過分散なカバレッジでのコピー数判定の検証
    /// </summary>
    public class CopyNumberOverdispersionTests
    {
        #region 公開メソッド

        /// <summary>
        /// 単一コピー集団のばらつきが大きいと、基準値の 1.5 倍を少し超えるだけの unitig は単一コピーのまま扱う
        /// </summary>
        [Fact]
        public void V_過分散なら少し高いだけのunitigを反復と判定しない()
        {
            Dictionary<int, double> l_カバレッジ = new() { [1] = 50D, [2] = 60D, [3] = 70D, [4] = 80D, [5] = 90D, [6] = 100D, [7] = 110D, [8] = 120D, [9] = 130D, [10] = 140D, [11] = 155D, [12] = 300D };
            var l_長さ一覧 = l_カバレッジ.Keys.ToDictionary(x => x, _ => 600);

            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧, p_基準の出所: コピー数基準の出所.Weighted);

            Assert.Equal(1, l_結果.A_コピー数[11]);
            Assert.Equal(3, l_結果.A_コピー数[12]);
        }

        /// <summary>
        /// ばらつきが小さければ、従来どおり基準値の 1.5 倍から多コピーとみなす
        /// </summary>
        [Fact]
        public void V_ばらつきが小さければ1_5倍から多コピーとみなす()
        {
            Dictionary<int, double> l_カバレッジ = new() { [1] = 97D, [2] = 98D, [3] = 99D, [4] = 100D, [5] = 100D, [6] = 100D, [7] = 101D, [8] = 102D, [9] = 103D, [10] = 104D, [11] = 155D, [12] = 300D };
            var l_長さ一覧 = l_カバレッジ.Keys.ToDictionary(x => x, _ => 600);

            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ一覧, p_基準の出所: コピー数基準の出所.Weighted);

            Assert.Equal(2, l_結果.A_コピー数[11]);
            Assert.Equal(3, l_結果.A_コピー数[12]);
        }

        #endregion
    }
}
