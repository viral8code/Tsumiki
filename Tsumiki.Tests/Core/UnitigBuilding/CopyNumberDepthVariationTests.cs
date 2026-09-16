using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Models.UnitigBuilding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 被覆のばらつきだけでコピー数を増やさないことの検証
    /// </summary>
    public class CopyNumberDepthVariationTests
    {
        #region 定数

        /// <summary>
        /// 単一コピー unitig の本数
        /// </summary>
        private const int 単一コピーの本数 = 40;

        /// <summary>
        /// 単一コピー unitig の長さ
        /// </summary>
        private const int unitig長 = 2_000;

        /// <summary>
        /// 単一コピー unitig の最小カバレッジ
        /// </summary>
        private const double 最小カバレッジ = 25D;

        /// <summary>
        /// 単一コピー unitig の最大カバレッジ
        /// </summary>
        /// <remarks>
        /// 最小の 4 倍
        /// </remarks>
        private const double 最大カバレッジ = 100D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 単一コピーの被覆が 4 倍まで変動しても、どれも 2 コピー以上にしない
        /// </summary>
        [Fact]
        public void V_被覆が4倍変動しても単一コピーのまま()
        {
            var (l_カバレッジ, l_長さ) = Get_単一コピー群();

            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ, p_基準の出所: コピー数基準の出所.Weighted);

            Assert.All(l_カバレッジ.Keys, x => Assert.Equal(1, l_結果.A_コピー数.GetValueOrDefault(x, 1)));
        }

        /// <summary>
        /// 被覆が 4 倍変動する中でも、基準の 3 倍を超える unitig は多コピーと判定する
        /// </summary>
        /// <remarks>
        /// 上の検証が「何でも 1 にする」ことで通っていないことを確かめる
        /// </remarks>
        [Fact]
        public void V_被覆が4倍変動しても明確な反復は多コピーにする()
        {
            var (l_カバレッジ, l_長さ) = Get_単一コピー群();
            var l_反復ID = 単一コピーの本数 + 1;
            l_カバレッジ[l_反復ID] = 3.5D * l_カバレッジ.Values.Order().ElementAt(単一コピーの本数 / 2);
            l_長さ[l_反復ID] = unitig長;

            var l_結果 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_長さ, p_基準の出所: コピー数基準の出所.Weighted);

            Assert.True(l_結果.A_コピー数.GetValueOrDefault(l_反復ID, 1) >= 2);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// カバレッジが最小から最大まで等間隔に並ぶ単一コピー unitig の集まり
        /// </summary>
        /// <returns>unitig ID ごとのカバレッジと長さ</returns>
        private static (Dictionary<int, double> A_カバレッジ, Dictionary<int, int> A_長さ) Get_単一コピー群()
        {
            Dictionary<int, double> l_カバレッジ = [];
            Dictionary<int, int> l_長さ = [];
            for (var i = 0; i < 単一コピーの本数; i++)
            {
                l_カバレッジ[i + 1] = 最小カバレッジ + ((最大カバレッジ - 最小カバレッジ) * i / (単一コピーの本数 - 1));
                l_長さ[i + 1] = unitig長;
            }
            return (l_カバレッジ, l_長さ);
        }

        #endregion
    }
}
