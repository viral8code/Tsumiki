using Tsumiki.Cores.UnitigBuilding;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 入口 2 つ・出口 2 つの反復を、経路とペアの合算で対応づける判定の検証
    /// </summary>
    public class TwoByTwoRepeatPairingTests
    {
        #region 定数

        /// <summary>
        /// 優勢とみなす割合
        /// </summary>
        private const decimal C_優勢閾値 = 0.8M;

        /// <summary>
        /// 合算で要る証拠数
        /// </summary>
        private const ulong C_最小証拠数 = 10UL;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 経路の 1 行が僅差でも、ペアが同じ交差の対応を指し、合算ではっきり寄れば解く
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_同じ向きを指す証拠を合わせて交差の対応を返す()
        {
            var l_組 = UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 1UL, 91UL }, { 51UL, 16UL } }, new ulong[,] { { 0UL, 10UL }, { 0UL, 0UL } }, C_優勢閾値, C_最小証拠数);
            Assert.NotNull(l_組);
            Assert.Equal([(10, 40), (20, 30)], l_組);
        }

        /// <summary>
        /// まっすぐの対応も同じように返す
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_まっすぐの対応を返す()
        {
            var l_組 = UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 61UL, 9UL }, { 17UL, 31UL } }, new ulong[,] { { 18UL, 0UL }, { 0UL, 8UL } }, C_優勢閾値, C_最小証拠数);
            Assert.Equal([(10, 30), (20, 40)], l_組);
        }

        /// <summary>
        /// 経路とペアが逆の対応を指せば解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_証拠が逆の対応を指せば解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 7UL, 0UL }, { 24UL, 0UL } }, new ulong[,] { { 15UL, 0UL }, { 0UL, 0UL } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 合算しても一方に寄らなければ解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_合算で寄らなければ解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 38UL, 30UL }, { 43UL, 100UL } }, new ulong[,] { { 1UL, 0UL }, { 0UL, 1UL } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 片方の入口がどの出口にも同じほど繋がるなら、合算で寄っていても解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_組の行で僅差なら解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 200UL, 0UL }, { 20UL, 20UL } }, new ulong[,] { { 10UL, 0UL }, { 0UL, 1UL } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 選ばない側の組を経路とペアの両方が支持するなら、1 対 1 で表せない形とみなして解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_選ばない組を両方の証拠が支持すれば解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 0UL, 11UL }, { 9UL, 6UL } }, new ulong[,] { { 0UL, 18UL }, { 27UL, 5UL } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 経路かペアの片方にしか証拠が無ければ、合算では解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_片方の証拠しか無ければ解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 1UL, 9UL }, { 5UL, 2UL } }, new ulong[,] { { 0UL, 0UL }, { 0UL, 0UL } }, C_優勢閾値, C_最小証拠数));
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 0UL, 0UL }, { 0UL, 0UL } }, new ulong[,] { { 20UL, 0UL }, { 0UL, 20UL } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 証拠が少なければ解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_証拠が少なければ解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 4UL, 0UL }, { 0UL, 4UL } }, new ulong[,] { { 1UL, 0UL }, { 0UL, 0UL } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 入口か出口が 2 つでなければ扱わない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_二対二でなければ扱わない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20, 50], [30, 40], new ulong[,] { { 50UL, 0UL }, { 0UL, 50UL }, { 0UL, 0UL } }, new ulong[,] { { 0UL, 0UL }, { 0UL, 0UL }, { 0UL, 0UL } }, C_優勢閾値, C_最小証拠数));
        }

        #endregion
    }
}
