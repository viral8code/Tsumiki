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
            var l_組 = UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 1, 91 }, { 51, 16 } }, new ulong[,] { { 0, 10 }, { 0, 0 } }, C_優勢閾値, C_最小証拠数);

            Assert.NotNull(l_組);
            Assert.Equal([(10, 40), (20, 30)], l_組);
        }

        /// <summary>
        /// まっすぐの対応も同じように返す
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_まっすぐの対応を返す()
        {
            var l_組 = UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 61, 9 }, { 17, 31 } }, new ulong[,] { { 18, 0 }, { 0, 8 } }, C_優勢閾値, C_最小証拠数);

            Assert.Equal([(10, 30), (20, 40)], l_組);
        }

        /// <summary>
        /// 経路とペアが逆の対応を指せば解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_証拠が逆の対応を指せば解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 7, 0 }, { 24, 0 } }, new ulong[,] { { 15, 0 }, { 0, 0 } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 合算しても一方に寄らなければ解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_合算で寄らなければ解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 38, 30 }, { 43, 100 } }, new ulong[,] { { 1, 0 }, { 0, 1 } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 片方の入口がどの出口にも同じほど繋がるなら、合算で寄っていても解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_組の行で僅差なら解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 200, 0 }, { 20, 20 } }, new ulong[,] { { 10, 0 }, { 0, 1 } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 選ばない側の組を経路とペアの両方が支持するなら、1 対 1 で表せない形とみなして解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_選ばない組を両方の証拠が支持すれば解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 0, 11 }, { 9, 6 } }, new ulong[,] { { 0, 18 }, { 27, 5 } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 経路かペアの片方にしか証拠が無ければ、合算では解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_片方の証拠しか無ければ解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 1, 9 }, { 5, 2 } }, new ulong[,] { { 0, 0 }, { 0, 0 } }, C_優勢閾値, C_最小証拠数));
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 0, 0 }, { 0, 0 } }, new ulong[,] { { 20, 0 }, { 0, 20 } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 証拠が少なければ解かない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_証拠が少なければ解かない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20], [30, 40], new ulong[,] { { 4, 0 }, { 0, 4 } }, new ulong[,] { { 1, 0 }, { 0, 0 } }, C_優勢閾値, C_最小証拠数));
        }

        /// <summary>
        /// 入口か出口が 2 つでなければ扱わない
        /// </summary>
        [Fact]
        public void Get_二対二の合算で決着した組_二対二でなければ扱わない()
        {
            Assert.Null(UnitigGraph.Get_二対二の合算で決着した組([10, 20, 50], [30, 40], new ulong[,] { { 50, 0 }, { 0, 50 }, { 0, 0 } }, new ulong[,] { { 0, 0 }, { 0, 0 }, { 0, 0 } }, C_優勢閾値, C_最小証拠数));
        }

        #endregion
    }
}
