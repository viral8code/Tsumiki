using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// multi-k で試す k の一覧の決め方。
    ///
    /// 実測では最適な k が範囲の両端に現れている(反復が少ない検体では上限の
    /// k=63、反復が 11% を占める検体では下限側の k=31)。したがって候補は
    /// 片側に寄せず、上限とその半分あたりの両方を含んでいる必要がある。
    /// </summary>
    public class MultiKAssemblerTests
    {
        private static Parameters Get_引数() => new();

        [Fact]
        public void CandidateList_ForTypicalIlluminaReads_SpansTheRangeWhereOptimaWereObserved()
        {
            var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), p_リード長: 150);

            // 実測で最適だった 63(Achromobacter)と 31 付近(R. sphaeroides)の
            // 両方が射程に入っていること。
            Assert.Equal(Consts.自動k長の上限, l_候補[^1]);
            Assert.True(l_候補[0] <= 33, $"lower end was {l_候補[0]}, too high to reach the repeat-rich optimum");
            Assert.Equal(Consts.マルチkで試す個数, l_候補.Count);
        }

        [Fact]
        public void CandidateList_IsSortedAscendingWithoutDuplicates()
        {
            foreach (var l_リード長 in new[] { 50, 75, 100, 150, 250, 300 })
            {
                var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), l_リード長);
                Assert.Equal(l_候補.OrderBy(x => x).ToList(), l_候補);
                Assert.Equal(l_候補.Distinct().Count(), l_候補.Count);
            }
        }

        /// <summary>
        /// k が偶数だと k-mer 自身がその逆相補と一致しうる(回文)ため、
        /// 正規形が縮退して隣接判定が壊れる。どの候補も奇数であること。
        /// </summary>
        [Fact]
        public void CandidateList_ContainsOnlyOddValues()
        {
            for (var l_リード長 = 40; l_リード長 <= 300; l_リード長++)
            {
                foreach (var l_k in MultiKAssembler.Get_k候補一覧(Get_引数(), l_リード長))
                {
                    Assert.True(l_k % 2 == 1, $"k={l_k} for read length {l_リード長} is even");
                }
            }
        }

        /// <summary>
        /// 候補はリード長より短くなければならない。そうでないと
        /// その k では k-mer が1つも取れない。
        /// </summary>
        [Fact]
        public void CandidateList_StaysBelowTheReadLength()
        {
            for (var l_リード長 = 40; l_リード長 <= 300; l_リード長++)
            {
                foreach (var l_k in MultiKAssembler.Get_k候補一覧(Get_引数(), l_リード長))
                {
                    Assert.True(l_k < l_リード長, $"k={l_k} is not shorter than read length {l_リード長}");
                }
            }
        }

        /// <summary>
        /// -k が明示指定されている場合は、その値を超えて試さないこと。
        /// 「これ以上は上げるな」という利用者の判断を multi-k が覆してはいけない。
        /// </summary>
        [Fact]
        public void CandidateList_WhenKmerLengthWasGivenExplicitly_NeverExceedsIt()
        {
            var l_引数 = new Parameters { A_k長 = 41 };

            var l_候補 = MultiKAssembler.Get_k候補一覧(l_引数, p_リード長: 250);

            Assert.All(l_候補, l_k => Assert.True(l_k <= 41, $"k={l_k} exceeded the explicit -k 41"));
            Assert.Equal(41, l_候補[^1]);
        }

        /// <summary>
        /// リード長が分からない場合でも一覧が作れること(既定値を上限に使う)。
        /// </summary>
        [Fact]
        public void CandidateList_WhenReadLengthIsUnknown_FallsBackToTheDefault()
        {
            var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), p_リード長: null);

            Assert.NotEmpty(l_候補);
            Assert.Equal(Consts.k長の既定値, l_候補[^1]);
        }

        /// <summary>
        /// 上限と下限が重なるほど短いリードでは、候補が1つに縮退すること
        /// (無理に複数試しても意味がない)。
        /// </summary>
        [Fact]
        public void CandidateList_ForVeryShortReads_CollapsesToASingleValue()
        {
            var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), p_リード長: 40);

            Assert.NotEmpty(l_候補);
            Assert.All(l_候補, l_k => Assert.True(l_k < 40));
        }
    }
}
