using Tsumiki.Commons;
using Tsumiki.Cores.Pipeline;
using Tsumiki.Core;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// multi-k で試す k の一覧の決め方
    /// </summary>
    /// <remarks>
    /// 実測では最適な k が範囲の両端に現れている (反復が少ない検体では上限の k=63、反復が 11% を占める検体では下限側の k=31) <br/>
    /// したがって候補は片側に寄せず、上限とその半分あたりの両方を含んでいる必要がある
    /// </remarks>
    public class MultiKAssemblerTests
    {
        #region 公開メソッド

        /// <summary>
        /// 一般的な Illumina リードでは、実測で最適だった範囲の両端を候補が含むことを確かめる
        /// </summary>
        [Fact]
        public void V_候補一覧_一般的なIlluminaリードでは最適だった範囲の両端を含む()
        {
            var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), p_リード長: 150);

            // 実測で最適だった 63 と 31 付近の両方が射程に入っていること
            Assert.True(l_候補[0] <= 33, $"lower end was {l_候補[0]}, too high to reach the repeat-rich optimum");
            Assert.Contains(l_候補, l_k長 => l_k長 is >= 55 and <= 71);
            // リード長に近い k まで届いていること
            // カバレッジが十分あれば
            // そちらのほうが良い場合があり、試さないと分からない
            Assert.True(l_候補[^1] >= 120, $"upper end was {l_候補[^1]}, too low to reach the high-k regime");
            Assert.Equal(Consts.マルチkで試す個数, l_候補.Count);
        }

        /// <summary>
        /// 候補一覧が重複なく昇順に並んでいることを確かめる
        /// </summary>
        [Fact]
        public void V_候補一覧_昇順で重複なく並ぶ()
        {
            foreach (var l_リード長 in new[] { 50, 75, 100, 150, 250, 300 })
            {
                var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), l_リード長);
                Assert.Equal(l_候補.OrderBy(x => x).ToList(), l_候補);
                Assert.Equal(l_候補.Distinct().Count(), l_候補.Count);
            }
        }

        /// <summary>
        /// k が偶数だと k-mer 自身がその逆相補と一致しうる (回文) ため、正規形が縮退して隣接判定が壊れる
        /// </summary>
        /// <remarks>
        /// どの候補も奇数であること
        /// </remarks>
        [Fact]
        public void V_候補一覧_奇数のみを含む()
        {
            for (var l_リード長 = 40; l_リード長 <= 300; l_リード長++)
            {
                foreach (var l_k長 in MultiKAssembler.Get_k候補一覧(Get_引数(), l_リード長))
                {
                    Assert.True(l_k長 % 2 == 1, $"k={l_k長} for read length {l_リード長} is even");
                }
            }
        }

        /// <summary>
        /// 候補はリード長より短くなければならない
        /// </summary>
        /// <remarks>
        /// そうでないとその k では k-mer が 1 つも取れない
        /// </remarks>
        [Fact]
        public void V_候補一覧_リード長より短い()
        {
            for (var l_リード長 = 40; l_リード長 <= 300; l_リード長++)
            {
                foreach (var l_k長 in MultiKAssembler.Get_k候補一覧(Get_引数(), l_リード長))
                {
                    Assert.True(l_k長 < l_リード長, $"k={l_k長} is not shorter than read length {l_リード長}");
                }
            }
        }

        /// <summary>
        /// -k に一覧が指定された場合は、それをそのまま使うこと
        /// </summary>
        /// <remarks>
        /// 利用者が試す値を選んだのに、自動の刻みで置き換えてはいけない
        /// </remarks>
        [Fact]
        public void V_候補一覧_k長一覧が指定されていればそのまま使う()
        {
            var l_引数 = new Parameters();
            l_引数.Set_k長一覧([95, 31, 63]);

            var l_候補 = MultiKAssembler.Get_k候補一覧(l_引数, p_リード長: 250);

            Assert.Equal([31, 63, 95], l_候補);
        }

        /// <summary>
        /// -k を 1 つだけ指定した場合は、その 1 つだけを候補にすること
        /// </summary>
        [Fact]
        public void V_候補一覧_k長を1つだけ指定すればそれだけを候補にする()
        {
            var l_引数 = new Parameters();
            l_引数.Set_k長一覧([41]);

            Assert.Equal([41], MultiKAssembler.Get_k候補一覧(l_引数, p_リード長: 250));
        }

        /// <summary>
        /// 予測 k-mer カバレッジは、1 リードから取れる k-mer の本数の比で縮むこと
        /// </summary>
        /// <remarks>
        /// カバレッジの薄いデータで高い k を試すのは時間を捨てるだけになる
        /// </remarks>
        [Fact]
        public void V_予測カバレッジ_1リードから取れるkmer数の比で縮む()
        {
            // リード長 150、k=31 で 25 x
            // k=135 なら 1 リードあたり 120 本から
            // 16 本へ減るので、25 * 16 / 120 = 3.33 x
            var l_予測 = MultiKAssembler.Get_予測kmerカバレッジ(p_直前の基準値: 25.0D, p_直前のk長: 31, p_次のk長: 135, p_リード長: 150);

            Assert.Equal(25.0D * 16D / 120D, l_予測, 3);
            Assert.True(l_予測 < Consts.マルチkの最小kmerカバレッジ);
        }

        /// <summary>
        /// k がリード長を超えると k-mer が 1 本も取れないので 0 になること
        /// </summary>
        /// <remarks>
        /// k = リード長 のときは 1 本だけ取れるので 0 にはならない
        /// </remarks>
        [Fact]
        public void V_予測カバレッジ_kがリード長を超えると0になる()
        {
            Assert.Equal(0D, MultiKAssembler.Get_予測kmerカバレッジ(p_直前の基準値: 25.0D, p_直前のk長: 31, p_次のk長: 151, p_リード長: 150));
            Assert.True(MultiKAssembler.Get_予測kmerカバレッジ(p_直前の基準値: 25.0D, p_直前のk長: 31, p_次のk長: 150, p_リード長: 150) > 0D);
        }

        /// <summary>
        /// リード長が分からない場合でも一覧が作れること (既定値を上限に使う)
        /// </summary>
        [Fact]
        public void V_候補一覧_リード長が不明なら既定値を使う()
        {
            var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), p_リード長: null);

            Assert.NotEmpty(l_候補);
            Assert.Equal(31, l_候補[^1]);
        }

        /// <summary>
        /// 上限と下限が重なるほど短いリードでは、候補が 1 つに縮退すること (無理に複数試しても意味がない)
        /// </summary>
        [Fact]
        public void V_候補一覧_非常に短いリードでは1つに縮退する()
        {
            var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), p_リード長: 40);

            Assert.NotEmpty(l_候補);
            Assert.All(l_候補, l_k長 => Assert.True(l_k長 < 40));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 既定のままの実行時引数
        /// </summary>
        /// <returns>実行時引数</returns>
        private static Parameters Get_引数() => new();

        #endregion

    }
}
