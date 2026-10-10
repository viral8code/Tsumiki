using Tsumiki.Commons;
using Tsumiki.Cores.Pipeline;
using Tsumiki.Core;
using Tsumiki.Models.Correction;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// multi-k で試す k の一覧の決め方
    /// </summary>
    public class MultiKAssemblerTests
    {
        #region 公開メソッド

        /// <summary>
        /// 最後の k と引き継ぎ無効時には次段用処理を要求しないことを検証する
        /// </summary>
        [Fact]
        public void V_最後のkではSuperRead用の引き継ぎを要求しない()
        {
            var l_引数 = new Parameters();
            Assert.False(MultiKAssembler.Is次段引き継ぎ必要(l_引数, [49], 49));
            Assert.True(MultiKAssembler.Is次段引き継ぎ必要(l_引数, [21, 49], 21));
            Assert.False(MultiKAssembler.Is次段引き継ぎ必要(l_引数, [21, 49], 49));
            l_引数.A_Is引き継ぎ = false;
            Assert.False(MultiKAssembler.Is次段引き継ぎ必要(l_引数, [21, 49], 21));
        }

        /// <summary>
        /// 一般的な Illumina リードでは、実測で最適だった範囲の両端を候補が含むことを確かめる
        /// </summary>
        [Fact]
        public void V_候補一覧_一般的なIlluminaリードでは最適だった範囲の両端を含む()
        {
            var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), p_リード長: 150);
            Assert.True(l_候補[0] <= 33, $"lower end was {l_候補[0]}, too high to reach the repeat-rich optimum");
            Assert.Contains(l_候補, l_k長 => l_k長 is >= 55 and <= 71);
            Assert.True(l_候補[^1] >= 120, $"upper end was {l_候補[^1]}, too low to reach the high-k regime");
            Assert.Equal(Consts.マルチkで試す個数, l_候補.Count);
        }

        /// <summary>
        /// 候補一覧が重複なく昇順に並んでいることを確かめる
        /// </summary>
        [Fact]
        public void V_候補一覧_昇順で重複なく並ぶ()
        {
            foreach (var l_リード長 in new[]
            {
                50,
                75,
                100,
                150,
                250,
                300
            }

            )
            {
                var l_候補 = MultiKAssembler.Get_k候補一覧(Get_引数(), l_リード長);
                Assert.Equal(l_候補.OrderBy(x => x).ToList(), l_候補);
                Assert.Equal(l_候補.Distinct().Count(), l_候補.Count);
            }
        }

        /// <summary>
        /// k が偶数だと k-mer 自身がその逆相補と一致しうる (回文) ため、正規形が縮退して隣接判定が壊れる
        /// </summary>
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
        [Fact]
        public void V_候補一覧_k長一覧が指定されていればそのまま使う()
        {
            var l_引数 = new Parameters();
            l_引数.Set_k長一覧([95, 31, 63]);
            var l_候補 = MultiKAssembler.Get_k候補一覧(l_引数, p_リード長: 250);
            Assert.Equal([31, 63, 95], l_候補);
        }

        /// <summary>
        /// エラー訂正で数えた区間があれば、予測カバレッジが下限を保てる最大の k を上限にして等比で刻むこと
        /// </summary>
        [Fact]
        public void V_候補一覧_訂正後のカバレッジから上限を決める()
        {
            var l_度数 = new 無誤り区間の度数(63)
            {
                A_単一コピー平均 = 40D
            };
            var l_信頼状況 = Enumerable.Repeat(true, 251 - 63 + 1).ToArray();
            l_度数.V_追加(l_信頼状況, false);
            l_度数.V_追加(l_信頼状況, true);
            var l_引数 = Get_引数();
            l_引数.A_無誤り区間の度数群.Add(l_度数);
            var l_候補 = MultiKAssembler.Get_k候補一覧(l_引数, p_リード長: 251);
            Assert.Equal(21, l_候補[0]);
            Assert.Equal(203, l_候補[^1]);
            Assert.Equal(Consts.マルチkで試す個数, l_候補.Count);
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
        [Fact]
        public void V_予測カバレッジ_1リードから取れるkmer数の比で縮む()
        {
            var l_予測 = MultiKAssembler.Get_予測kmerカバレッジ(p_直前の基準値: 25.0D, p_直前のk長: 31, p_次のk長: 135, p_リード長: 150);
            Assert.Equal(25.0D * 16D / 120D, l_予測, 3);
            Assert.True(l_予測 < Consts.マルチkの最小kmerカバレッジ);
        }

        /// <summary>
        /// k がリード長を超えると k-mer が 1 本も取れないので 0 になること
        /// </summary>
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
