using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 複数のアセンブリ候補から1つを選ぶ規則の検証。
    ///
    /// 単一のスコアに畳む方式を採らなかった経緯がそのままここの主題になる。
    /// 「NG50 × 完全性 × 正確性」で選ぶ実装を試したところ、反復配列を飛ばして
    /// 中間を落としたキメラ(完全性 0.675、NG50 16,300)が、正直に途切れた答え
    /// (完全性 0.933、NG50 8,300)より高い点になった。連続性の利得が
    /// 完全性の損失を上回るためで、指数を調整して隠すのではなく
    /// 「まず完全性で足切りしてから連続性を見る」という順序にした。
    /// </summary>
    public class AssemblySelectorTests
    {
        private static (アセンブリ実行結果, アセンブリ評価) Get_候補(
            int p_k長, long p_NG50, double p_完全性, double p_正確性 = 1.0)
        {
            // 期待延べ数を固定し、そこから逆算して欠損・過剰を決める。
            const long l_期待延べ数 = 1_000_000;
            var l_実行結果 = new アセンブリ実行結果(
                p_k長, $"k{p_k長}_unitigs.fasta", $"k{p_k長}_contigs.fasta",
                $"k{p_k長}_scaffolds.fasta", 2, 20.0);
            var l_評価 = new アセンブリ評価(
                A_期待延べ数: l_期待延べ数,
                A_欠損延べ数: (long)(l_期待延べ数 * (1 - p_完全性)),
                A_過剰延べ数: (long)(l_期待延べ数 * (1 - p_正確性)),
                A_総延長: 5_000_000,
                A_本数: 100,
                A_NG50: p_NG50);
            return (l_実行結果, l_評価);
        }

        [Fact]
        public void Select_EmptyCandidates_ReturnsNull()
        {
            Assert.Null(AssemblySelector.Get_最良([]));
        }

        [Fact]
        public void Select_SingleCandidate_ReturnsIt()
        {
            var 候補 = Get_候補(31, 50_000, 0.97);

            var 選択 = AssemblySelector.Get_最良([候補]);

            Assert.NotNull(選択);
            Assert.Equal(31, 選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 完全性が同程度なら、連続性が高いほうを採る。
        /// Axy の実データ(どの k でも配列は落ちず、k=63 が最も繋がる)がこの形。
        /// </summary>
        [Fact]
        public void Select_SimilarCompleteness_PicksTheMostContiguous()
        {
            var 選択 = AssemblySelector.Get_最良([
                Get_候補(31, 63_058, 0.980),
                Get_候補(45, 151_085, 0.981),
                Get_候補(63, 175_674, 0.979),
            ]);

            Assert.NotNull(選択);
            Assert.Equal(63, 選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// R. sphaeroides の実測値そのもの。完全性の差(97.1% と 93.2%)は
        /// 許容差に収まるため両方が残り、連続性で k=31 が選ばれる。
        /// </summary>
        [Fact]
        public void Select_RealWorldSpread_PicksTheKnownBestK()
        {
            var 選択 = AssemblySelector.Get_最良([
                Get_候補(31, 53_893, 0.9708),
                Get_候補(45, 36_387, 0.9684),
                Get_候補(55, 19_347, 0.9528),
                Get_候補(63, 15_750, 0.9316),
            ]);

            Assert.NotNull(選択);
            Assert.Equal(31, 選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// これが二段階にした理由。連続性では圧倒的に上でも、
        /// 完全性が許容差を超えて落ちている候補は採らない。
        /// 掛け算で選んでいたらこちらが選ばれていた。
        /// </summary>
        [Fact]
        public void Select_MuchMoreContiguousButIncomplete_IsRejected()
        {
            var 選択 = AssemblySelector.Get_最良([
                Get_候補(31, 8_300, 0.933),
                Get_候補(63, 16_300, 0.675),
            ]);

            Assert.NotNull(選択);
            Assert.Equal(31, 選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 同じ領域を重複して出している候補も、連続性が高くても採らない。
        /// </summary>
        [Fact]
        public void Select_MoreContiguousButDuplicated_IsRejected()
        {
            var 選択 = AssemblySelector.Get_最良([
                Get_候補(31, 50_000, 0.97, p_正確性: 0.99),
                Get_候補(63, 90_000, 0.97, p_正確性: 0.60),
            ]);

            Assert.NotNull(選択);
            Assert.Equal(31, 選択.Value.A_実行結果.A_k長);
        }

        /// <summary>
        /// 完全性の差がちょうど許容差の内側なら足切りされないこと
        /// (境界の振る舞いを固定しておく)。
        /// </summary>
        [Fact]
        public void Select_CompletenessGapJustInsideTheTolerance_KeepsTheCandidate()
        {
            var 選択 = AssemblySelector.Get_最良([
                Get_候補(31, 10_000, 0.99),
                Get_候補(63, 90_000, 0.99 - AssemblySelector.完全性の許容差 + 0.001),
            ]);

            Assert.NotNull(選択);
            Assert.Equal(63, 選択.Value.A_実行結果.A_k長);
        }

        [Fact]
        public void Select_CompletenessGapJustOutsideTheTolerance_RejectsTheCandidate()
        {
            var 選択 = AssemblySelector.Get_最良([
                Get_候補(31, 10_000, 0.99),
                Get_候補(63, 90_000, 0.99 - AssemblySelector.完全性の許容差 - 0.001),
            ]);

            Assert.NotNull(選択);
            Assert.Equal(31, 選択.Value.A_実行結果.A_k長);
        }
    }
}
