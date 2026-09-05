using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 複数の k のアセンブリを統合する処理の検証。
    ///
    /// 統合は誤った連結を持ち込みうる操作なので、繋ぐべきときに繋ぐことと
    /// 同じくらい、根拠が無いときに繋がないことを固定しておく必要がある。
    /// </summary>
    public class AssemblyMergerTests : IDisposable
    {
        private readonly string _tempDir;

        public AssemblyMergerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_merger_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private const int アンカーk長 = 31;

        // 既定では2つ以上の k による裏付けを求める。以下の多くのテストは
        // 証拠源が1つの状況を見たいので、明示的に1を渡している。
        // 照合そのものは専用のテストで確かめる。

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        private アセンブリ実行結果 WriteAssembly(string name, int kmerLength, params string[] sequences)
        {
            var path = Path.Combine(this._tempDir, name);
            using (var writer = new FastaWriter(path))
            {
                var id = 1;
                foreach (var seq in sequences)
                {
                    writer.V_書き込み($"NODE{id++}", seq);
                }
            }
            return new アセンブリ実行結果(kmerLength, path, path, null, 2, 20.0);
        }

        private static List<string> ReadSequences(string path)
        {
            List<string> result = [];
            using var reader = new FastaReader(path);
            while (reader.Get_続きがあるか())
            {
                result.Add(reader.Get_次の配列().A_配列);
            }
            return result;
        }

        /// <summary>
        /// 骨格が途切れている箇所を、別の k の配列が跨いでいる場合。
        /// 繋いだ結果が元のゲノムそのものに戻ること。
        /// </summary>
        [Fact]
        public void Merge_OtherKSpansABackboneJunction_JoinsThemBackIntoTheTruth()
        {
            var 左 = RandomSequence(5_000, seed: 601);
            var 中間 = RandomSequence(400, seed: 602);
            var 右 = RandomSequence(5_000, seed: 603);
            var truth = 左 + 中間 + 右;

            // 骨格は中間で切れている。
            var 骨格 = this.WriteAssembly("backbone.fasta", 63, 左, 右);
            // 別の k は切れ目を跨いでいる(両端に十分なアンカーを持つ)。
            var 他 = this.WriteAssembly("other.fasta", 31, truth);

            var 出力 = Path.Combine(this._tempDir, "merged.fasta");
            var 繋いだか = AssemblyMerger.V_統合(骨格, [骨格, 他], アンカーk長, 出力, p_必要な独立支持数: 1);

            Assert.True(繋いだか);
            var 結果 = ReadSequences(出力);
            _ = Assert.Single(結果);
            Assert.Equal(truth, 結果[0]);
        }

        /// <summary>
        /// 骨格側の片方が逆向きに出力されていても、向きを揃えて繋げること。
        /// </summary>
        [Fact]
        public void Merge_BackbonePieceIsReverseComplemented_StillJoinsCorrectly()
        {
            var 左 = RandomSequence(5_000, seed: 611);
            var 中間 = RandomSequence(300, seed: 612);
            var 右 = RandomSequence(5_000, seed: 613);
            var truth = 左 + 中間 + 右;

            var 骨格 = this.WriteAssembly("backbone_rc.fasta", 63, 左, Util.V_逆相補(右));
            var 他 = this.WriteAssembly("other_rc.fasta", 31, truth);

            var 出力 = Path.Combine(this._tempDir, "merged_rc.fasta");
            var 繋いだか = AssemblyMerger.V_統合(骨格, [骨格, 他], アンカーk長, 出力, p_必要な独立支持数: 1);

            Assert.True(繋いだか);
            var 結果 = ReadSequences(出力);
            _ = Assert.Single(結果);
            Assert.True(結果[0] == truth || 結果[0] == Util.V_逆相補(truth),
                "merged sequence should be the truth in one orientation or the other");
        }

        /// <summary>
        /// 跨いでいる配列が無ければ何もしないこと。
        /// 根拠が無いのに繋ぐのが最も避けたい失敗。
        /// </summary>
        [Fact]
        public void Merge_NoOtherAssemblySpansAnything_DoesNothing()
        {
            var 左 = RandomSequence(5_000, seed: 621);
            var 右 = RandomSequence(5_000, seed: 622);

            var 骨格 = this.WriteAssembly("backbone_none.fasta", 63, 左, 右);
            // 別の k も同じところで切れている。
            var 他 = this.WriteAssembly("other_none.fasta", 31, 左, 右);

            var 出力 = Path.Combine(this._tempDir, "merged_none.fasta");
            var 繋いだか = AssemblyMerger.V_統合(骨格, [骨格, 他], アンカーk長, 出力, p_必要な独立支持数: 1);

            Assert.False(繋いだか);
            Assert.False(File.Exists(出力));
        }

        /// <summary>
        /// 反復配列のせいで行き先が2つある場合は繋がないこと。
        /// 片方を選ぶ根拠が無く、選べば誤アセンブリになる。
        /// </summary>
        [Fact]
        public void Merge_AmbiguousDestination_RefusesToJoin()
        {
            var 共通の左 = RandomSequence(5_000, seed: 631);
            var 右候補1 = RandomSequence(5_000, seed: 632);
            var 右候補2 = RandomSequence(5_000, seed: 633);
            var 中間 = RandomSequence(200, seed: 634);

            var 骨格 = this.WriteAssembly("backbone_amb.fasta", 63, 共通の左, 右候補1, 右候補2);
            // 同じ左から2つの異なる右へ繋がる証拠が両方ある。
            var 他 = this.WriteAssembly("other_amb.fasta", 31,
                共通の左 + 中間 + 右候補1,
                共通の左 + 中間 + 右候補2);

            var 出力 = Path.Combine(this._tempDir, "merged_amb.fasta");
            var 繋いだか = AssemblyMerger.V_統合(骨格, [骨格, 他], アンカーk長, 出力, p_必要な独立支持数: 1);

            Assert.False(繋いだか);
        }

        /// <summary>
        /// 3本を2箇所で繋ぐ連鎖。1回の統合で最後まで繋がること。
        /// </summary>
        [Fact]
        public void Merge_ChainOfThreePieces_JoinsAllOfThem()
        {
            var a = RandomSequence(4_000, seed: 641);
            var g1 = RandomSequence(200, seed: 642);
            var b = RandomSequence(4_000, seed: 643);
            var g2 = RandomSequence(200, seed: 644);
            var c = RandomSequence(4_000, seed: 645);
            var truth = a + g1 + b + g2 + c;

            var 骨格 = this.WriteAssembly("backbone_chain.fasta", 63, a, b, c);
            var 他 = this.WriteAssembly("other_chain.fasta", 31, truth);

            var 出力 = Path.Combine(this._tempDir, "merged_chain.fasta");
            var 繋いだか = AssemblyMerger.V_統合(骨格, [骨格, 他], アンカーk長, 出力, p_必要な独立支持数: 1);

            Assert.True(繋いだか);
            var 結果 = ReadSequences(出力);
            _ = Assert.Single(結果);
            Assert.Equal(truth, 結果[0]);
        }

        /// <summary>
        /// 繋がらなかった骨格配列も、統合結果から失われないこと。
        /// </summary>
        [Fact]
        public void Merge_UnjoinedBackbonePieces_AreStillEmitted()
        {
            var 左 = RandomSequence(4_000, seed: 651);
            var 中間 = RandomSequence(200, seed: 652);
            var 右 = RandomSequence(4_000, seed: 653);
            var 孤立 = RandomSequence(3_000, seed: 654);

            var 骨格 = this.WriteAssembly("backbone_iso.fasta", 63, 左, 右, 孤立);
            var 他 = this.WriteAssembly("other_iso.fasta", 31, 左 + 中間 + 右);

            var 出力 = Path.Combine(this._tempDir, "merged_iso.fasta");
            var 繋いだか = AssemblyMerger.V_統合(骨格, [骨格, 他], アンカーk長, 出力, p_必要な独立支持数: 1);

            Assert.True(繋いだか);
            var 結果 = ReadSequences(出力);
            Assert.Equal(2, 結果.Count);
            Assert.Contains(結果, x => x == 左 + 中間 + 右);
            Assert.Contains(結果, x => x == 孤立 || x == Util.V_逆相補(孤立));
        }

        /// <summary>
        /// 既定では、1つの k だけが主張する隣接は採らないこと。
        ///
        /// 骨格が途切れているのは繋ぐ根拠が足りないと判断した結果であることが多く、
        /// それを1本の配列で覆すと、その配列自身が誤アセンブリだった場合に
        /// そのまま持ち込む。実データでは、証拠に使ったアセンブリ由来の
        /// 誤アセンブリが骨格の21箇所から60箇所へ増えた。
        /// </summary>
        [Fact]
        public void Merge_OnlyOneKSupportsTheJoin_IsNotAcceptedByDefault()
        {
            var 左 = RandomSequence(5_000, seed: 671);
            var 中間 = RandomSequence(300, seed: 672);
            var 右 = RandomSequence(5_000, seed: 673);

            var 骨格 = this.WriteAssembly("backbone_sup.fasta", 63, 左, 右);
            var 他 = this.WriteAssembly("other_sup.fasta", 31, 左 + 中間 + 右);

            var 出力 = Path.Combine(this._tempDir, "merged_sup.fasta");

            Assert.False(AssemblyMerger.V_統合(骨格, [骨格, 他], アンカーk長, 出力));
        }

        /// <summary>
        /// 2つの k が同じ隣接を主張していれば採ること。
        /// </summary>
        [Fact]
        public void Merge_TwoIndependentKsAgree_IsAccepted()
        {
            var 左 = RandomSequence(5_000, seed: 681);
            var 中間 = RandomSequence(300, seed: 682);
            var 右 = RandomSequence(5_000, seed: 683);
            var truth = 左 + 中間 + 右;

            var 骨格 = this.WriteAssembly("backbone_two.fasta", 63, 左, 右);
            var 他1 = this.WriteAssembly("other_two_a.fasta", 31, truth);
            var 他2 = this.WriteAssembly("other_two_b.fasta", 41, truth);

            var 出力 = Path.Combine(this._tempDir, "merged_two.fasta");

            Assert.True(AssemblyMerger.V_統合(骨格, [骨格, 他1, 他2], アンカーk長, 出力));
            Assert.Equal(truth, ReadSequences(出力)[0]);
        }

        /// <summary>
        /// 統合の総延長が、骨格の総延長を下回らないこと。
        /// 配列を落とすなら統合しないほうがましなので、これは不変条件。
        /// </summary>
        [Fact]
        public void Merge_NeverLosesSequence()
        {
            var a = RandomSequence(4_000, seed: 661);
            var g = RandomSequence(150, seed: 662);
            var b = RandomSequence(4_000, seed: 663);
            var 孤立 = RandomSequence(2_000, seed: 664);

            var 骨格 = this.WriteAssembly("backbone_len.fasta", 63, a, b, 孤立);
            var 他 = this.WriteAssembly("other_len.fasta", 31, a + g + b);

            var 出力 = Path.Combine(this._tempDir, "merged_len.fasta");
            _ = AssemblyMerger.V_統合(骨格, [骨格, 他], アンカーk長, 出力, p_必要な独立支持数: 1);

            var 骨格の総延長 = a.Length + b.Length + 孤立.Length;
            Assert.True(ReadSequences(出力).Sum(x => x.Length) >= 骨格の総延長);
        }
    }
}
