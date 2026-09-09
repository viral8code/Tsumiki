using Tsumiki.Common;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 短い反復解決の拒否権(-rv)が使う r-mer 検証器そのものの検証<br/>
    /// head→repeat→tail の接合点を実際に跨いだリードが無ければ支持は
    /// 得られず、跨ぐリードがあれば支持が得られること、接合点を跨がない
    /// (=各配列の内部だけに収まる)リードだけでは支持にならないことを確認する<br/>
    /// フィクスチャの head/repeat/tail は(このアセンブリの)k=8 で
    /// 隣接する unitig 同士なので k-1=7 塩基を共有している<br/>
    /// r をこの
    /// 重なりより確実に長く取らないと、跨ぐ窓も共有区間の内側に収まって
    /// しまい判定にならない(RepeatRMerVerifier のクラスコメント参照)ため、
    /// ここでは r=18(=k+10、AssemblyPipeline の既定の決め方と同じ)を使う
    /// </summary>
    public class RepeatRMerVerifierTests : IDisposable
    {
        private const int R = 18;

        /// <summary>
        /// アセンブリ側の k(=head/repeat/tail が共有する重なりの長さ+1)
        /// </summary>
        private const int AssemblyK = 8;

        private readonly string _tempDir;

        public RepeatRMerVerifierTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_rmer_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
            // Get_接合点の支持数 は head/tail から共有重なり(k-1 塩基)を除くのに
            // 現在の実行時引数の k 長を参照するため、フィクスチャの重なり長
            // (AssemblyK-1=7)に合わせておく
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = AssemblyK, A_スレッド数 = 1 };
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private const string Head = "ACAGTTCGCGAGCCCTCCGTC"; // 21bp
        private const string Repeat = "CTCCGTCAGCTTGTTTGGAGCAGA"; // 24bp
        private const string Tail = "GAGCAGAGTCGTTCTGCGAGG"; // 21bp
        private const string OtherTail = "TTTTTTTTTTTTTTTTTTTTT"; // 21bp、Repeatと接合しない無関係な配列

        private string WriteFastq(string name, IEnumerable<string> reads)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new StreamWriter(path);
            var idx = 0;
            foreach (var seq in reads)
            {
                writer.WriteLine($"@r{idx}");
                writer.WriteLine(seq);
                writer.WriteLine("+");
                writer.WriteLine(new string('I', seq.Length));
                idx++;
            }
            return path;
        }

        private static IEnumerable<string> SlidingReads(string sequence, int readLength)
        {
            for (var i = 0; i + readLength <= sequence.Length; i++)
            {
                yield return sequence.Substring(i, readLength);
            }
        }

        [Fact]
        public void Get_接合点の支持数_IsPositive_WhenReadsActuallyCrossBothJunctions()
        {
            var walk = Head + Repeat[(AssemblyK - 1)..] + Tail[(AssemblyK - 1)..];
            var path = this.WriteFastq("cross.fq", SlidingReads(walk, 25));

            var verifier = RepeatRMerVerifier.V_構築([path, string.Empty], R);

            var support = verifier.Get_接合点の支持数(Head, Repeat, Tail);

            Assert.True(support >= Consts.r_mer接合点支持の閾値の既定値,
                $"expected support ({support}) to reach the default threshold when reads truly cross the junctions");
            Assert.True(verifier.Get_接合点に支持があるか(Head, Repeat, Tail, Consts.r_mer接合点支持の閾値の既定値));
        }

        /// <summary>
        /// リードが Head・Repeat・Tail それぞれの内部だけに収まり、
        /// どちらの接合点も跨がない場合、支持は既定の閾値に届かないはず
        /// (各配列は単独でも実在する配列なので、接合点を跨がない限り
        /// 「その組み合わせが正しい」証拠にはならない)
        /// </summary>
        [Fact]
        public void Get_接合点の支持数_StaysBelowThreshold_WhenReadsNeverCrossEitherJunction()
        {
            // 各配列を丸ごと読んだリードを与える
            // r より短いリードだと r-mer が
            // 1 つも作られず、何を数えても 0 になって検定にならない
            List<string> reads = [Head, Repeat, Tail];
            var path = this.WriteFastq("internal_only.fq", reads);

            var verifier = RepeatRMerVerifier.V_構築([path, string.Empty], R);

            var support = verifier.Get_接合点の支持数(Head, Repeat, Tail);

            Assert.True(support < Consts.r_mer接合点支持の閾値の既定値,
                $"expected support ({support}) to stay below the threshold when no read actually crosses a junction");
            Assert.False(verifier.Get_接合点に支持があるか(Head, Repeat, Tail, Consts.r_mer接合点支持の閾値の既定値));
        }

        /// <summary>
        /// Head-Repeat 接合点は本物のリードに跨がれているが、
        /// Repeat-Tail 側は無関係な配列(OtherTail)であり跨ぐリードが無い場合でも
        /// 支持は得られる(Head-Repeat側の支持だけでカウントされるため)<br/>
        /// この支持数は、両方の接合点が本物のリードに跨がれている場合の
        /// 支持数を超えないはず
        /// </summary>
        [Fact]
        public void Get_接合点の支持数_CountsOnlyTheJunctionThatIsActuallyCrossed()
        {
            var walk = Head + Repeat[(AssemblyK - 1)..];
            var path = this.WriteFastq("one_side.fq", SlidingReads(walk, 25));

            var verifier = RepeatRMerVerifier.V_構築([path, string.Empty], R);

            var withOnlyHeadSideCrossable = verifier.Get_接合点の支持数(Head, Repeat, OtherTail);

            Assert.True(withOnlyHeadSideCrossable > 0);

            var bothWalk = Head + Repeat[(AssemblyK - 1)..] + Tail[(AssemblyK - 1)..];
            var bothPath = this.WriteFastq("both_sides.fq", SlidingReads(bothWalk, 25));
            var bothVerifier = RepeatRMerVerifier.V_構築([bothPath, string.Empty], R);
            var withBothCrossable = bothVerifier.Get_接合点の支持数(Head, Repeat, Tail);

            Assert.True(withBothCrossable > withOnlyHeadSideCrossable);
        }

        [Fact]
        public void Get_接合点の支持数_MatchesReverseComplementReadsToo()
        {
            var walk = Head + Repeat[(AssemblyK - 1)..] + Tail[(AssemblyK - 1)..];
            var rc = Util.V_逆相補(walk);
            var path = this.WriteFastq("rc.fq", SlidingReads(rc, 25));

            var verifier = RepeatRMerVerifier.V_構築([path, string.Empty], R);

            Assert.True(verifier.Get_接合点に支持があるか(Head, Repeat, Tail, Consts.r_mer接合点支持の閾値の既定値));
        }

        [Fact]
        public void V_構築_IgnoresMissingOrEmptyPaths()
        {
            var verifier = RepeatRMerVerifier.V_構築([string.Empty, Path.Combine(this._tempDir, "does_not_exist.fq")], R);

            Assert.False(verifier.Get_接合点に支持があるか(Head, Repeat, Tail, 1));
        }

        [Fact]
        public void V_構築_RejectsNonPositiveRLength()
        {
            _ = Assert.Throws<ArgumentException>(() => RepeatRMerVerifier.V_構築([], 0));
            _ = Assert.Throws<ArgumentException>(() => RepeatRMerVerifier.V_構築([], -1));
        }

        /// <summary>
        /// 2bit パックが ulong に収まらない長さ (33 以上) でも、ふるいへ
        /// 切り替えて同じ判定ができること<br/>
        /// 跨いだリードがあれば支持が出て、
        /// 無ければ出ない
        /// </summary>
        [Fact]
        public void Get_接合点の支持数_WorksBeyondThePackableRLength()
        {
            const int l_長いR = 40;
            var l_head = RandomSequence(120, seed: 20260922);
            var l_repeat = l_head[^(AssemblyK - 1)..] + RandomSequence(120, seed: 20260923);
            var l_tail = l_repeat[^(AssemblyK - 1)..] + RandomSequence(120, seed: 20260924);

            var l_跨ぐ = l_head + l_repeat[(AssemblyK - 1)..] + l_tail[(AssemblyK - 1)..];
            var l_跨ぐパス = this.WriteFastq("long_cross.fq", SlidingReads(l_跨ぐ, 100));
            var l_跨がないパス = this.WriteFastq(
                "long_apart.fq", SlidingReads(l_head, 100).Concat(SlidingReads(l_tail, 100)));

            var l_跨ぐ検証器 = RepeatRMerVerifier.V_構築([l_跨ぐパス, string.Empty], l_長いR);
            var l_跨がない検証器 = RepeatRMerVerifier.V_構築([l_跨がないパス, string.Empty], l_長いR);

            Assert.True(l_跨ぐ検証器.Get_接合点に支持があるか(
                l_head, l_repeat, l_tail, Consts.r_mer接合点支持の閾値の既定値));
            Assert.Equal(0, l_跨がない検証器.Get_接合点の支持数(l_head, l_repeat, l_tail));
        }

        private static string RandomSequence(int p_長さ, int seed)
        {
            var l_乱数 = new Random(seed);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }
    }
}
