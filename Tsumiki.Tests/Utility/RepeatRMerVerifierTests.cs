using Tsumiki.Common;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 短い反復解決の拒否権 (-rv) が使う r-mer 検証器そのものの検証
    /// </summary>
    /// <remarks>
    /// head→repeat→tail の接合点を実際に跨いだリードが無ければ支持は
    /// 得られず、跨ぐリードがあれば支持が得られること、接合点を跨がない
    /// (=各配列の内部だけに収まる) リードだけでは支持にならないことを確認する<br/>
    /// フィクスチャの head/repeat/tail は (このアセンブリの) k=8 で
    /// 隣接する unitig 同士なので k-1=7 塩基を共有している<br/>
    /// r をこの
    /// 重なりより確実に長く取らないと、跨ぐ窓も共有区間の内側に収まって
    /// しまい判定にならない (RepeatRMerVerifier のクラスコメント参照) ため、
    /// ここでは r=18(=k+10、AssemblyPipeline の既定の決め方と同じ) を使う
    /// </remarks>
    public class RepeatRMerVerifierTests : IDisposable
    {
        /// <summary>
        /// この検証で使う r-mer 長
        /// </summary>
        private const int r長 = 18;

        /// <summary>
        /// アセンブリ側の k(=head/repeat/tail が共有する重なりの長さ+1)
        /// </summary>
        private const int アセンブリk長 = 8;

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public RepeatRMerVerifierTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_rmer_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
            // Get_接合点の支持数 は head/tail から共有重なり (k-1 塩基) を除くのに
            // 現在の実行時引数の k 長を参照するため、フィクスチャの重なり長
            // (アセンブリk長-1=7) に合わせておく
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = アセンブリk長, A_スレッド数 = 1 };
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        /// <summary>
        /// 反復の手前にある配列
        /// </summary>
        private const string 反復前配列 = "ACAGTTCGCGAGCCCTCCGTC";

        /// <summary>
        /// 手前と先に挟まれた反復配列
        /// </summary>
        private const string 反復配列 = "CTCCGTCAGCTTGTTTGGAGCAGA";

        /// <summary>
        /// 反復の先にある配列
        /// </summary>
        private const string 反復後配列 = "GAGCAGAGTCGTTCTGCGAGG";

        /// <summary>
        /// 反復配列と接合しない無関係な配列
        /// </summary>
        private const string 無関係な後続配列 = "TTTTTTTTTTTTTTTTTTTTT";

        /// <summary>
        /// リードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_name">ファイル名</param>
        /// <param name="p_reads">書き出すリード</param>
        /// <returns>書き出したパス</returns>
        private string V_書き込み_fastq(string p_name, IEnumerable<string> p_reads)
        {
            var l_path = Path.Combine(this._tempDir, p_name);
            using var l_writer = new StreamWriter(l_path);
            var l_idx = 0;
            foreach (var l_seq in p_reads)
            {
                l_writer.WriteLine($"@r{l_idx}");
                l_writer.WriteLine(l_seq);
                l_writer.WriteLine("+");
                l_writer.WriteLine(new string('I', l_seq.Length));
                l_idx++;
            }
            return l_path;
        }

        /// <summary>
        /// 配列を 1 塩基ずつずらして覆うリードを作る
        /// </summary>
        /// <param name="p_sequence">覆う配列</param>
        /// <param name="p_readLength">リード長</param>
        /// <returns>リード</returns>
        private static IEnumerable<string> V_生成_スライドリード(string p_sequence, int p_readLength)
        {
            for (var i = 0; i + p_readLength <= p_sequence.Length; i++)
            {
                yield return p_sequence.Substring(i, p_readLength);
            }
        }

        /// <summary>
        /// 両方の接合点を実際に跨ぐリードがあれば、支持数は正になる
        /// </summary>
        [Fact]
        public void 両方の接合点を実際に跨ぐリードがあれば支持数は正になる()
        {
            var l_walk = 反復前配列 + 反復配列[(アセンブリk長 - 1)..] + 反復後配列[(アセンブリk長 - 1)..];
            var l_path = this.V_書き込み_fastq("cross.fq", V_生成_スライドリード(l_walk, 25));

            var l_verifier = RepeatRMerVerifier.V_構築([l_path, string.Empty], r長);

            var l_support = l_verifier.Get_接合点の支持数(反復前配列, 反復配列, 反復後配列);

            Assert.True(l_support >= Consts.r_mer接合点支持の閾値の既定値,
                $"expected support ({l_support}) to reach the default threshold when reads truly cross the junctions");
            Assert.True(l_verifier.Get_接合点に支持があるか(反復前配列, 反復配列, 反復後配列, Consts.r_mer接合点支持の閾値の既定値));
        }

        /// <summary>
        /// どちらの接合点も跨がないリードだけでは、支持が既定の閾値に届かないこと
        /// </summary>
        /// <remarks>
        /// 各配列は単独でも実在する配列なので、接合点を跨がない限り
        /// その組み合わせが正しい証拠にはならない
        /// </remarks>
        [Fact]
        public void どちらの接合点も跨がないリードだけでは支持が閾値に届かない()
        {
            // 各配列を丸ごと読んだリードを与える
            // r より短いリードだと r-mer が
            // 1 つも作られず、何を数えても 0 になって検定にならない
            List<string> l_reads = [反復前配列, 反復配列, 反復後配列];
            var l_path = this.V_書き込み_fastq("internal_only.fq", l_reads);

            var l_verifier = RepeatRMerVerifier.V_構築([l_path, string.Empty], r長);

            var l_support = l_verifier.Get_接合点の支持数(反復前配列, 反復配列, 反復後配列);

            Assert.True(l_support < Consts.r_mer接合点支持の閾値の既定値,
                $"expected support ({l_support}) to stay below the threshold when no read actually crosses a junction");
            Assert.False(l_verifier.Get_接合点に支持があるか(反復前配列, 反復配列, 反復後配列, Consts.r_mer接合点支持の閾値の既定値));
        }

        /// <summary>
        /// Head-Repeat 接合点は本物のリードに跨がれているが、
        /// Repeat-Tail 側は無関係な配列 (OtherTail) であり跨ぐリードが無い場合でも
        /// 支持は得られる (Head-Repeat 側の支持だけでカウントされるため)
        /// </summary>
        /// <remarks>
        /// この支持数は、両方の接合点が本物のリードに跨がれている場合の
        /// 支持数を超えないはず
        /// </remarks>
        [Fact]
        public void 実際に跨いだ接合点だけが支持数に数えられる()
        {
            var l_walk = 反復前配列 + 反復配列[(アセンブリk長 - 1)..];
            var l_path = this.V_書き込み_fastq("one_side.fq", V_生成_スライドリード(l_walk, 25));

            var l_verifier = RepeatRMerVerifier.V_構築([l_path, string.Empty], r長);

            var l_withOnlyHeadSideCrossable = l_verifier.Get_接合点の支持数(反復前配列, 反復配列, 無関係な後続配列);

            Assert.True(l_withOnlyHeadSideCrossable > 0);

            var l_bothWalk = 反復前配列 + 反復配列[(アセンブリk長 - 1)..] + 反復後配列[(アセンブリk長 - 1)..];
            var l_bothPath = this.V_書き込み_fastq("both_sides.fq", V_生成_スライドリード(l_bothWalk, 25));
            var l_bothVerifier = RepeatRMerVerifier.V_構築([l_bothPath, string.Empty], r長);
            var l_withBothCrossable = l_bothVerifier.Get_接合点の支持数(反復前配列, 反復配列, 反復後配列);

            Assert.True(l_withBothCrossable > l_withOnlyHeadSideCrossable);
        }

        /// <summary>
        /// 逆相補のリードでも、接合点の支持が得られる
        /// </summary>
        [Fact]
        public void 逆相補のリードでも接合点の支持が得られる()
        {
            var l_walk = 反復前配列 + 反復配列[(アセンブリk長 - 1)..] + 反復後配列[(アセンブリk長 - 1)..];
            var l_rc = Util.V_逆相補(l_walk);
            var l_path = this.V_書き込み_fastq("rc.fq", V_生成_スライドリード(l_rc, 25));

            var l_verifier = RepeatRMerVerifier.V_構築([l_path, string.Empty], r長);

            Assert.True(l_verifier.Get_接合点に支持があるか(反復前配列, 反復配列, 反復後配列, Consts.r_mer接合点支持の閾値の既定値));
        }

        /// <summary>
        /// 存在しないパスや空のパスは無視して構築できる
        /// </summary>
        [Fact]
        public void 存在しないパスや空のパスは無視して構築できる()
        {
            var l_verifier = RepeatRMerVerifier.V_構築([string.Empty, Path.Combine(this._tempDir, "does_not_exist.fq")], r長);

            Assert.False(l_verifier.Get_接合点に支持があるか(反復前配列, 反復配列, 反復後配列, 1));
        }

        /// <summary>
        /// r 長が 0 以下なら構築を拒否する
        /// </summary>
        [Fact]
        public void r長が0以下なら構築を拒否する()
        {
            _ = Assert.Throws<ArgumentException>(() => RepeatRMerVerifier.V_構築([], 0));
            _ = Assert.Throws<ArgumentException>(() => RepeatRMerVerifier.V_構築([], -1));
        }

        /// <summary>
        /// 2 bit パックが ulong に収まらない長さ (33 以上) でも、ふるいへ
        /// 切り替えて同じ判定ができること
        /// </summary>
        /// <remarks>
        /// 跨いだリードがあれば支持が出て、
        /// 無ければ出ない
        /// </remarks>
        [Fact]
        public void rmer長が2bitパックに収まらなくても判定できる()
        {
            const int l_長いR = 40;
            var l_head = V_生成_乱数配列(120, p_seed: 20260922);
            var l_repeat = l_head[^(アセンブリk長 - 1)..] + V_生成_乱数配列(120, p_seed: 20260923);
            var l_tail = l_repeat[^(アセンブリk長 - 1)..] + V_生成_乱数配列(120, p_seed: 20260924);

            var l_跨ぐ = l_head + l_repeat[(アセンブリk長 - 1)..] + l_tail[(アセンブリk長 - 1)..];
            var l_跨ぐパス = this.V_書き込み_fastq("long_cross.fq", V_生成_スライドリード(l_跨ぐ, 100));
            var l_跨がないパス = this.V_書き込み_fastq(
                "long_apart.fq", V_生成_スライドリード(l_head, 100).Concat(V_生成_スライドリード(l_tail, 100)));

            var l_跨ぐ検証器 = RepeatRMerVerifier.V_構築([l_跨ぐパス, string.Empty], l_長いR);
            var l_跨がない検証器 = RepeatRMerVerifier.V_構築([l_跨がないパス, string.Empty], l_長いR);

            Assert.True(l_跨ぐ検証器.Get_接合点に支持があるか(
                l_head, l_repeat, l_tail, Consts.r_mer接合点支持の閾値の既定値));
            Assert.Equal(0, l_跨がない検証器.Get_接合点の支持数(l_head, l_repeat, l_tail));
        }

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_seed">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_長さ, int p_seed)
        {
            var l_乱数 = new Random(p_seed);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }
    }
}
