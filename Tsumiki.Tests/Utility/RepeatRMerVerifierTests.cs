using Tsumiki.Commons;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Tests.Utility
{
    /// <summary>
    /// 短い反復解決の拒否権 (-rv) が使う r-mer 検証器そのものの検証
    /// </summary>
    public class RepeatRMerVerifierTests : IDisposable
    {
        #region 定数

        /// <summary>
        /// この検証で使う r-mer 長
        /// </summary>
        private const int r長 = 18;

        /// <summary>
        /// アセンブリ側の k (=head/repeat/tail が共有する重なりの長さ+1)
        /// </summary>
        private const int アセンブリk長 = 8;

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

        #endregion

        #region 内部変数

        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _作業ディレクトリ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 検証用の状態を初期化する
        /// </summary>
        public RepeatRMerVerifierTests()
        {
            this._作業ディレクトリ = Path.Combine(Path.GetTempPath(), "tsumiki_rmer_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._作業ディレクトリ);
            // Get_接合点の支持数 は head/tail から共有重なり (k-1 塩基) を除くのに
            // 現在の実行時引数の k 長を参照するため、フィクスチャの重なり長
            // (アセンブリ k 長-1=7) に合わせておく
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = アセンブリk長, A_スレッド数 = 1 };
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// パック幅の境界でも全塩基を保持し、逆相補を同一視する
        /// </summary>
        /// <param name="p_長さ">検査する窓の長さ</param>
        [Theory]
        [InlineData(32)]
        [InlineData(33)]
        [InlineData(64)]
        [InlineData(65)]
        [InlineData(128)]
        [InlineData(129)]
        public void V_厳密キー_全塩基を区別(int p_長さ)
        {
            var l_配列 = "A" + new string('C', p_長さ - 2) + "G";
            var l_キー = RepeatRMerVerifier.Get_正準値(l_配列);
            Assert.Equal(l_キー, RepeatRMerVerifier.Get_正準値(Util.V_逆相補(l_配列)));
            for (var i = 0; i < p_長さ; i++)
            {
                var l_変更 = l_配列.ToCharArray();
                l_変更[i] = 'T';
                Assert.NotEqual(l_キー, RepeatRMerVerifier.Get_正準値(l_変更));
            }
        }

        /// <summary>
        /// 一時ディレクトリを片付ける
        /// </summary>
        public void Dispose()
        {
            if (Directory.Exists(this._作業ディレクトリ))
            {
                Directory.Delete(this._作業ディレクトリ, recursive: true);
            }
        }

        /// <summary>
        /// 両方の接合点を実際に跨ぐリードがあれば、支持数は正になる
        /// </summary>
        [Fact]
        public void V_両方の接合点を実際に跨ぐリードがあれば支持数は正になる()
        {
            var l_経路 = 反復前配列 + 反復配列[(アセンブリk長 - 1)..] + 反復後配列[(アセンブリk長 - 1)..];
            var l_パス = this.V_書き込み_fastq("cross.fq", V_生成_スライドリード(l_経路, 25));

            var l_検証器 = RepeatRMerVerifier.V_構築([l_パス, string.Empty], r長);

            var l_支持 = l_検証器.Get_接合点の支持数(反復前配列, 反復配列, 反復後配列);

            Assert.True(l_支持 >= Consts.r_mer接合点支持の閾値の既定値, $"expected support ({l_支持}) to reach the default threshold when reads truly cross the junctions");
            Assert.True(l_検証器.Has接合点支持(反復前配列, 反復配列, 反復後配列, Consts.r_mer接合点支持の閾値の既定値));
        }

        /// <summary>
        /// どちらの接合点も跨がないリードだけでは、支持が既定の閾値に届かないこと
        /// </summary>
        /// <remarks>
        /// 各配列は単独でも実在する配列なので、接合点を跨がない限りその組み合わせが正しい証拠にはならない
        /// </remarks>
        [Fact]
        public void V_どちらの接合点も跨がないリードだけでは支持が閾値に届かない()
        {
            // 各配列を丸ごと読んだリードを与える
            // r より短いリードだと r-mer が
            // 1 つも作られず、何を数えても 0 になって検定にならない
            List<string> l_リード群 = [反復前配列, 反復配列, 反復後配列];
            var l_パス = this.V_書き込み_fastq("internal_only.fq", l_リード群);

            var l_検証器 = RepeatRMerVerifier.V_構築([l_パス, string.Empty], r長);

            var l_支持 = l_検証器.Get_接合点の支持数(反復前配列, 反復配列, 反復後配列);

            Assert.True(l_支持 < Consts.r_mer接合点支持の閾値の既定値, $"expected support ({l_支持}) to stay below the threshold when no read actually crosses a junction");
            Assert.False(l_検証器.Has接合点支持(反復前配列, 反復配列, 反復後配列, Consts.r_mer接合点支持の閾値の既定値));
        }

        /// <summary>
        /// Head-Repeat 接合点は本物のリードに跨がれているが、Repeat-Tail 側は無関係な配列 (OtherTail) であり跨ぐリードが無い場合でも支持は得られる (Head-Repeat 側の支持だけでカウントされるため)
        /// </summary>
        /// <remarks>
        /// この支持数は、両方の接合点が本物のリードに跨がれている場合の支持数を超えないはず
        /// </remarks>
        [Fact]
        public void V_実際に跨いだ接合点だけが支持数に数えられる()
        {
            var l_経路 = 反復前配列 + 反復配列[(アセンブリk長 - 1)..];
            var l_パス = this.V_書き込み_fastq("one_side.fq", V_生成_スライドリード(l_経路, 25));

            var l_検証器 = RepeatRMerVerifier.V_構築([l_パス, string.Empty], r長);

            var l_withOnlyHeadSideCrossable = l_検証器.Get_接合点の支持数(反復前配列, 反復配列, 無関係な後続配列);

            Assert.True(l_withOnlyHeadSideCrossable > 0);
            Assert.False(l_検証器.Has接合点支持(反復前配列, 反復配列, 無関係な後続配列, 1));

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
        public void V_逆相補のリードでも接合点の支持が得られる()
        {
            var l_経路 = 反復前配列 + 反復配列[(アセンブリk長 - 1)..] + 反復後配列[(アセンブリk長 - 1)..];
            var l_rc = Util.V_逆相補(l_経路);
            var l_パス = this.V_書き込み_fastq("rc.fq", V_生成_スライドリード(l_rc, 25));

            var l_検証器 = RepeatRMerVerifier.V_構築([l_パス, string.Empty], r長);

            Assert.True(l_検証器.Has接合点支持(反復前配列, 反復配列, 反復後配列, Consts.r_mer接合点支持の閾値の既定値));
        }

        /// <summary>
        /// 存在しないパスや空のパスは無視して構築できる
        /// </summary>
        [Fact]
        public void V_存在しないパスや空のパスは無視して構築できる()
        {
            var l_検証器 = RepeatRMerVerifier.V_構築([string.Empty, Path.Combine(this._作業ディレクトリ, "does_not_exist.fq")], r長);

            Assert.False(l_検証器.Has接合点支持(反復前配列, 反復配列, 反復後配列, 1));
        }

        /// <summary>
        /// r 長が 0 以下なら構築を拒否する
        /// </summary>
        [Fact]
        public void V_r長が0以下なら構築を拒否する()
        {
            _ = Assert.Throws<ArgumentException>(() => RepeatRMerVerifier.V_構築([], 0));
            _ = Assert.Throws<ArgumentException>(() => RepeatRMerVerifier.V_構築([], -1));
        }

        /// <summary>
        /// 2 bit パックが ulong に収まらない長さ (33 以上) でも、ふるいへ切り替えて同じ判定ができること
        /// </summary>
        /// <remarks>
        /// 跨いだリードがあれば支持が出て、無ければ出ない
        /// </remarks>
        [Fact]
        public void V_rmer長が2bitパックに収まらなくても判定できる()
        {
            const int l_長いR = 40;
            var l_先頭 = V_生成_乱数配列(120, p_乱数種: 20_260_922);
            var l_反復配列 = l_先頭[^(アセンブリk長 - 1)..] + V_生成_乱数配列(120, p_乱数種: 20_260_923);
            var l_末尾 = l_反復配列[^(アセンブリk長 - 1)..] + V_生成_乱数配列(120, p_乱数種: 20_260_924);

            var l_跨ぐ = l_先頭 + l_反復配列[(アセンブリk長 - 1)..] + l_末尾[(アセンブリk長 - 1)..];
            var l_跨ぐパス = this.V_書き込み_fastq("long_cross.fq", V_生成_スライドリード(l_跨ぐ, 100));
            var l_跨がないパス = this.V_書き込み_fastq("long_apart.fq", V_生成_スライドリード(l_先頭, 100).Concat(V_生成_スライドリード(l_末尾, 100)));

            var l_跨ぐ検証器 = RepeatRMerVerifier.V_構築([l_跨ぐパス, string.Empty], l_長いR);
            var l_跨がない検証器 = RepeatRMerVerifier.V_構築([l_跨がないパス, string.Empty], l_長いR);

            Assert.True(l_跨ぐ検証器.Has接合点支持(l_先頭, l_反復配列, l_末尾, Consts.r_mer接合点支持の閾値の既定値));
            Assert.Equal(0, l_跨がない検証器.Get_接合点の支持数(l_先頭, l_反復配列, l_末尾));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// リードを FASTQ として書き出す
        /// </summary>
        /// <param name="p_名前">ファイル名</param>
        /// <param name="p_リード群">書き出すリード</param>
        /// <returns>書き出したパス</returns>
        private string V_書き込み_fastq(string p_名前, IEnumerable<string> p_リード群)
        {
            var l_パス = Path.Combine(this._作業ディレクトリ, p_名前);
            using var l_書き込み = new StreamWriter(l_パス);
            var l_位置 = 0;
            foreach (var l_配列 in p_リード群)
            {
                l_書き込み.WriteLine($"@r{l_位置}");
                l_書き込み.WriteLine(l_配列);
                l_書き込み.WriteLine("+");
                l_書き込み.WriteLine(new string('I', l_配列.Length));
                l_位置++;
            }
            return l_パス;
        }

        /// <summary>
        /// 配列を 1 塩基ずつずらして覆うリードを作る
        /// </summary>
        /// <param name="p_配列">覆う配列</param>
        /// <param name="p_リード長">リード長</param>
        /// <returns>リード</returns>
        private static IEnumerable<string> V_生成_スライドリード(string p_配列, int p_リード長)
        {
            for (var i = 0; i + p_リード長 <= p_配列.Length; i++)
            {
                yield return p_配列.Substring(i, p_リード長);
            }
        }

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_乱数種">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_生成_乱数配列(int p_長さ, int p_乱数種)
        {
            var l_乱数 = new Random(p_乱数種);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_乱数.Next(4)]));
        }

        #endregion

    }
}
