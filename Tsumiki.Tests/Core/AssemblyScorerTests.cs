using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// リファレンス無しでアセンブリの良さを測る評価器の検証。
    ///
    /// multi-k で複数のアセンブリから1つを選ぶには、リファレンスを使わずに
    /// 良し悪しを決められなければならない。連続性(N50)だけで選ぶと
    /// 誤って繋いだものほど高く出るため、完全性と正確性を併せて見る必要がある。
    /// ここではその「誤って繋いだものが落ちる」ことを主に固定する。
    /// </summary>
    public class AssemblyScorerTests : IDisposable
    {
        private readonly string _tempDir;

        public AssemblyScorerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_scorer_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private const int K = 21;
        private const int 深さ = 20;

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        /// <summary>
        /// 与えた配列群から k-mer インデックスを作る。深さは一律なので
        /// 単一コピー基準値は 深さ そのものになる。
        /// </summary>
        private TrustedKmerIndex BuildIndex(params string[] sequences)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };
            var index = new TrustedKmerIndex(this._tempDir);
            foreach (var seq in sequences)
            {
                var bytes = seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + K <= bytes.Length; i++)
                {
                    for (var rep = 0; rep < 深さ; rep++)
                    {
                        index.V_登録(bytes.AsSpan(i, K), p_ワーカー番号: 0);
                    }
                }
            }
            _ = index.V_カットオフ(p_カットオフ: 2);
            return index;
        }

        private string WriteFasta(string name, params string[] sequences)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new FastaWriter(path);
            var id = 1;
            foreach (var seq in sequences)
            {
                writer.V_書き込み($"NODE{id++}", seq);
            }
            return path;
        }

        [Fact]
        public void Score_PerfectAssembly_HasFullCompletenessAndAccuracy()
        {
            var truth = RandomSequence(20_000, seed: 501);
            using var index = this.BuildIndex(truth);

            var path = this.WriteFasta("perfect.fasta", truth);
            var 評価 = AssemblyScorer.Get_評価(path, index, K, 深さ, truth.Length);

            Assert.NotNull(評価);
            Assert.Equal(0, 評価.A_欠損延べ数);
            Assert.Equal(0, 評価.A_過剰延べ数);
            Assert.Equal(1.0, 評価.A_完全性, 6);
            Assert.Equal(1.0, 評価.A_正確性, 6);
        }

        /// <summary>
        /// 断片化しているが取りこぼしの無いアセンブリ。完全性は満点のまま、
        /// 連続性だけが落ちること。
        /// </summary>
        [Fact]
        public void Score_FragmentedButComplete_LosesContiguityOnly()
        {
            var truth = RandomSequence(20_000, seed: 502);
            using var index = this.BuildIndex(truth);

            // k-1 塩基重ねて切ると、境界の k-mer も失われない。
            var 断片 = new List<string>();
            for (var i = 0; i < truth.Length; i += 4_000)
            {
                var 終端 = Math.Min(truth.Length, i + 4_000 + K - 1);
                断片.Add(truth[i..終端]);
            }

            var 一本 = this.WriteFasta("whole.fasta", truth);
            var 断片化 = this.WriteFasta("fragmented.fasta", [.. 断片]);

            var 評価一本 = AssemblyScorer.Get_評価(一本, index, K, 深さ, truth.Length);
            var 評価断片 = AssemblyScorer.Get_評価(断片化, index, K, 深さ, truth.Length);

            Assert.NotNull(評価一本);
            Assert.NotNull(評価断片);
            Assert.Equal(0, 評価断片.A_欠損延べ数);
            Assert.True(評価断片.A_NG50 < 評価一本.A_NG50);
        }

        /// <summary>
        /// これが評価器の存在意義。反復配列を通り抜けて中間を飛ばした
        /// 誤アセンブリは、素の連続性では「改善」に見える(実際、過去に
        /// N50 が 99,974 から 199,945 へ伸びた誤アセンブリがあった)。
        /// 飛ばした領域の k-mer が欠損として現れるため、完全性で見抜ける。
        ///
        /// 重要なのは、連続性ではこのキメラのほうが上だという点である。
        /// だからこそ選択規則は「まず完全性で足切りし、そのあとで連続性を見る」
        /// という順序でなければならない(掛け算にすると連続性の利得が勝ってしまう)。
        /// </summary>
        [Fact]
        public void Score_ChimeraThatSkipsSequence_ScoresBelowTheFragmentedButHonestAssembly()
        {
            // A-R-B-R-C。R は2回現れる反復配列。
            var a = RandomSequence(8_000, seed: 511);
            var r = RandomSequence(300, seed: 512);
            var b = RandomSequence(8_000, seed: 513);
            var c = RandomSequence(8_000, seed: 514);
            var truth = a + r + b + r + c;

            using var index = this.BuildIndex(truth);

            // 正直な答え: R で切れているが、A も B も C も出ている。
            var 正直 = this.WriteFasta("honest.fasta", a + r, r + b + r, r + c);
            // 誤アセンブリ: R を1回通り抜けて B を丸ごと飛ばした A-R-C。
            var キメラ = this.WriteFasta("chimera.fasta", a + r + c);

            var 評価正直 = AssemblyScorer.Get_評価(正直, index, K, 深さ, truth.Length);
            var 評価キメラ = AssemblyScorer.Get_評価(キメラ, index, K, 深さ, truth.Length);

            Assert.NotNull(評価正直);
            Assert.NotNull(評価キメラ);

            // 飛ばした B のぶんだけキメラ側に欠損が出る。
            Assert.True(評価キメラ.A_欠損延べ数 > 評価正直.A_欠損延べ数,
                $"chimera missing={評価キメラ.A_欠損延べ数}, honest missing={評価正直.A_欠損延べ数}");
            Assert.True(評価キメラ.A_完全性 < 評価正直.A_完全性);
            // 完全性の差は足切りの許容差を大きく超えていること。
            Assert.True(評価正直.A_完全性 - 評価キメラ.A_完全性 > AssemblySelector.完全性の許容差,
                $"chimera completeness={評価キメラ.A_完全性:F3}, honest={評価正直.A_完全性:F3}");
            // 連続性だけを見るとキメラのほうが良く見えることを明示しておく。
            Assert.True(評価キメラ.A_NG50 > 評価正直.A_NG50,
                "this test is only meaningful while the chimera looks better on contiguity alone");
        }

        /// <summary>
        /// 同じ配列を2回出した水増しは、正確性が落ちて総合点が下がること。
        /// 連続性(NG50)はむしろ上がるため、この判定が無いと選んでしまう。
        /// </summary>
        [Fact]
        public void Score_DuplicatedSequence_IsPenalisedByAccuracy()
        {
            var truth = RandomSequence(20_000, seed: 521);
            using var index = this.BuildIndex(truth);

            var 正常 = this.WriteFasta("single.fasta", truth);
            var 水増し = this.WriteFasta("doubled.fasta", truth, truth);

            var 評価正常 = AssemblyScorer.Get_評価(正常, index, K, 深さ, truth.Length);
            var 評価水増し = AssemblyScorer.Get_評価(水増し, index, K, 深さ, truth.Length);

            Assert.NotNull(評価正常);
            Assert.NotNull(評価水増し);
            Assert.True(評価水増し.A_過剰延べ数 > 0);
            Assert.True(評価水増し.A_正確性 < 評価正常.A_正確性);
        }

        /// <summary>
        /// 2コピーの反復配列を2回出すのは正しい。水増しと誤判定しないこと。
        /// </summary>
        [Fact]
        public void Score_TwoCopyRepeatEmittedTwice_IsNotPenalised()
        {
            var 単一 = RandomSequence(10_000, seed: 531);
            var 反復 = RandomSequence(500, seed: 532);
            // 反復が2回現れるゲノム。
            var truth = 単一 + 反復 + RandomSequence(5_000, seed: 533) + 反復;

            using var index = this.BuildIndex(truth);

            var path = this.WriteFasta("repeat.fasta", truth);
            var 評価 = AssemblyScorer.Get_評価(path, index, K, 深さ, truth.Length);

            Assert.NotNull(評価);
            Assert.Equal(0, 評価.A_過剰延べ数);
            Assert.Equal(0, 評価.A_欠損延べ数);
        }

        /// <summary>
        /// NG50 は自分の総延長ではなく推定ゲノムサイズを分母にすること。
        /// 素の N50 だと「配列を落として短くなったアセンブリ」ほど有利になり、
        /// k を跨いだ比較に使えない。
        /// </summary>
        [Fact]
        public void Score_NG50_UsesTheGenomeSizeAsDenominator_NotTheAssemblyLength()
        {
            var truth = RandomSequence(20_000, seed: 541);
            using var index = this.BuildIndex(truth);

            // ゲノムの4割だけを1本で出したアセンブリ。
            // 自分の総延長を分母にすれば N50 は 8,000 になるが、
            // ゲノムサイズを分母にすると半分に届かないので 0 になる。
            var 一部 = this.WriteFasta("partial.fasta", truth[..8_000]);

            var 評価 = AssemblyScorer.Get_評価(一部, index, K, 深さ, truth.Length);

            Assert.NotNull(評価);
            Assert.Equal(0, 評価.A_NG50);
            Assert.InRange(評価.A_完全性, 0.35, 0.45);
        }
    }
}
