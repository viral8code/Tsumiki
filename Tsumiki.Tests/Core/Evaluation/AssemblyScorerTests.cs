using Tsumiki.Common;
using Tsumiki.Core.Evaluation;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// リファレンス無しでアセンブリの良さを測る評価器の検証
    /// </summary>
    /// <remarks>
    /// multi-k で複数のアセンブリから 1 つを選ぶには、リファレンスを使わずに
    /// 良し悪しを決められなければならない<br/>
    /// 連続性 (N50) だけで選ぶと
    /// 誤って繋いだものほど高く出るため、完全性と正確性を併せて見る必要がある<br/>
    /// ここではその「誤って繋いだものが落ちる」ことを主に固定する
    /// </remarks>
    public class AssemblyScorerTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public AssemblyScorerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_scorer_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
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
        /// この検証で使う k 長
        /// </summary>
        private const int K = 21;

        /// <summary>
        /// 深さ
        /// </summary>
        private const int 深さ = 20;

        /// <summary>
        /// 種を決めた乱数から塩基配列を作る
        /// </summary>
        /// <param name="length">作る長さ</param>
        /// <param name="seed">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        /// <summary>
        /// 与えた配列群から k-mer インデックスを作る
        /// </summary>
        /// <remarks>
        /// 深さは一律なので
        /// 単一コピー基準値は 深さ そのものになる
        /// </remarks>
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

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="name">ファイル名</param>
        /// <param name="sequences">書き出す配列</param>
        /// <returns>書き出したパス</returns>
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

        /// <summary>
        /// 名前を付けた配列を FASTA として書き出す
        /// </summary>
        /// <param name="name">ファイル名</param>
        /// <param name="entries">書き出す名前と配列</param>
        /// <returns>書き出したパス</returns>
        private string WriteFastaWithNames(string name, params (string A_名前, string A_配列)[] entries)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new FastaWriter(path);
            foreach (var (l_名前, l_配列) in entries)
            {
                writer.V_書き込み(l_名前, l_配列);
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
        /// 断片化しているが取りこぼしの無いアセンブリ
        /// </summary>
        /// <remarks>
        /// 完全性は満点のまま、
        /// 連続性だけが落ちること
        /// </remarks>
        [Fact]
        public void Score_FragmentedButComplete_LosesContiguityOnly()
        {
            var truth = RandomSequence(20_000, seed: 502);
            using var index = this.BuildIndex(truth);

            // k-1 塩基重ねて切ると、境界の k-mer も失われない
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
        /// これが評価器の存在意義
        /// </summary>
        /// <remarks>
        /// 反復配列を通り抜けて中間を飛ばした
        /// 誤アセンブリは、素の連続性では「改善」に見える (実際、過去に
        /// N50 が 99,974 から 199,945 へ伸びた誤アセンブリがあった)<br/>
        /// 飛ばした領域の k-mer が欠損として現れるため、完全性で見抜ける<br/>
        /// 重要なのは、連続性ではこのキメラのほうが上だという点である<br/>
        /// だからこそ選択規則は「まず完全性で足切りし、そのあとで連続性を見る」
        /// という順序でなければならない (掛け算にすると連続性の利得が勝ってしまう)
        /// </remarks>
        [Fact]
        public void Score_ChimeraThatSkipsSequence_ScoresBelowTheFragmentedButHonestAssembly()
        {
            // A-R-B-R-C
            // R は 2 回現れる反復配列
            var a = RandomSequence(8_000, seed: 511);
            var r = RandomSequence(300, seed: 512);
            var b = RandomSequence(8_000, seed: 513);
            var c = RandomSequence(8_000, seed: 514);
            var truth = a + r + b + r + c;

            using var index = this.BuildIndex(truth);

            // 正直な答え: R で切れているが、A も B も C も出ている
            var 正直 = this.WriteFasta("honest.fasta", a + r, r + b + r, r + c);
            // 誤アセンブリ: R を 1 回通り抜けて B を丸ごと飛ばした A-R-C
            var キメラ = this.WriteFasta("chimera.fasta", a + r + c);

            var 評価正直 = AssemblyScorer.Get_評価(正直, index, K, 深さ, truth.Length);
            var 評価キメラ = AssemblyScorer.Get_評価(キメラ, index, K, 深さ, truth.Length);

            Assert.NotNull(評価正直);
            Assert.NotNull(評価キメラ);

            // 飛ばした B のぶんだけキメラ側に欠損が出る
            Assert.True(評価キメラ.A_欠損延べ数 > 評価正直.A_欠損延べ数,
                $"chimera missing={評価キメラ.A_欠損延べ数}, honest missing={評価正直.A_欠損延べ数}");
            Assert.True(評価キメラ.A_完全性 < 評価正直.A_完全性);
            // 完全性の差は足切りの許容差を大きく超えていること
            Assert.True(評価正直.A_完全性 - 評価キメラ.A_完全性 > AssemblySelector.完全性の許容差,
                $"chimera completeness={評価キメラ.A_完全性:F3}, honest={評価正直.A_完全性:F3}");
            // 連続性だけを見るとキメラのほうが良く見えることを明示しておく
            Assert.True(評価キメラ.A_NG50 > 評価正直.A_NG50,
                "this test is only meaningful while the chimera looks better on contiguity alone");
        }

        /// <summary>
        /// 同じ配列を 2 回出した水増しは、正確性が落ちて総合点が下がること
        /// </summary>
        /// <remarks>
        /// 連続性 (NG50) はむしろ上がるため、この判定が無いと選んでしまう
        /// </remarks>
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
        /// 2 コピーの反復配列を 2 回出すのは正しい
        /// </summary>
        /// <remarks>
        /// 水増しと誤判定しないこと
        /// </remarks>
        [Fact]
        public void Score_TwoCopyRepeatEmittedTwice_IsNotPenalised()
        {
            var 単一 = RandomSequence(10_000, seed: 531);
            var 反復 = RandomSequence(500, seed: 532);
            // 反復が 2 回現れるゲノム
            var truth = 単一 + 反復 + RandomSequence(5_000, seed: 533) + 反復;

            using var index = this.BuildIndex(truth);

            var path = this.WriteFasta("repeat.fasta", truth);
            var 評価 = AssemblyScorer.Get_評価(path, index, K, 深さ, truth.Length);

            Assert.NotNull(評価);
            Assert.Equal(0, 評価.A_過剰延べ数);
            Assert.Equal(0, 評価.A_欠損延べ数);
        }

        /// <summary>
        /// NG50 は自分の総延長ではなく推定ゲノムサイズを分母にすること
        /// </summary>
        /// <remarks>
        /// 素の N50 だと「配列を落として短くなったアセンブリ」ほど有利になり、
        /// k を跨いだ比較に使えない
        /// </remarks>
        [Fact]
        public void Score_NG50_UsesTheGenomeSizeAsDenominator_NotTheAssemblyLength()
        {
            var truth = RandomSequence(20_000, seed: 541);
            using var index = this.BuildIndex(truth);

            // ゲノムの 4 割だけを 1 本で出したアセンブリ
            // 自分の総延長を分母にすれば N50 は 8,000 になるが、
            // ゲノムサイズを分母にすると半分に届かないので 0 になる
            var 一部 = this.WriteFasta("partial.fasta", truth[..8_000]);

            var 評価 = AssemblyScorer.Get_評価(一部, index, K, 深さ, truth.Length);

            Assert.NotNull(評価);
            Assert.Equal(0, 評価.A_NG50);
            Assert.InRange(評価.A_完全性, 0.35, 0.45);
        }

        /// <summary>
        /// 提案 H: 環状に閉じた complicon(ContigMaker が名前に"circular"を
        /// 付けたもの) は本数・総延長に数えること
        /// </summary>
        [Fact]
        public void Score_CircularContig_IsCountedInCircularStats()
        {
            var truth = RandomSequence(20_000, seed: 601);
            using var index = this.BuildIndex(truth);

            var path = this.WriteFastaWithNames("circular.fasta", ("NODE1_circular", truth));
            var 評価 = AssemblyScorer.Get_評価(path, index, K, 深さ, truth.Length);

            Assert.NotNull(評価);
            Assert.Equal(1, 評価.A_環状本数);
            Assert.Equal(1.0, 評価.A_環状化率, 6);
        }

        [Fact]
        public void Score_LinearContig_IsNotCountedAsCircular()
        {
            var truth = RandomSequence(20_000, seed: 602);
            using var index = this.BuildIndex(truth);

            var path = this.WriteFastaWithNames("linear.fasta", ("NODE1", truth));
            var 評価 = AssemblyScorer.Get_評価(path, index, K, 深さ, truth.Length);

            Assert.NotNull(評価);
            Assert.Equal(0, 評価.A_環状本数);
            Assert.Equal(0.0, 評価.A_環状化率);
        }

        /// <summary>
        /// 染色体よりはるかに小さいプラスミドでも、複製単位として数えられる
        /// 長さがあれば環状化率に反映されること
        /// </summary>
        /// <remarks>
        /// 閉じた複製単位は
        /// 「完全長を組み上げられた」ことの核心なので、連続性向けの物差しで
        /// 落としてはいけない
        /// </remarks>
        [Fact]
        public void Score_SmallCircularPlasmid_CountsTowardCircularFraction()
        {
            var chromosome = RandomSequence(20_000, seed: 603);
            var plasmid = RandomSequence(2_000, seed: 604);
            using var index = this.BuildIndex(chromosome, plasmid);

            var genomeSize = chromosome.Length + plasmid.Length;
            var path = this.WriteFastaWithNames(
                "with_plasmid.fasta", ("NODE1_circular", chromosome), ("NODE2_circular", plasmid));
            var 評価 = AssemblyScorer.Get_評価(path, index, K, 深さ, genomeSize);

            Assert.NotNull(評価);
            Assert.Equal(2, 評価.A_環状本数);
            Assert.InRange(評価.A_環状化率, 0.99, 1.0);
        }

        /// <summary>
        /// 評価に含める最小長 (500 bp) を下回る配列は、環状の目印が付いていても
        /// 数えない
        /// </summary>
        /// <remarks>
        /// de Bruijn グラフにはホモポリマーや短いタンデム反復に由来する
        /// 極小の閉路が多数あり、実データではこれが k あたり 10 本前後現れて
        /// 環状本数を埋め尽くした<br/>
        /// 環状本数は候補選択の最優先キーなので、
        /// 数えてしまうと k の選択がその雑音で決まる
        /// </remarks>
        [Fact]
        public void Score_TooShortSequences_AreNotCountedAtAll()
        {
            var chromosome = RandomSequence(20_000, seed: 605);
            using var index = this.BuildIndex(chromosome);

            var path = this.WriteFastaWithNames(
                "with_artefacts.fasta",
                ("NODE1_circular", chromosome),
                ("NODE2_circular", chromosome[..30]),
                ("NODE3_circular", chromosome[100..106]));
            var 評価 = AssemblyScorer.Get_評価(path, index, K, 深さ, chromosome.Length);

            Assert.NotNull(評価);
            Assert.Equal(1, 評価.A_環状本数);
            Assert.Equal(1, 評価.A_本数);
        }
    }
}
