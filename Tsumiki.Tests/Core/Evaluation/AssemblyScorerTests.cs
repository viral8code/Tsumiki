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
        /// <param name="p_長さ">作る長さ</param>
        /// <param name="p_シード">乱数の種</param>
        /// <returns>塩基配列</returns>
        private static string V_乱数配列(int p_長さ, int p_シード)
        {
            var l_rng = new Random(p_シード);
            return string.Concat(Enumerable.Range(0, p_長さ).Select(_ => "ACGT"[l_rng.Next(4)]));
        }

        /// <summary>
        /// 与えた配列群から k-mer インデックスを作る
        /// </summary>
        /// <param name="p_配列群">元になる配列</param>
        /// <returns>信頼できる k-mer 集合</returns>
        /// <remarks>
        /// 深さは一律なので
        /// 単一コピー基準値は 深さ そのものになる
        /// </remarks>
        private TrustedKmerIndex V_構築_索引(params string[] p_配列群)
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 1 };
            var l_index = new TrustedKmerIndex(this._tempDir);
            foreach (var l_seq in p_配列群)
            {
                var l_bytes = l_seq.Select(Util.Get_塩基ID).ToArray();
                for (var i = 0; i + K <= l_bytes.Length; i++)
                {
                    for (var l_rep = 0; l_rep < 深さ; l_rep++)
                    {
                        l_index.V_登録(l_bytes.AsSpan(i, K), p_ワーカー番号: 0);
                    }
                }
            }
            _ = l_index.V_カットオフ(p_カットオフ: 2);
            return l_index;
        }

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_配列群">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_Fasta(string p_ファイル名, params string[] p_配列群)
        {
            var l_path = Path.Combine(this._tempDir, p_ファイル名);
            using var l_writer = new FastaWriter(l_path);
            var l_id = 1;
            foreach (var l_seq in p_配列群)
            {
                l_writer.V_書き込み($"NODE{l_id++}", l_seq);
            }
            return l_path;
        }

        /// <summary>
        /// 名前を付けた配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_ファイル名">ファイル名</param>
        /// <param name="p_エントリ群">書き出す名前と配列</param>
        /// <returns>書き出したパス</returns>
        private string V_書き出し_Fasta_名前付き(string p_ファイル名, params (string A_名前, string A_配列)[] p_エントリ群)
        {
            var l_path = Path.Combine(this._tempDir, p_ファイル名);
            using var l_writer = new FastaWriter(l_path);
            foreach (var (l_名前, l_配列) in p_エントリ群)
            {
                l_writer.V_書き込み(l_名前, l_配列);
            }
            return l_path;
        }

        /// <summary>
        /// 完全なアセンブリは完全性も正確性も満点になることを確かめる
        /// </summary>
        [Fact]
        public void 完全なアセンブリは完全性も正確性も満点になる()
        {
            var l_truth = V_乱数配列(20_000, p_シード: 501);
            using var l_index = this.V_構築_索引(l_truth);

            var l_path = this.V_書き出し_Fasta("perfect.fasta", l_truth);
            var l_評価 = AssemblyScorer.Get_評価(l_path, l_index, K, 深さ, l_truth.Length);

            Assert.NotNull(l_評価);
            Assert.Equal(0, l_評価.A_欠損延べ数);
            Assert.Equal(0, l_評価.A_過剰延べ数);
            Assert.Equal(1.0, l_評価.A_完全性, 6);
            Assert.Equal(1.0, l_評価.A_正確性, 6);
        }

        /// <summary>
        /// 断片化しているが取りこぼしの無いアセンブリは連続性だけが落ちることを確かめる
        /// </summary>
        /// <remarks>
        /// 完全性は満点のまま、
        /// 連続性だけが落ちること
        /// </remarks>
        [Fact]
        public void 断片化しているが取りこぼしの無いアセンブリは連続性だけが落ちる()
        {
            var l_truth = V_乱数配列(20_000, p_シード: 502);
            using var l_index = this.V_構築_索引(l_truth);

            // k-1 塩基重ねて切ると、境界の k-mer も失われない
            var l_断片 = new List<string>();
            for (var i = 0; i < l_truth.Length; i += 4_000)
            {
                var l_終端 = Math.Min(l_truth.Length, i + 4_000 + K - 1);
                l_断片.Add(l_truth[i..l_終端]);
            }

            var l_一本 = this.V_書き出し_Fasta("whole.fasta", l_truth);
            var l_断片化 = this.V_書き出し_Fasta("fragmented.fasta", [.. l_断片]);

            var l_評価一本 = AssemblyScorer.Get_評価(l_一本, l_index, K, 深さ, l_truth.Length);
            var l_評価断片 = AssemblyScorer.Get_評価(l_断片化, l_index, K, 深さ, l_truth.Length);

            Assert.NotNull(l_評価一本);
            Assert.NotNull(l_評価断片);
            Assert.Equal(0, l_評価断片.A_欠損延べ数);
            Assert.True(l_評価断片.A_NG50 < l_評価一本.A_NG50);
        }

        /// <summary>
        /// 配列を飛ばしたキメラは断片化しているが正直なアセンブリより低く評価されることを確かめる
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
        public void 配列を飛ばしたキメラは断片化しているが正直なアセンブリより低く評価される()
        {
            // A-R-B-R-C
            // R は 2 回現れる反復配列
            var l_a = V_乱数配列(8_000, p_シード: 511);
            var l_r = V_乱数配列(300, p_シード: 512);
            var l_b = V_乱数配列(8_000, p_シード: 513);
            var l_c = V_乱数配列(8_000, p_シード: 514);
            var l_truth = l_a + l_r + l_b + l_r + l_c;

            using var l_index = this.V_構築_索引(l_truth);

            // 正直な答え: R で切れているが、A も B も C も出ている
            var l_正直 = this.V_書き出し_Fasta("honest.fasta", l_a + l_r, l_r + l_b + l_r, l_r + l_c);
            // 誤アセンブリ: R を 1 回通り抜けて B を丸ごと飛ばした A-R-C
            var l_キメラ = this.V_書き出し_Fasta("chimera.fasta", l_a + l_r + l_c);

            var l_評価正直 = AssemblyScorer.Get_評価(l_正直, l_index, K, 深さ, l_truth.Length);
            var l_評価キメラ = AssemblyScorer.Get_評価(l_キメラ, l_index, K, 深さ, l_truth.Length);

            Assert.NotNull(l_評価正直);
            Assert.NotNull(l_評価キメラ);

            // 飛ばした B のぶんだけキメラ側に欠損が出る
            Assert.True(l_評価キメラ.A_欠損延べ数 > l_評価正直.A_欠損延べ数,
                $"chimera missing={l_評価キメラ.A_欠損延べ数}, honest missing={l_評価正直.A_欠損延べ数}");
            Assert.True(l_評価キメラ.A_完全性 < l_評価正直.A_完全性);
            // 完全性の差は足切りの許容差を大きく超えていること
            Assert.True(l_評価正直.A_完全性 - l_評価キメラ.A_完全性 > AssemblySelector.完全性の許容差,
                $"chimera completeness={l_評価キメラ.A_完全性:F3}, honest={l_評価正直.A_完全性:F3}");
            // 連続性だけを見るとキメラのほうが良く見えることを明示しておく
            Assert.True(l_評価キメラ.A_NG50 > l_評価正直.A_NG50,
                "this test is only meaningful while the chimera looks better on contiguity alone");
        }

        /// <summary>
        /// 配列を 2 回出した水増しは正確性が落ちて総合点が下がることを確かめる
        /// </summary>
        /// <remarks>
        /// 連続性 (NG50) はむしろ上がるため、この判定が無いと選んでしまう
        /// </remarks>
        [Fact]
        public void 配列を2回出した水増しは正確性が落ちて総合点が下がる()
        {
            var l_truth = V_乱数配列(20_000, p_シード: 521);
            using var l_index = this.V_構築_索引(l_truth);

            var l_正常 = this.V_書き出し_Fasta("single.fasta", l_truth);
            var l_水増し = this.V_書き出し_Fasta("doubled.fasta", l_truth, l_truth);

            var l_評価正常 = AssemblyScorer.Get_評価(l_正常, l_index, K, 深さ, l_truth.Length);
            var l_評価水増し = AssemblyScorer.Get_評価(l_水増し, l_index, K, 深さ, l_truth.Length);

            Assert.NotNull(l_評価正常);
            Assert.NotNull(l_評価水増し);
            Assert.True(l_評価水増し.A_過剰延べ数 > 0);
            Assert.True(l_評価水増し.A_正確性 < l_評価正常.A_正確性);
        }

        /// <summary>
        /// 2 コピーの反復配列を 2 回出しても減点されないことを確かめる
        /// </summary>
        /// <remarks>
        /// 水増しと誤判定しないこと
        /// </remarks>
        [Fact]
        public void 二コピーの反復配列を2回出しても減点されない()
        {
            var l_単一 = V_乱数配列(10_000, p_シード: 531);
            var l_反復 = V_乱数配列(500, p_シード: 532);
            // 反復が 2 回現れるゲノム
            var l_truth = l_単一 + l_反復 + V_乱数配列(5_000, p_シード: 533) + l_反復;

            using var l_index = this.V_構築_索引(l_truth);

            var l_path = this.V_書き出し_Fasta("repeat.fasta", l_truth);
            var l_評価 = AssemblyScorer.Get_評価(l_path, l_index, K, 深さ, l_truth.Length);

            Assert.NotNull(l_評価);
            Assert.Equal(0, l_評価.A_過剰延べ数);
            Assert.Equal(0, l_評価.A_欠損延べ数);
        }

        /// <summary>
        /// NG50 は自分の総延長ではなく推定ゲノムサイズを分母にすることを確かめる
        /// </summary>
        /// <remarks>
        /// 素の N50 だと「配列を落として短くなったアセンブリ」ほど有利になり、
        /// k を跨いだ比較に使えない
        /// </remarks>
        [Fact]
        public void NG50は自分の総延長ではなく推定ゲノムサイズを分母にする()
        {
            var l_truth = V_乱数配列(20_000, p_シード: 541);
            using var l_index = this.V_構築_索引(l_truth);

            // ゲノムの 4 割だけを 1 本で出したアセンブリ
            // 自分の総延長を分母にすれば N50 は 8,000 になるが、
            // ゲノムサイズを分母にすると半分に届かないので 0 になる
            var l_一部 = this.V_書き出し_Fasta("partial.fasta", l_truth[..8_000]);

            var l_評価 = AssemblyScorer.Get_評価(l_一部, l_index, K, 深さ, l_truth.Length);

            Assert.NotNull(l_評価);
            Assert.Equal(0, l_評価.A_NG50);
            Assert.InRange(l_評価.A_完全性, 0.35, 0.45);
        }

        /// <summary>
        /// 環状に閉じた complicon(ContigMaker が名前に "circular" を
        /// 付けたもの)は本数・総延長に数えられることを確かめる
        /// </summary>
        [Fact]
        public void 環状contigは環状本数に数えられる()
        {
            var l_truth = V_乱数配列(20_000, p_シード: 601);
            using var l_index = this.V_構築_索引(l_truth);

            var l_path = this.V_書き出し_Fasta_名前付き("circular.fasta", ("NODE1_circular", l_truth));
            var l_評価 = AssemblyScorer.Get_評価(l_path, l_index, K, 深さ, l_truth.Length);

            Assert.NotNull(l_評価);
            Assert.Equal(1, l_評価.A_環状本数);
            Assert.Equal(1.0, l_評価.A_環状化率, 6);
        }

        /// <summary>
        /// 線状の contig は環状として数えられないことを確かめる
        /// </summary>
        [Fact]
        public void 線状contigは環状として数えられない()
        {
            var l_truth = V_乱数配列(20_000, p_シード: 602);
            using var l_index = this.V_構築_索引(l_truth);

            var l_path = this.V_書き出し_Fasta_名前付き("linear.fasta", ("NODE1", l_truth));
            var l_評価 = AssemblyScorer.Get_評価(l_path, l_index, K, 深さ, l_truth.Length);

            Assert.NotNull(l_評価);
            Assert.Equal(0, l_評価.A_環状本数);
            Assert.Equal(0.0, l_評価.A_環状化率);
        }

        /// <summary>
        /// 染色体よりはるかに小さいプラスミドでも、複製単位として数えられる
        /// 長さがあれば環状化率に反映されることを確かめる
        /// </summary>
        /// <remarks>
        /// 閉じた複製単位は
        /// 「完全長を組み上げられた」ことの核心なので、連続性向けの物差しで
        /// 落としてはいけない
        /// </remarks>
        [Fact]
        public void 小さい環状プラスミドも環状化率に反映される()
        {
            var l_chromosome = V_乱数配列(20_000, p_シード: 603);
            var l_plasmid = V_乱数配列(2_000, p_シード: 604);
            using var l_index = this.V_構築_索引(l_chromosome, l_plasmid);

            var l_genomeSize = l_chromosome.Length + l_plasmid.Length;
            var l_path = this.V_書き出し_Fasta_名前付き(
                "with_plasmid.fasta", ("NODE1_circular", l_chromosome), ("NODE2_circular", l_plasmid));
            var l_評価 = AssemblyScorer.Get_評価(l_path, l_index, K, 深さ, l_genomeSize);

            Assert.NotNull(l_評価);
            Assert.Equal(2, l_評価.A_環状本数);
            Assert.InRange(l_評価.A_環状化率, 0.99, 1.0);
        }

        /// <summary>
        /// 評価に含める最小長 (500 bp) を下回る配列は、環状の目印が付いていても
        /// 数えられないことを確かめる
        /// </summary>
        /// <remarks>
        /// de Bruijn グラフにはホモポリマーや短いタンデム反復に由来する
        /// 極小の閉路が多数あり、実データではこれが k あたり 10 本前後現れて
        /// 環状本数を埋め尽くした<br/>
        /// 環状本数は候補選択の最優先キーなので、
        /// 数えてしまうと k の選択がその雑音で決まる
        /// </remarks>
        [Fact]
        public void 評価に含める最小長を下回る配列は数えられない()
        {
            var l_chromosome = V_乱数配列(20_000, p_シード: 605);
            using var l_index = this.V_構築_索引(l_chromosome);

            var l_path = this.V_書き出し_Fasta_名前付き(
                "with_artefacts.fasta",
                ("NODE1_circular", l_chromosome),
                ("NODE2_circular", l_chromosome[..30]),
                ("NODE3_circular", l_chromosome[100..106]));
            var l_評価 = AssemblyScorer.Get_評価(l_path, l_index, K, 深さ, l_chromosome.Length);

            Assert.NotNull(l_評価);
            Assert.Equal(1, l_評価.A_環状本数);
            Assert.Equal(1, l_評価.A_本数);
        }
    }
}
