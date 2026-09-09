using Tsumiki.Common;
using Tsumiki.Core.Scaffolding;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// GapFiller が埋められなかったスキャフォールドのギャップを、その両端に
    /// 実際にマップされた局所リードだけで再アセンブリして埋める処理の検証
    /// </summary>
    /// <remarks>
    /// AssemblyMerger(-mg) と違い、他の k の「既に確定した結論」を持ち込むの
    /// ではなく、生リードから新しく証拠を集める<br/>
    /// ここでは生の FASTQ を
    /// 直接与えて、局所アセンブリだけでギャップが埋まる/埋まらないことを
    /// 検証する (GapFiller 側は経由しない)
    /// </remarks>
    public class LocalAssemblerTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public LocalAssemblerTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_localasm_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
            // LocalAssembler は内部で TrustedKmerIndex/KmerKey を使うため、
            // 現在の実行時引数の k 長を、テストで使うk(21) に合わせておく
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 21, A_スレッド数 = 1 };
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
        /// スキャフォールドを FASTA として書き出す
        /// </summary>
        /// <param name="name">ファイル名</param>
        /// <param name="sequence">書き出す配列</param>
        /// <returns>書き出したパス</returns>
        private string WriteScaffold(string name, string sequence)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new FastaWriter(path);
            writer.V_書き込み("SCAFFOLD1", sequence);
            return path;
        }

        /// <summary>
        /// ゲノムを覆うリードを FASTQ として書き出す
        /// </summary>
        /// <param name="name">ファイル名</param>
        /// <param name="genomes">元になるゲノム</param>
        /// <param name="readLength">リード長</param>
        /// <returns>書き出したパス</returns>
        private string WriteReads(string name, IEnumerable<string> genomes, int readLength)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new StreamWriter(path);
            var idx = 0;
            foreach (var genome in genomes)
            {
                for (var i = 0; i + readLength <= genome.Length; i++)
                {
                    writer.WriteLine($"@r{idx}");
                    writer.WriteLine(genome.Substring(i, readLength));
                    writer.WriteLine("+");
                    writer.WriteLine(new string('I', readLength));
                    idx++;
                }
            }
            return path;
        }

        /// <summary>
        /// FASTA から 1 本だけの配列を読み込む
        /// </summary>
        /// <param name="path">読み込むパス</param>
        /// <returns>読み込んだ配列</returns>
        private static string ReadSingleSequence(string path)
        {
            using var reader = new FastaReader(path);
            Assert.True(reader.Get_続きがあるか());
            return reader.Get_次の配列().A_配列;
        }

        [Fact]
        public void FillGaps_ReadsActuallySpanTheGap_RestoresTheTrueSequence()
        {
            const int k = 21;
            var prefix = RandomSequence(300, seed: 1);
            var fill = RandomSequence(60, seed: 2);
            var suffix = RandomSequence(300, seed: 3);
            var truth = prefix + fill + suffix;

            var scaffoldPath = this.WriteScaffold("scaffold.fasta", prefix + new string('N', fill.Length) + suffix);
            var readsPath = this.WriteReads("reads.fq", [truth], readLength: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(1, stats.A_対象ギャップ数);
            Assert.Equal(1, stats.A_埋めたギャップ数);
            Assert.Equal(fill.Length, stats.A_埋めた塩基数);
            Assert.Equal(truth, ReadSingleSequence(scaffoldPath));
        }

        [Fact]
        public void FillGaps_NoLocalReadMapsNearEitherAnchor_LeavesTheGapAsN()
        {
            const int k = 21;
            var prefix = RandomSequence(300, seed: 10);
            var fill = RandomSequence(60, seed: 11);
            var suffix = RandomSequence(300, seed: 12);

            var scaffoldPath = this.WriteScaffold("scaffold_noreads.fasta", prefix + new string('N', fill.Length) + suffix);
            // まったく無関係な配列からリードを取る
            var unrelated = RandomSequence(500, seed: 999);
            var readsPath = this.WriteReads("reads_unrelated.fq", [unrelated], readLength: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(1, stats.A_対象ギャップ数);
            Assert.Equal(0, stats.A_埋めたギャップ数);
            Assert.Equal(1, stats.A_局所リードが集まらなかった数);
            Assert.Contains('N', ReadSingleSequence(scaffoldPath));
        }

        /// <summary>
        /// アンカー付近のリードはあるが、両者を橋渡しする配列 (ギャップの中身) を
        /// 読んだリードが無い場合、経路が繋がらないので埋められない
        /// </summary>
        [Fact]
        public void FillGaps_LocalReadsNeverBridgeTheGap_LeavesItAsNWithNoPathVerdict()
        {
            const int k = 21;
            var prefix = RandomSequence(300, seed: 20);
            var fill = RandomSequence(60, seed: 21);
            var suffix = RandomSequence(300, seed: 22);

            var scaffoldPath = this.WriteScaffold("scaffold_nopath.fasta", prefix + new string('N', fill.Length) + suffix);
            // prefix と suffix それぞれの内部だけを読んだリード (橋渡しは無い)
            var readsPath = this.WriteReads("reads_nopath.fq", [prefix, suffix], readLength: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(1, stats.A_対象ギャップ数);
            Assert.Equal(0, stats.A_埋めたギャップ数);
            Assert.Equal(1, stats.A_到達できなかった数);
            Assert.Contains('N', ReadSingleSequence(scaffoldPath));
        }

        /// <summary>
        /// 橋渡しの配列が 2 通りとも読まれている場合、どちらが正しいか
        /// 決められないので N のまま残す (誤った配列で埋めるより安全)
        /// </summary>
        [Fact]
        public void FillGaps_TwoEquallySupportedFills_LeavesItAsNRatherThanGuessing()
        {
            const int k = 21;
            var prefix = RandomSequence(300, seed: 30);
            var suffix = RandomSequence(300, seed: 31);
            var fillA = RandomSequence(60, seed: 32);
            var fillB = RandomSequence(60, seed: 33);

            var scaffoldPath = this.WriteScaffold("scaffold_ambiguous.fasta", prefix + new string('N', fillA.Length) + suffix);
            var readsPath = this.WriteReads(
                "reads_ambiguous.fq", [prefix + fillA + suffix, prefix + fillB + suffix], readLength: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(1, stats.A_対象ギャップ数);
            Assert.Equal(0, stats.A_埋めたギャップ数);
            Assert.Equal(1, stats.A_一意に定まらなかった数);
            Assert.Contains('N', ReadSingleSequence(scaffoldPath));
        }

        [Fact]
        public void FillGaps_NoGapsInScaffold_ReportsNoEligibleGap()
        {
            const int k = 21;
            var truth = RandomSequence(300, seed: 40);
            var scaffoldPath = this.WriteScaffold("scaffold_nogap.fasta", truth);
            var readsPath = this.WriteReads("reads_nogap.fq", [truth], readLength: 100);

            var stats = LocalAssembler.V_充填_ギャップ(scaffoldPath, readsPath, string.Empty, k, this._tempDir);

            Assert.Equal(0, stats.A_対象ギャップ数);
            Assert.Equal(truth, ReadSingleSequence(scaffoldPath));
        }
    }
}
