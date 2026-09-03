using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// ペアエンドから推定される「インサートサイズ」が、リードに挟まれた内側の
    /// 未読区間ではなく、真のフラグメント長(左リードの5'端から右リードの
    /// 3'端まで)の単位になっていることを検証する。
    ///
    /// 実データ(150bpリード・IS350ライブラリ)で同一unitig由来サンプルの
    /// 中央値が58と報告されていた。リード長150bpより短いフラグメントは
    /// physically ありえないため、これは単位の取り違えを示していた。
    /// 内側距離58に両リード長を足すと358となりライブラリ名と一致する。
    /// この取り違えはギャップ長推定(ギャップ = インサートサイズ - 既知長)にも
    /// そのまま伝播するため、単位を明示的に固定しておく。
    /// </summary>
    public class InsertSizeEstimationTests : IDisposable
    {
        private readonly string _tempDir;

        public InsertSizeEstimationTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_insertsize_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        /// <summary>
        /// 決定的な擬似乱数で非反復的な塩基配列を作る。k=21 では
        /// この長さの乱数配列に重複k-merが現れる確率は無視できる。
        /// </summary>
        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        private static void WriteFastq(string path, IEnumerable<(string Id, string Seq)> reads)
        {
            using var writer = new StreamWriter(path);
            foreach (var (id, seq) in reads)
            {
                writer.WriteLine($"@{id}");
                writer.WriteLine(seq);
                writer.WriteLine("+");
                writer.WriteLine(new string('I', seq.Length)); // Q40相当
            }
        }

        [Fact]
        public void SameUnitigSamples_MeasureFullFragmentLength_NotTheInnerDistance()
        {
            const int k = 21;
            const int unitigLength = 600;
            const int readLength = 50;
            const int fragmentStart = 100;
            const int trueFragmentLength = 350;

            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = 1 };

            var unitigSeq = RandomSequence(unitigLength, seed: 12345);
            var unitigsPath = Path.Combine(this._tempDir, "unitigs.fasta");
            File.WriteAllText(unitigsPath, $">1\n{unitigSeq}\n");

            // FR配置: read1 はフラグメント左端から順鎖方向、
            // read2 はフラグメント右端から逆鎖方向に読まれる。
            var read1 = unitigSeq.Substring(fragmentStart, readLength);
            var read2 = Util.ReverseComprement(
                unitigSeq.Substring(fragmentStart + trueFragmentLength - readLength, readLength));

            var path1 = Path.Combine(this._tempDir, "r1.fq");
            var path2 = Path.Combine(this._tempDir, "r2.fq");
            // 中央値を安定させるため同一ペアを複数本入れる。
            var pairs = Enumerable.Range(0, 5).ToList();
            WriteFastq(path1, pairs.Select(i => ($"pair{i}/1", read1)));
            WriteFastq(path2, pairs.Select(i => ($"pair{i}/2", read2)));

            var contigMaker = new ContigMaker(unitigsPath);
            contigMaker.MappingPairedReads(path1, path2);

            Assert.NotEmpty(contigMaker.SameUnitigInsertSizeSamples);
            // 内側距離(= 350 - 50 - 50 = 250)ではなく、フラグメント長 350 が
            // 得られなければならない。
            Assert.All(
                contigMaker.SameUnitigInsertSizeSamples,
                sample => Assert.Equal(trueFragmentLength, sample));
        }

        /// <summary>
        /// フラグメント長を変えたときに推定値が同じだけ動くこと(定数ぶんの
        /// ずれではなく、単位そのものが一致していること)を確認する。
        /// </summary>
        [Theory]
        [InlineData(200)]
        [InlineData(350)]
        [InlineData(500)]
        public void SameUnitigSamples_TrackTheActualFragmentLength(int trueFragmentLength)
        {
            const int k = 21;
            const int unitigLength = 900;
            const int readLength = 50;
            const int fragmentStart = 120;

            ConfigurationManager.Arguments = new Parameters { Kmer = k, ThreadCount = 1 };

            var unitigSeq = RandomSequence(unitigLength, seed: 777);
            var unitigsPath = Path.Combine(this._tempDir, $"unitigs_{trueFragmentLength}.fasta");
            File.WriteAllText(unitigsPath, $">1\n{unitigSeq}\n");

            var read1 = unitigSeq.Substring(fragmentStart, readLength);
            var read2 = Util.ReverseComprement(
                unitigSeq.Substring(fragmentStart + trueFragmentLength - readLength, readLength));

            var path1 = Path.Combine(this._tempDir, $"r1_{trueFragmentLength}.fq");
            var path2 = Path.Combine(this._tempDir, $"r2_{trueFragmentLength}.fq");
            WriteFastq(path1, [("pair/1", read1)]);
            WriteFastq(path2, [("pair/2", read2)]);

            var contigMaker = new ContigMaker(unitigsPath);
            contigMaker.MappingPairedReads(path1, path2);

            Assert.NotEmpty(contigMaker.SameUnitigInsertSizeSamples);
            Assert.All(
                contigMaker.SameUnitigInsertSizeSamples,
                sample => Assert.Equal(trueFragmentLength, sample));
        }
    }
}
