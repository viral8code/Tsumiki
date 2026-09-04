using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// リード長が k より短いリードが混ざっていても処理が破綻しないことを固定する。
    ///
    /// トリミング済みのデータではリード長がばらつく。GAGE-B の
    /// R. sphaeroides MiSeq(trimmed)では 755,847 本のうち 8% 以上が
    /// k=63 未満で、最短は 19bp だった。マッピング側に長さの判定が無く、
    /// 19bp のリードに対して添字 62 までアクセスして例外になっていた。
    ///
    /// しかもその例外はワーカースレッドの中で起き、キューが満杯になった
    /// プロデューサーが永久に待ち続けたため、ログも例外も出ないまま
    /// 2時間以上プロセスが停止した。長さの判定と、
    /// ワーカーの例外を伝える仕組み(ReadPipelineTests)の両方が要る。
    /// </summary>
    public class VariableLengthReadTests : IDisposable
    {
        private readonly string _tempDir;

        public VariableLengthReadTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_varlen_tests_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(this._tempDir);
        }

        public void Dispose()
        {
            if (Directory.Exists(this._tempDir))
            {
                Directory.Delete(this._tempDir, recursive: true);
            }
        }

        private const int K = 31;

        private static string RandomSequence(int length, int seed)
        {
            var rng = new Random(seed);
            return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]));
        }

        private string WriteFastq(string name, IEnumerable<(string A_ID, string A_配列)> reads)
        {
            var path = Path.Combine(this._tempDir, name);
            using var writer = new StreamWriter(path);
            foreach (var (id, seq) in reads)
            {
                writer.WriteLine($"@{id}");
                writer.WriteLine(seq);
                writer.WriteLine("+");
                writer.WriteLine(new string('I', seq.Length)); // Q40相当
            }
            return path;
        }

        /// <summary>
        /// k より短いリードと十分長いリードが混ざったペアエンド入力。
        /// 短いリードは黙って読み飛ばされ、長いリード由来の隣接だけが残ること。
        /// </summary>
        [Fact]
        public void MapPairedReads_ReadsShorterThanK_AreSkippedWithoutFailing()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 4 };

            var unitigSeq = RandomSequence(600, seed: 987);
            var unitigsPath = Path.Combine(this._tempDir, "unitigs.fasta");
            File.WriteAllText(unitigsPath, $">1\n{unitigSeq}\n");

            // 19bp(最短の実例と同じ長さ)から 200bp まで、k をまたぐ長さを混ぜる。
            var lengths = new[] { 19, 30, K - 1, K, K + 1, 120, 200 };
            var reads1 = new List<(string, string)>();
            var reads2 = new List<(string, string)>();
            for (var i = 0; i < lengths.Length; i++)
            {
                var length = lengths[i];
                reads1.Add(($"pair{i}/1", unitigSeq[..length]));
                reads2.Add(($"pair{i}/2", Util.V_逆相補(unitigSeq[^length..])));
            }

            var path1 = this.WriteFastq("short.1.fq", reads1);
            var path2 = this.WriteFastq("short.2.fq", reads2);

            var contigMaker = new ContigMaker(unitigsPath);

            // 例外を投げずに完走すること。対策前はここで
            // IndexOutOfRangeException がワーカー内で起き、
            // そのままハングしていた。
            contigMaker.V_マッピング_ペアリード(path1, path2);

            // k 以上のリードからは標本が取れていること
            // (短いリードのせいで全部落ちてしまっていないことの確認)。
            Assert.NotEmpty(contigMaker.A_インサートサイズ標本);
        }

        [Fact]
        public void MapSingleReads_ReadsShorterThanK_AreSkippedWithoutFailing()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 4 };

            var unitigSeq = RandomSequence(400, seed: 654);
            var unitigsPath = Path.Combine(this._tempDir, "unitigs_single.fasta");
            File.WriteAllText(unitigsPath, $">1\n{unitigSeq}\n");

            var reads = new List<(string, string)>();
            for (var i = 0; i < 50; i++)
            {
                // 半分を k 未満にする。
                var length = i % 2 == 0 ? 19 : 150;
                reads.Add(($"read{i}", unitigSeq[..length]));
            }
            var path = this.WriteFastq("short_single.fq", reads);

            var contigMaker = new ContigMaker(unitigsPath);
            contigMaker.V_マッピング_リード(path);
        }

        /// <summary>
        /// すべてのリードが k 未満でも、例外にならず単に何も得られないこと。
        /// </summary>
        [Fact]
        public void MapPairedReads_EveryReadShorterThanK_CompletesWithNoSamples()
        {
            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = K, A_スレッド数 = 4 };

            var unitigSeq = RandomSequence(400, seed: 321);
            var unitigsPath = Path.Combine(this._tempDir, "unitigs_allshort.fasta");
            File.WriteAllText(unitigsPath, $">1\n{unitigSeq}\n");

            var reads1 = Enumerable.Range(0, 40).Select(i => ($"pair{i}/1", unitigSeq[..19]));
            var reads2 = Enumerable.Range(0, 40).Select(i => ($"pair{i}/2", unitigSeq[..20]));

            var path1 = this.WriteFastq("allshort.1.fq", reads1);
            var path2 = this.WriteFastq("allshort.2.fq", reads2);

            var contigMaker = new ContigMaker(unitigsPath);
            contigMaker.V_マッピング_ペアリード(path1, path2);

            Assert.Empty(contigMaker.A_インサートサイズ標本);
        }
    }
}
