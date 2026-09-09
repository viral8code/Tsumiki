using Tsumiki.IO;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// 最終成果物から、リード長より短い配列を落とす処理の検証<br/>
    /// リード 1 本に収まる配列はリードそのもの以上の情報を運ばず、
    /// この帯には繋ぎ間違えた断片が集まりやすい
    /// </summary>
    public class ShortSequenceFilterTests : IDisposable
    {
        /// <summary>
        /// 一時ディレクトリのパス
        /// </summary>
        private readonly string _tempDir;

        public ShortSequenceFilterTests()
        {
            this._tempDir = Path.Combine(Path.GetTempPath(), "tsumiki_shortfilter_tests_" + Guid.NewGuid().ToString("N"));
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
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// 配列を FASTA として書き出す
        /// </summary>
        /// <param name="p_全件">書き出す名前と配列</param>
        /// <returns>書き出したパス</returns>
        private string Get_書き出し(params (string A_ID, string A_配列)[] p_全件)
        {
            var l_パス = Path.Combine(this._tempDir, "scaffolds.fasta");
            using (var l_書き込み = new FastaWriter(l_パス))
            {
                foreach (var (l_ID, l_配列) in p_全件)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }
            return l_パス;
        }

        [Fact]
        public void V_除外_短い配列_DropsOnlyWhatIsShorterThanTheRead()
        {
            var l_パス = this.Get_書き出し(
                ("SCAFFOLD1", new string('A', 300)),
                ("SCAFFOLD2", new string('C', 150)),
                ("SCAFFOLD3", new string('G', 149)),
                ("SCAFFOLD4", new string('T', 1)));

            Tsumiki.Program.V_除外_短い配列(l_パス, 150);

            var l_残り = FastaReader.Get_全エントリ(l_パス);
            Assert.Equal(["SCAFFOLD1", "SCAFFOLD2"], l_残り.Select(x => x.A_ID));
            Assert.Equal(300, l_残り[0].A_配列.Length);
            Assert.Equal(150, l_残り[1].A_配列.Length);
        }

        [Fact]
        public void V_除外_短い配列_LeavesTheFileAloneWhenNothingIsShort()
        {
            var l_パス = this.Get_書き出し(
                ("SCAFFOLD1", new string('A', 300)),
                ("SCAFFOLD2", new string('C', 200)));
            var l_元 = File.ReadAllBytes(l_パス);

            Tsumiki.Program.V_除外_短い配列(l_パス, 150);

            Assert.Equal(l_元, File.ReadAllBytes(l_パス));
        }

        [Fact]
        public void V_除外_短い配列_DoesNothingWhenTheReadLengthIsUnknown()
        {
            var l_パス = this.Get_書き出し(("SCAFFOLD1", new string('A', 10)));
            var l_元 = File.ReadAllBytes(l_パス);

            Tsumiki.Program.V_除外_短い配列(l_パス, null);

            Assert.Equal(l_元, File.ReadAllBytes(l_パス));
        }
    
        [Fact]
        public void V_削除_中間ファイル_KeepsTheFinalProductsAndTheLog()
        {
            foreach (var l_名前 in new[]
            {
                "unitigs.fasta", "contigs.fasta", "scaffolds.fasta", "assembly.gfa",
                "assembly.report.json", "assembly.ambiguous.tsv", "Tsumiki.log",
                "corrected.1.fq", "preprocessed.1.fq", "polished.fasta", "merged_scaffolds.fasta",
            })
            {
                File.WriteAllText(Path.Combine(this._tempDir, l_名前), l_名前);
            }
            _ = Directory.CreateDirectory(Path.Combine(this._tempDir, "k93"));
            File.WriteAllText(Path.Combine(this._tempDir, "k93", "unitigs.fasta"), "x");

            Tsumiki.Program.V_削除_中間ファイル(this._tempDir);

            var l_残り = Directory.EnumerateFiles(this._tempDir)
                .Select(Path.GetFileName).Order(StringComparer.Ordinal).ToList();
            Assert.Equal(
                [
                    "Tsumiki.log", "assembly.ambiguous.tsv", "assembly.gfa", "assembly.report.json",
                    "contigs.fasta", "scaffolds.fasta", "unitigs.fasta",
                ],
                l_残り);
            Assert.Empty(Directory.EnumerateDirectories(this._tempDir));
        }
}
}
