using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;

namespace Tsumiki.Tests.Core
{
    /// <summary>
    /// tools/simulate_reads.py が出力した合成データ(正解のerrors.tsv付き)に対して
    /// ErrorCorrector.V_訂正_リードファイル を実際に走らせ、注入したエラーのうち
    /// 何割を正しく真の塩基へ戻せたか(recall)、逆に正しかった塩基を
    /// 誤って書き換えてしまった割合(誤訂正率)を測定する検証用テスト。
    /// 合成データが存在しない場合はスキップする(通常のCI/dotnet testの対象外、
    /// 手動でtools/simulate_reads.pyを実行した後に手動で実行する想定)。
    /// </summary>
    public class ErrorCorrectorGroundTruthValidation
    {
        // Bash tool経由(Git Bash/MSYS)で python tools/simulate_reads.py --out-dir /tmp/tsumiki_synth
        // を実行した場合の実際の出力先(MSYSが/tmpをこのWindowsパスへ解決する)。
        // .NET のファイルAPIはMSYSのパス変換を経由しないため、Windows形式で直接指定する。
        private static readonly string SynthDir = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "Temp", "tsumiki_synth");

        [Fact]
        public void MeasureCorrectionAccuracyAgainstGroundTruth()
        {
            var refPath = Path.Combine(SynthDir, "reference.fasta");
            var read1Path = Path.Combine(SynthDir, "reads.1.fq");
            var read2Path = Path.Combine(SynthDir, "reads.2.fq");
            var errorsPath = Path.Combine(SynthDir, "errors.tsv");
            if (!File.Exists(refPath) || !File.Exists(read1Path) || !File.Exists(errorsPath))
            {
                return; // 合成データ未生成。tools/simulate_reads.py --out-dir /tmp/tsumiki_synth で生成してから実行する。
            }

            ConfigurationManager.A_実行時引数 = new Parameters { A_k長 = 31, A_kmerカットオフ = 2, A_スレッド数 = 8 };

            var outDir = Path.Combine(Path.GetTempPath(), "tsumiki_ec_validation_" + Guid.NewGuid().ToString("N"));
            _ = Directory.CreateDirectory(outDir);
            var corrected1 = Path.Combine(outDir, "corrected.1.fq");
            var corrected2 = Path.Combine(outDir, "corrected.2.fq");

            try
            {
                ErrorCorrector.V_訂正_リードファイル(read1Path, read2Path, outDir, corrected1, corrected2);

                // read_id -> mate -> position -> true_base (注入されたエラーの正解)
                var trueErrors = new Dictionary<(string ReadId, int Mate, int Position), char>();
                foreach (var line in File.ReadLines(errorsPath).Skip(1))
                {
                    var parts = line.Split('\t');
                    var readId = parts[0];
                    var mate = int.Parse(parts[1]);
                    var position = int.Parse(parts[2]);
                    var trueBase = parts[3][0];
                    trueErrors[(readId, mate, position)] = trueBase;
                }

                var originalReads = LoadReadsById(read1Path, 1);
                foreach (var kv in LoadReadsById(read2Path, 2))
                {
                    originalReads[kv.Key] = kv.Value;
                }

                var fixedCount = 0;
                var stillWrongCount = 0;
                var newlyWrongCount = 0; // 元々正しかった塩基を誤って書き換えてしまった数
                var totalChangedPositions = 0;

                ValidateFile(corrected1, 1, originalReads, trueErrors, ref fixedCount, ref stillWrongCount, ref newlyWrongCount, ref totalChangedPositions);
                ValidateFile(corrected2, 2, originalReads, trueErrors, ref fixedCount, ref stillWrongCount, ref newlyWrongCount, ref totalChangedPositions);

                var totalInjectedErrors = trueErrors.Count;
                var recall = totalInjectedErrors == 0 ? 0.0 : (double)fixedCount / totalInjectedErrors;
                var falseCorrectionRate = totalChangedPositions == 0 ? 0.0 : (double)newlyWrongCount / totalChangedPositions;

                Console.WriteLine($"Injected errors: {totalInjectedErrors}");
                Console.WriteLine($"Fixed back to true base (recall): {fixedCount} ({recall:P2})");
                Console.WriteLine($"Still wrong (not fixed, or fixed to a different wrong base): {stillWrongCount}");
                Console.WriteLine($"Total positions changed by corrector: {totalChangedPositions}");
                Console.WriteLine($"Of those, changed a previously-CORRECT base to something wrong (false corrections): {newlyWrongCount} ({falseCorrectionRate:P2})");

                // 大まかな健全性チェック: recallは意味のある水準まで達し、
                // 誤訂正率は低く抑えられているべき。
                Assert.True(recall > 0.5, $"Expected recall > 50%, got {recall:P2}");
                Assert.True(falseCorrectionRate < 0.05, $"Expected false-correction rate < 5%, got {falseCorrectionRate:P2}");
            }
            finally
            {
                if (Directory.Exists(outDir))
                {
                    Directory.Delete(outDir, recursive: true);
                }
            }
        }

        private static Dictionary<(string, int), string> LoadReadsById(string path, int mate)
        {
            var result = new Dictionary<(string, int), string>();
            using var reader = new 簡易FASTQ読み込み(path);
            while (reader.Get_続きがあるか())
            {
                var (rawId, seq) = reader.Get_次のリード();
                var id = rawId.TrimStart('@').Split('/')[0];
                result[(id, mate)] = seq;
            }
            return result;
        }

        private static void ValidateFile(
            string correctedPath, int mate,
            Dictionary<(string, int), string> originalReads,
            Dictionary<(string ReadId, int Mate, int Position), char> trueErrors,
            ref int fixedCount, ref int stillWrongCount, ref int newlyWrongCount, ref int totalChangedPositions)
        {
            using var reader = new 簡易FASTQ読み込み(correctedPath);
            while (reader.Get_続きがあるか())
            {
                var (rawId, correctedSeq) = reader.Get_次のリード();
                var readId = rawId.TrimStart('@').Split('/')[0];
                if (!originalReads.TryGetValue((readId, mate), out var originalSeq))
                {
                    continue;
                }

                for (var pos = 0; pos < correctedSeq.Length && pos < originalSeq.Length; pos++)
                {
                    var wasError = trueErrors.TryGetValue((readId, mate, pos), out var trueBase);
                    var changed = correctedSeq[pos] != originalSeq[pos];

                    if (changed)
                    {
                        totalChangedPositions++;
                    }

                    if (wasError)
                    {
                        if (correctedSeq[pos] == trueBase)
                        {
                            fixedCount++;
                        }
                        else
                        {
                            stillWrongCount++;
                        }
                    }
                    else if (changed)
                    {
                        // 元々エラーではなかった(=正しかった)位置を書き換えてしまった。
                        newlyWrongCount++;
                    }
                }
            }
        }

        /// <summary>FASTQを「id, 配列」の2行単位として読むだけの軽量リーダー(品質行は無視)。</summary>
        private sealed class 簡易FASTQ読み込み(string p_パス) : IDisposable
        {
            private readonly StreamReader _reader = new(p_パス);

            public bool Get_続きがあるか() => !this._reader.EndOfStream;

            public (string A_ID, string A_配列) Get_次のリード()
            {
                var id = this._reader.ReadLine()!;
                var seq = this._reader.ReadLine()!;
                _ = this._reader.ReadLine(); // '+'
                _ = this._reader.ReadLine(); // quality
                return (id, seq);
            }

            public void Dispose() => this._reader.Dispose();
        }
    }
}
