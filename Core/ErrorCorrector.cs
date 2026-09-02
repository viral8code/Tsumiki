using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// k-merスペクトラムに基づく、Quake/BayesHammer類似の簡易リードエラー訂正。
    ///
    /// 考え方: あるk-merが十分な回数(カットオフ以上)出現していれば「信頼できる」
    /// (真のゲノム配列由来である可能性が高い)とみなし、逆に出現回数が少ない
    /// k-merは配列決定エラーに由来する可能性が高いとみなす。1本のリード中に
    /// 信頼できないk-merが含まれる場合、その周辺の1塩基を置換することで
    /// 信頼できるk-merに変わるかどうかを試し、最も多くの「信頼できないk-mer窓」を
    /// 信頼できる状態に変える置換を貪欲に選んで適用する(1リードにつき複数回、
    /// 改善が見込めなくなるまで反復する)。
    ///
    /// 曖昧塩基(N等、Consts.InvalidBase)を含む位置は書き換えない
    /// (置換候補にも含めない)。そのような位置を含む窓は常に「信頼できない」
    /// 扱いとし、スコアリングの対象からも除外する。
    /// </summary>
    internal static class ErrorCorrector
    {
        public readonly record struct CorrectionResult(byte[] Read, int CorrectionCount);

        /// <summary>
        /// readPath1(+readPath2)を読み込んでエラー訂正を行い、
        /// 結果を outPath1(+outPath2)に書き出す。
        /// 「信頼できるk-mer」の判定には、本アセンブリと同じ -kc のカットオフ値を
        /// 使って構築した専用の TrustedKmerIndex(このメソッド内で完結し、
        /// 本パイプライン用のインデックスとは独立)を用いる。
        /// </summary>
        public static void CorrectReadFiles(string readPath1, string? readPath2, string tempDir, string outPath1, string? outPath2)
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            var cutoff = ConfigurationManager.Arguments.KmerCutoff;

            var ecTempDir = Path.Combine(tempDir, "error_correction");
            _ = Directory.CreateDirectory(ecTempDir);

            Console.WriteLine("[ErrorCorrection] Building k-mer spectrum...");
            using (var trustedIndex = new TrustedKmerIndex(ecTempDir))
            {
                KmerCounting.LoadReadFile(readPath1, trustedIndex);
                if (readPath2 != null)
                {
                    KmerCounting.LoadReadFile(readPath2, trustedIndex);
                }
                _ = trustedIndex.Cutoff(cutoff);

                Console.WriteLine("[ErrorCorrection] Correcting reads...");
                var stats1 = CorrectFile(readPath1, outPath1, trustedIndex, kmerLength);
                Console.WriteLine($"[ErrorCorrection] {Path.GetFileName(readPath1)}: " +
                    $"{stats1.CorrectedReads}/{stats1.TotalReads} reads corrected ({stats1.TotalCorrections} base corrections total).");

                if (readPath2 != null && outPath2 != null)
                {
                    var stats2 = CorrectFile(readPath2, outPath2, trustedIndex, kmerLength);
                    Console.WriteLine($"[ErrorCorrection] {Path.GetFileName(readPath2)}: " +
                        $"{stats2.CorrectedReads}/{stats2.TotalReads} reads corrected ({stats2.TotalCorrections} base corrections total).");
                }
            }

            Directory.Delete(ecTempDir, recursive: true);
        }

        private readonly record struct FileCorrectionStats(int TotalReads, int CorrectedReads, int TotalCorrections);

        private static FileCorrectionStats CorrectFile(string inPath, string outPath, TrustedKmerIndex trustedIndex, int kmerLength)
        {
            var totalReads = 0;
            var correctedReads = 0;
            var totalCorrections = 0;

            using var reader = new FastqReader(inPath);
            using var writer = new FastqWriter(outPath);
            while (reader.HasNext())
            {
                var readData = reader.NextReadSimple();
                totalReads++;
                var result = CorrectRead(readData.SimpleRead!, trustedIndex, kmerLength);
                if (result.CorrectionCount > 0)
                {
                    correctedReads++;
                    totalCorrections += result.CorrectionCount;
                }
                var correctedSeq = string.Join(string.Empty, result.Read.Select(Util.ByteToBaseString));
                writer.Write(readData.ID, correctedSeq, readData.Quality);
            }

            return new FileCorrectionStats(totalReads, correctedReads, totalCorrections);
        }

        /// <summary>
        /// 1リード(Consts.NucleotideID空間のバイト列、曖昧塩基はInvalidBase)を
        /// 貪欲法で訂正する。副作用のない純粋関数(入力readは変更しない)。
        /// </summary>
        public static CorrectionResult CorrectRead(ReadOnlySpan<byte> read, TrustedKmerIndex trustedIndex, int kmerLength, int maxIterations = 10)
        {
            if (read.Length < kmerLength)
            {
                return new CorrectionResult(read.ToArray(), 0);
            }

            var bases = read.ToArray();
            var numWindows = bases.Length - kmerLength + 1;
            var corrections = 0;

            for (var iteration = 0; iteration < maxIterations; iteration++)
            {
                var trusted = ComputeWindowTrust(bases, kmerLength, trustedIndex);
                if (Array.TrueForAll(trusted, t => t))
                {
                    break;
                }

                var bestPosition = -1;
                byte bestBase = 0;
                var bestGain = 0;

                for (var pos = 0; pos < bases.Length; pos++)
                {
                    if (bases[pos] == Consts.InvalidBase)
                    {
                        continue;
                    }

                    var windowStart = Math.Max(0, pos - kmerLength + 1);
                    var windowEnd = Math.Min(numWindows - 1, pos);

                    var anyUntrusted = false;
                    for (var w = windowStart; w <= windowEnd; w++)
                    {
                        if (!trusted[w])
                        {
                            anyUntrusted = true;
                            break;
                        }
                    }
                    if (!anyUntrusted)
                    {
                        continue;
                    }

                    var currentBase = bases[pos];
                    for (byte candidate = 1; candidate <= 4; candidate++)
                    {
                        if (candidate == currentBase)
                        {
                            continue;
                        }

                        var gain = EvaluateSubstitution(bases, pos, candidate, windowStart, windowEnd, kmerLength, trusted, trustedIndex);
                        if (gain > bestGain)
                        {
                            bestGain = gain;
                            bestPosition = pos;
                            bestBase = candidate;
                        }
                    }
                }

                if (bestPosition < 0)
                {
                    // これ以上、信頼できる窓を純増させる置換が見つからない
                    // (=残った信頼できない窓は、単発の置換では解決できない)。
                    break;
                }

                bases[bestPosition] = bestBase;
                corrections++;
            }

            return new CorrectionResult(bases, corrections);
        }

        private static bool[] ComputeWindowTrust(byte[] bases, int kmerLength, TrustedKmerIndex trustedIndex)
        {
            var numWindows = bases.Length - kmerLength + 1;
            var trusted = new bool[numWindows];
            for (var w = 0; w < numWindows; w++)
            {
                trusted[w] = IsWindowTrusted(bases, w, kmerLength, trustedIndex);
            }
            return trusted;
        }

        private static bool IsWindowTrusted(byte[] bases, int windowStart, int kmerLength, TrustedKmerIndex trustedIndex)
        {
            for (var i = windowStart; i < windowStart + kmerLength; i++)
            {
                if (bases[i] == Consts.InvalidBase)
                {
                    return false;
                }
            }
            return trustedIndex.Contains(bases.AsSpan(windowStart, kmerLength));
        }

        /// <summary>
        /// position を candidate に置換した場合の「信頼できる窓の純増数」を計算する
        /// (windowStart..windowEnd の範囲、すなわち position を含みうる窓のみが
        /// 影響を受けるため、その範囲だけを再評価すれば十分)。bases は評価後、
        /// 呼び出し前の状態に戻す(副作用を残さない)。
        /// </summary>
        private static int EvaluateSubstitution(
            byte[] bases, int position, byte candidate, int windowStart, int windowEnd, int kmerLength,
            bool[] trustedBefore, TrustedKmerIndex trustedIndex)
        {
            var original = bases[position];
            bases[position] = candidate;

            var gain = 0;
            for (var w = windowStart; w <= windowEnd; w++)
            {
                var nowTrusted = IsWindowTrusted(bases, w, kmerLength, trustedIndex);
                if (nowTrusted && !trustedBefore[w])
                {
                    gain++;
                }
                else if (!nowTrusted && trustedBefore[w])
                {
                    gain--;
                }
            }

            bases[position] = original;
            return gain;
        }
    }
}
