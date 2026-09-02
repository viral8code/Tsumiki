namespace Tsumiki.IO
{
    /// <summary>
    /// FASTQ のクオリティ文字列から、現在有効な Phred オフセット(33 or 64)が
    /// もっともらしいかどうかを推測する。オフセットの自動切り替えは行わず、
    /// 怪しい場合に警告を出すだけに留める(判断はユーザーに委ねる)。
    /// </summary>
    internal readonly record struct PhredSample(int MinAscii, int MaxAscii, int SampledReads, int SampledChars)
    {
        /// <summary>
        /// サンプル全体を通してASCIIコードが一切変化しなかったか
        /// (実機のシーケンサ出力では通常あり得ない、人工的/ビニング済みの
        /// クオリティである可能性を示す)。
        /// </summary>
        public bool IsUniform => this.SampledChars > 0 && this.MinAscii == this.MaxAscii;
    }

    internal static class PhredSniffer
    {
        // 実データで現実的にありうる最大のPhredスコア(Illumina/MGI/BGI いずれも
        // 通常は40強が上限)。これを大きく超えるスコアが観測された場合は、
        // オフセットの取り違えを疑う。
        private const int PlausibleMaxQ = 45;

        public static PhredSample Sample(IEnumerable<string> qualityLines, int maxReadsToSample = 20_000)
        {
            var minAscii = int.MaxValue;
            var maxAscii = int.MinValue;
            var readCount = 0;
            var charCount = 0;

            foreach (var quality in qualityLines)
            {
                if (readCount >= maxReadsToSample)
                {
                    break;
                }
                readCount++;
                foreach (var c in quality)
                {
                    charCount++;
                    if (c < minAscii)
                    {
                        minAscii = c;
                    }
                    if (c > maxAscii)
                    {
                        maxAscii = c;
                    }
                }
            }

            return charCount == 0
                ? new PhredSample(0, 0, readCount, 0)
                : new PhredSample(minAscii, maxAscii, readCount, charCount);
        }

        /// <summary>
        /// サンプル結果が phredOffsetInEffect(現在有効な -p 値)と矛盾していそうな
        /// 場合に警告文を返す。問題なさそうな場合は null を返す。
        /// </summary>
        public static string? BuildWarning(PhredSample sample, int phredOffsetInEffect)
        {
            if (sample.SampledChars == 0)
            {
                return null;
            }

            var minQ = sample.MinAscii - phredOffsetInEffect;
            var maxQ = sample.MaxAscii - phredOffsetInEffect;
            var otherOffset = phredOffsetInEffect == 33 ? 64 : 33;

            List<string> issues = [];
            if (minQ < 0 || maxQ > PlausibleMaxQ)
            {
                issues.Add(
                    $"observed quality ASCII range [{sample.MinAscii}, {sample.MaxAscii}] decodes to Q[{minQ}, {maxQ}] " +
                    $"under Phred{phredOffsetInEffect}, which is implausible for real sequencing data " +
                    $"(negative or > {PlausibleMaxQ}). This data may actually be Phred{otherOffset} " +
                    $"-- consider re-running with -p {otherOffset} if so.");
            }
            if (sample.IsUniform)
            {
                issues.Add(
                    $"quality is completely uniform (every sampled base is ASCII {sample.MinAscii}) across " +
                    $"{sample.SampledReads} sampled read(s) -- this is unusual for real sequencer output and " +
                    "may indicate a placeholder/binned quality scheme rather than a genuine Phred offset mismatch.");
            }

            return issues.Count == 0 ? null : string.Join(" ", issues);
        }

        /// <summary>
        /// filePath をサンプリングし、疑わしい場合はコンソールへ警告を出す。
        /// </summary>
        public static void WarnIfImplausible(string filePath, int phredOffsetInEffect, int maxReadsToSample = 20_000)
        {
            var sample = Sample(SampleQualityLines(filePath, maxReadsToSample), maxReadsToSample);
            var warning = BuildWarning(sample, phredOffsetInEffect);
            if (warning != null)
            {
                Console.WriteLine($"[Warning] Phred encoding check for {Path.GetFileName(filePath)}: {warning}");
            }
        }

        private static IEnumerable<string> SampleQualityLines(string filePath, int maxReadsToSample)
        {
            using var reader = new FastqReader(filePath);
            var count = 0;
            while (count < maxReadsToSample && reader.HasNext())
            {
                yield return reader.NextRead().Quality;
                count++;
            }
        }
    }
}
