using Tsumiki.Model;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTQ のクオリティ文字列から Phred オフセット(33 or 64)を推定する。
    /// -p が明示指定されていない場合に限り推定値を自動採用し、明示指定されて
    /// いる場合は(ユーザーの判断を尊重して)警告のみに留める。
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
        /// サンプルから、どちらのオフセットが妥当かを判定する。
        /// 片方だけが「Q が負にならず、かつ現実的な上限を超えない」を満たす場合に
        /// そのオフセットを返す。両方満たす/両方満たさない場合は判別できないため
        /// null を返す。
        /// </summary>
        public static int? InferOffset(PhredSample sample)
        {
            if (sample.SampledChars == 0)
            {
                return null;
            }

            bool IsPlausible(int offset)
            {
                return sample.MinAscii - offset >= 0 && sample.MaxAscii - offset <= PlausibleMaxQ;
            }

            var plausible33 = IsPlausible(33);
            var plausible64 = IsPlausible(64);
            if (plausible33 == plausible64)
            {
                return null;
            }
            return plausible33 ? 33 : 64;
        }

        /// <summary>
        /// リードファイルをサンプリングして Phred オフセットを推定し、
        /// -p が明示指定されていなければ推定値を param に適用する。
        ///
        /// 自動適用する理由: クオリティによる k-mer 除外は
        /// 「quality - Phred - QualityCutoff が負なら捨てる」で判定するため、
        /// 実際は Phred64 のデータを Phred33 として読むと、すべての塩基のスコアが
        /// 31 以上に見えてしまい品質フィルタが事実上まったく効かなくなる
        /// (実データで ASCII 64 = Q0 の塩基がそのまま k-mer に使われていた)。
        /// 警告を出すだけでは静かに品質が落ちるため、判別がついた場合は
        /// 自動で正しい側に寄せる。
        ///
        /// read1 と read2 で推定結果が食い違う場合は自信が持てないため、
        /// 自動適用せず警告のみに留める。
        /// </summary>
        public static void ResolveOffset(Parameters param, string readPath1, string? readPath2, int maxReadsToSample = 20_000)
        {
            var sample1 = Sample(SampleQualityLines(readPath1, maxReadsToSample), maxReadsToSample);
            var inferred = InferOffset(sample1);

            if (!string.IsNullOrWhiteSpace(readPath2))
            {
                var sample2 = Sample(SampleQualityLines(readPath2, maxReadsToSample), maxReadsToSample);
                var inferred2 = InferOffset(sample2);
                if (inferred != inferred2)
                {
                    Console.WriteLine(
                        "[Warning] Phred offset inference disagreed between the two read files " +
                        $"(read1 -> {inferred?.ToString() ?? "undetermined"}, read2 -> {inferred2?.ToString() ?? "undetermined"}). " +
                        $"Keeping -p {param.Phred} as-is.");
                    inferred = null;
                }
            }

            if (inferred is { } offset && offset != param.Phred)
            {
                if (param.IsPhredExplicitlySet)
                {
                    Console.WriteLine(
                        $"[Warning] Quality strings look like Phred{offset}, but -p {param.Phred} was given explicitly. " +
                        $"Honouring the explicit value; re-run with -p {offset} if the data really is Phred{offset}.");
                }
                else
                {
                    param.SetInferredPhred(offset);
                    Console.WriteLine(
                        $"[Info] Phred offset auto-detected as {offset} from the quality strings " +
                        $"(observed ASCII range [{sample1.MinAscii}, {sample1.MaxAscii}]). Pass -p explicitly to override.");
                }
            }

            WarnIfImplausible(readPath1, param.Phred, maxReadsToSample);
            if (!string.IsNullOrWhiteSpace(readPath2))
            {
                WarnIfImplausible(readPath2!, param.Phred, maxReadsToSample);
            }
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
