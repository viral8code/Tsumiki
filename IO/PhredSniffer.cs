using Tsumiki.Model;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTQ のクオリティ文字列から Phred オフセット(33 or 64)を推定する。
    /// -p が明示指定されていない場合に限り推定値を自動採用し、明示指定されて
    /// いる場合は(ユーザーの判断を尊重して)警告のみに留める。
    /// </summary>
    internal static class PhredSniffer
    {
        // 実データで現実的にありうる最大のPhredスコア(Illumina/MGI/BGI いずれも
        // 通常は40強が上限)。これを大きく超えるスコアが観測された場合は、
        // オフセットの取り違えを疑う。
        private const int 現実的なQ上限 = 45;

        public static Phred標本 Get_標本(IEnumerable<string> p_クオリティ行, int p_標本上限 = 20_000)
        {
            var l_最小ASCII = int.MaxValue;
            var l_最大ASCII = int.MinValue;
            var l_リード数 = 0;
            var l_文字数 = 0;

            foreach (var l_クオリティ in p_クオリティ行)
            {
                if (l_リード数 >= p_標本上限)
                {
                    break;
                }
                l_リード数++;
                foreach (var l_文字 in l_クオリティ)
                {
                    l_文字数++;
                    if (l_文字 < l_最小ASCII)
                    {
                        l_最小ASCII = l_文字;
                    }
                    if (l_文字 > l_最大ASCII)
                    {
                        l_最大ASCII = l_文字;
                    }
                }
            }

            return l_文字数 == 0
                ? new Phred標本(0, 0, l_リード数, 0)
                : new Phred標本(l_最小ASCII, l_最大ASCII, l_リード数, l_文字数);
        }

        /// <summary>
        /// 標本が p_有効オフセット(現在有効な -p 値)と矛盾していそうな
        /// 場合に警告文を返す。問題なさそうな場合は null を返す。
        /// </summary>
        public static string? Get_警告文(Phred標本 p_標本, int p_有効オフセット)
        {
            if (p_標本.A_標本文字数 == 0)
            {
                return null;
            }

            var l_最小Q = p_標本.A_最小ASCII - p_有効オフセット;
            var l_最大Q = p_標本.A_最大ASCII - p_有効オフセット;
            var l_別のオフセット = p_有効オフセット == 33 ? 64 : 33;

            List<string> l_指摘 = [];
            if (l_最小Q < 0 || l_最大Q > 現実的なQ上限)
            {
                l_指摘.Add(
                    $"observed quality ASCII range [{p_標本.A_最小ASCII}, {p_標本.A_最大ASCII}] decodes to Q[{l_最小Q}, {l_最大Q}] " +
                    $"under Phred{p_有効オフセット}, which is implausible for real sequencing data " +
                    $"(negative or > {現実的なQ上限}). This data may actually be Phred{l_別のオフセット} " +
                    $"-- consider re-running with -p {l_別のオフセット} if so.");
            }
            if (p_標本.A_一様か)
            {
                l_指摘.Add(
                    $"quality is completely uniform (every sampled base is ASCII {p_標本.A_最小ASCII}) across " +
                    $"{p_標本.A_標本リード数} sampled read(s) -- this is unusual for real sequencer output and " +
                    "may indicate a placeholder/binned quality scheme rather than a genuine Phred offset mismatch.");
            }

            return l_指摘.Count == 0 ? null : string.Join(" ", l_指摘);
        }

        /// <summary>
        /// その標本を p_オフセット で解釈したとき、Q が負にならず、かつ
        /// 現実的な上限を超えないかどうか。
        /// </summary>
        private static bool Get_妥当なオフセットか(Phred標本 p_標本, int p_オフセット)
        {
            return p_標本.A_最小ASCII - p_オフセット >= 0 && p_標本.A_最大ASCII - p_オフセット <= 現実的なQ上限;
        }

        /// <summary>
        /// 標本から、どちらのオフセットが妥当かを判定する。
        /// 片方だけが妥当な場合にそのオフセットを返す。
        /// 両方妥当/両方不当な場合は判別できないため null を返す。
        /// </summary>
        public static int? Get_推定オフセット(Phred標本 p_標本)
        {
            if (p_標本.A_標本文字数 == 0)
            {
                return null;
            }

            var l_33が妥当 = Get_妥当なオフセットか(p_標本, 33);
            var l_64が妥当 = Get_妥当なオフセットか(p_標本, 64);
            if (l_33が妥当 == l_64が妥当)
            {
                return null;
            }
            return l_33が妥当 ? 33 : 64;
        }

        /// <summary>
        /// リードファイルを標本抽出して Phred オフセットを推定し、
        /// -p が明示指定されていなければ推定値を適用する。
        ///
        /// 自動適用する理由: クオリティによる k-mer 除外は
        /// 「クオリティ - Phredオフセット - クオリティカットオフ が負なら捨てる」で
        /// 判定するため、実際は Phred64 のデータを Phred33 として読むと、すべての
        /// 塩基のスコアが31以上に見えてしまい品質フィルタが事実上まったく効かなくなる
        /// (実データで ASCII 64 = Q0 の塩基がそのまま k-mer に使われていた)。
        /// 警告を出すだけでは静かに品質が落ちるため、判別がついた場合は
        /// 自動で正しい側に寄せる。
        ///
        /// read1 と read2 で推定結果が食い違う場合は自信が持てないため、
        /// 自動適用せず警告のみに留める。
        /// </summary>
        public static void V_解決_Phredオフセット(Parameters p_引数, string p_リード1のパス, string? p_リード2のパス, int p_標本上限 = 20_000)
        {
            var l_標本1 = Get_標本(Get_クオリティ行(p_リード1のパス, p_標本上限), p_標本上限);
            var l_推定 = Get_推定オフセット(l_標本1);

            if (!string.IsNullOrWhiteSpace(p_リード2のパス))
            {
                var l_標本2 = Get_標本(Get_クオリティ行(p_リード2のパス, p_標本上限), p_標本上限);
                var l_推定2 = Get_推定オフセット(l_標本2);
                if (l_推定 != l_推定2)
                {
                    Console.WriteLine(
                        "[Warning] Phred offset inference disagreed between the two read files " +
                        $"(read1 -> {l_推定?.ToString() ?? "undetermined"}, read2 -> {l_推定2?.ToString() ?? "undetermined"}). " +
                        $"Keeping -p {p_引数.A_Phredオフセット} as-is.");
                    l_推定 = null;
                }
            }

            if (l_推定 is { } l_オフセット && l_オフセット != p_引数.A_Phredオフセット)
            {
                if (p_引数.A_Phredが明示指定されたか)
                {
                    Console.WriteLine(
                        $"[Warning] Quality strings look like Phred{l_オフセット}, but -p {p_引数.A_Phredオフセット} was given explicitly. " +
                        $"Honouring the explicit value; re-run with -p {l_オフセット} if the data really is Phred{l_オフセット}.");
                }
                else
                {
                    p_引数.Set_推定Phredオフセット(l_オフセット);
                    Console.WriteLine(
                        $"[Info] Phred offset auto-detected as {l_オフセット} from the quality strings " +
                        $"(observed ASCII range [{l_標本1.A_最小ASCII}, {l_標本1.A_最大ASCII}]). Pass -p explicitly to override.");
                }
            }

            V_警告_疑わしいオフセット(p_リード1のパス, p_引数.A_Phredオフセット, p_標本上限);
            if (!string.IsNullOrWhiteSpace(p_リード2のパス))
            {
                V_警告_疑わしいオフセット(p_リード2のパス!, p_引数.A_Phredオフセット, p_標本上限);
            }
        }

        /// <summary>
        /// ファイルを標本抽出し、疑わしい場合はコンソールへ警告を出す。
        /// </summary>
        public static void V_警告_疑わしいオフセット(string p_ファイルパス, int p_有効オフセット, int p_標本上限 = 20_000)
        {
            var l_標本 = Get_標本(Get_クオリティ行(p_ファイルパス, p_標本上限), p_標本上限);
            var l_警告 = Get_警告文(l_標本, p_有効オフセット);
            if (l_警告 != null)
            {
                Console.WriteLine($"[Warning] Phred encoding check for {Path.GetFileName(p_ファイルパス)}: {l_警告}");
            }
        }

        private static IEnumerable<string> Get_クオリティ行(string p_ファイルパス, int p_標本上限)
        {
            using var l_読み込み = new FastqReader(p_ファイルパス);
            var l_件数 = 0;
            while (l_件数 < p_標本上限 && l_読み込み.Get_続きがあるか())
            {
                yield return l_読み込み.Get_次のリード().A_クオリティ;
                l_件数++;
            }
        }
    }
}
