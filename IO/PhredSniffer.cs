using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.IO
{
    /// <summary>
    /// FASTQ のクオリティ文字列から Phred オフセット (33 or 64) を推定する
    /// </summary>
    /// <remarks>
    /// -p が明示指定されていない場合に限り推定値を自動採用し、明示指定されている場合は (ユーザーの判断を尊重して) 警告のみに留める
    /// </remarks>
    internal static class PhredSniffer
    {
        #region 定数

        /// <summary>
        /// 現実的な Q 上限
        /// </summary>
        /// <remarks>
        /// 実データで現実的にありうる最大の Phred スコア<br/>
        /// (Illumina/MGI/BGI いずれも通常は 40 強が上限) <br/>
        /// これを大きく超えるスコアが観測された場合は、オフセットの取り違えを疑う
        /// </remarks>
        private const int 現実的なQ上限 = 45;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// クオリティ行から、オフセットを見分けるための標本を集めて返す
        /// </summary>
        /// <param name="p_クオリティ行">クオリティ行</param>
        /// <param name="p_標本上限">見る行数の上限</param>
        /// <returns>集めた標本</returns>
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
        /// 標本が p_有効オフセット (現在有効な -p 値) と矛盾していそうな場合に警告文を返す
        /// </summary>
        /// <param name="p_標本"></param>
        /// <param name="p_有効オフセット"></param>
        /// <remarks>
        /// 問題なさそうな場合は null を返す
        /// </remarks>
        /// <returns></returns>
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
                l_指摘.Add($"observed quality ASCII range [{p_標本.A_最小ASCII}, {p_標本.A_最大ASCII}] decodes to Q[{l_最小Q}, {l_最大Q}] under Phred{p_有効オフセット}, which is implausible for real sequencing data (negative or > {現実的なQ上限}). This data may actually be Phred{l_別のオフセット} -- consider re-running with -p {l_別のオフセット} if so.");
            }

            if (p_標本.A_一様か)
            {
                l_指摘.Add($"quality is completely uniform (every sampled base is ASCII {p_標本.A_最小ASCII}) across {p_標本.A_標本リード数} sampled read(s) -- this is unusual for real sequencer output and may indicate a placeholder/binned quality scheme rather than a genuine Phred offset mismatch.");
            }

            return l_指摘.Count == 0 ? null : string.Join(" ", l_指摘);
        }

        /// <summary>
        /// 標本から、どちらのオフセットが妥当かを判定する
        /// </summary>
        /// <param name="p_標本"></param>
        /// <remarks>
        /// 片方だけが妥当な場合にそのオフセットを返す<br/>
        /// 両方妥当/両方不当な場合は判別できないため null を返す
        /// </remarks>
        /// <returns></returns>
        public static int? Get_推定オフセット(Phred標本 p_標本)
        {
            if (p_標本.A_標本文字数 == 0)
            {
                return null;
            }

            var l_33が妥当 = Get_妥当なオフセットか(p_標本, 33);
            var l_64が妥当 = Get_妥当なオフセットか(p_標本, 64);
            return l_33が妥当 == l_64が妥当 ? null : l_33が妥当 ? 33 : 64;
        }

        /// <summary>
        /// クオリティ文字から Phred オフセットを推定し、実行時引数へ反映する
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_リード1のパス">リード 1 のパス</param>
        /// <param name="p_リード2のパス">リード 2 のパス、単一リードなら null</param>
        /// <param name="p_標本上限">見る行数の上限</param>
        /// <remarks>
        /// -p で明示指定されている場合は、利用者の判断を優先して推定結果で上書きしない<br/>
        /// 警告だけでは足りない<br/>
        /// Phred64 のデータを Phred33 として読むとすべてのスコアが 31 以上に見え、品質フィルタが完全に無効化されるが、その事実はログを読まない限り気付けない<br/>
        /// read1 と read2 で推定が食い違う場合は自信が持てないため警告に留める
        /// </remarks>
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
                    Logger.V_出力(メッセージID.Phred_ファイル間で不一致, Get_表示用オフセット(l_推定), Get_表示用オフセット(l_推定2), p_引数.A_Phredオフセット);
                    l_推定 = null;
                }
            }

            if (l_推定 is { } l_オフセット && l_オフセット != p_引数.A_Phredオフセット)
            {
                if (p_引数.A_Phredが明示指定されたか)
                {
                    Logger.V_出力(メッセージID.Phred_明示指定と不一致, l_オフセット, p_引数.A_Phredオフセット, l_オフセット, l_オフセット);
                }
                else
                {
                    p_引数.Set_推定Phredオフセット(l_オフセット);
                    Logger.V_出力(メッセージID.Phred_自動判定, l_オフセット, l_標本1.A_最小ASCII, l_標本1.A_最大ASCII);
                }
            }

            V_警告_疑わしいオフセット(p_リード1のパス, p_引数.A_Phredオフセット, p_標本上限);
            if (!string.IsNullOrWhiteSpace(p_リード2のパス))
            {
                V_警告_疑わしいオフセット(p_リード2のパス!, p_引数.A_Phredオフセット, p_標本上限);
            }
        }

        /// <summary>
        /// ファイルを標本抽出し、疑わしい場合はコンソールへ警告を出す
        /// </summary>
        /// <param name="p_ファイルパス"></param>
        /// <param name="p_有効オフセット"></param>
        /// <param name="p_標本上限"></param>
        public static void V_警告_疑わしいオフセット(string p_ファイルパス, int p_有効オフセット, int p_標本上限 = 20_000)
        {
            var l_標本 = Get_標本(Get_クオリティ行(p_ファイルパス, p_標本上限), p_標本上限);
            var l_警告 = Get_警告文(l_標本, p_有効オフセット);
            if (l_警告 != null)
            {
                Logger.V_出力(メッセージID.Phred_検査の警告, Path.GetFileName(p_ファイルパス), l_警告);
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// その標本を p_オフセット で解釈したとき、Q が負にならず、かつ現実的な上限を超えないかどうか
        /// </summary>
        /// <param name="p_標本"></param>
        /// <param name="p_オフセット"></param>
        /// <returns></returns>
        private static bool Get_妥当なオフセットか(Phred標本 p_標本, int p_オフセット)
        {
            return p_標本.A_最小ASCII - p_オフセット >= 0 && p_標本.A_最大ASCII - p_オフセット <= 現実的なQ上限;
        }

        /// <summary>
        /// 推定できなかった場合も含めた、表示用のオフセット
        /// </summary>
        /// <param name="p_推定"></param>
        /// <returns></returns>
        private static string Get_表示用オフセット(int? p_推定)
        {
            return p_推定?.ToString() ?? Messages.Get_文言(メッセージID.Phred_未確定);
        }

        /// <summary>
        /// FASTQ からクオリティ行だけを取り出して返す
        /// </summary>
        /// <param name="p_ファイルパス">FASTQ のパス</param>
        /// <param name="p_標本上限">取り出す行数の上限</param>
        /// <returns>クオリティ行</returns>
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

        #endregion
    }
}
