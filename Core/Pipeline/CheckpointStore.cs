using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using Tsumiki.Common;
using Tsumiki.Model;

namespace Tsumiki.Core
{
    /// <summary>
    /// k ごとの成果を作業ディレクトリに残し、次回の実行でそこから再開できるようにする。
    ///
    /// multi-k は k の本数だけ全工程を繰り返すため、後半の k で落ちたときに
    /// 前半をやり直すのは丸ごと無駄になる。既に同じ条件で作り終えている k は
    /// 読み直して飛ばす。
    ///
    /// 条件が少しでも違えば再利用してはならない。入力ファイル・k・
    /// 結果に効くオプションから署名を作り、一致したときだけ再開する。
    /// 更新時刻ではなくファイル長を見るのは、前処理やエラー訂正を挟むと
    /// 中間ファイルが作り直されて時刻だけが変わるため。
    /// </summary>
    internal static class CheckpointStore
    {
        private const string チェックポイントファイル名 = "checkpoint.tsv";

        private const string 引き継ぎファイル名 = "carryover.tsv";

        /// <summary>成果物が無いことを示す印。</summary>
        private const string 無し = "-";

        /// <summary>
        /// この k の結果を再利用してよいかを決める署名。
        /// 結果を変えうる入力とオプションだけを混ぜる。
        /// </summary>
        public static string Get_署名(Parameters p_引数, int p_k長)
        {
            var l_材料 = string.Join(
                '\n',
                Consts.バージョン,
                Get_ファイル指紋(p_引数.A_リード1のパス),
                Get_ファイル指紋(p_引数.A_リード2のパス),
                $"k={p_k長}",
                $"kc={(p_引数.A_kmerカットオフが明示指定されたか ? p_引数.A_kmerカットオフ.ToString() : "auto")}",
                $"q={p_引数.A_クオリティカットオフ}",
                $"phred={p_引数.A_Phredオフセット}",
                $"pu={p_引数.A_ペア結合閾値.ToString(CultureInfo.InvariantCulture)}",
                $"pc={p_引数.A_ペア支持数閾値}",
                $"i={p_引数.A_インサートサイズ?.ToString() ?? 無し}",
                $"ab={p_引数.A_曖昧塩基を許容するか}",
                $"nc={!p_引数.A_引き継ぐか}",
                $"sr={p_引数.A_SuperReadを作るか}",
                $"rv={p_引数.A_反復をrMerで検証するか}",
                $"la={p_引数.A_局所アセンブリするか}",
                $"my={p_引数.A_救済kmerを使うか}",
                $"gfa={p_引数.A_GFAを出力するか}");
            return Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(l_材料)));
        }

        /// <summary>
        /// 作業ディレクトリに、同じ署名で作り終えた結果があれば読み直す。
        /// 無い・署名が違う・成果物が欠けている場合は null。
        /// </summary>
        public static アセンブリ実行結果? Get_再開結果(
            string p_作業ディレクトリ, string p_署名, List<引き継ぎ配列>? p_次への引き継ぎ)
        {
            var l_パス = Path.Combine(p_作業ディレクトリ, チェックポイントファイル名);
            if (!File.Exists(l_パス))
            {
                return null;
            }

            Dictionary<string, string> l_項目 = [];
            foreach (var l_行 in File.ReadLines(l_パス))
            {
                var l_区切り = l_行.IndexOf('\t');
                if (l_区切り > 0)
                {
                    l_項目[l_行[..l_区切り]] = l_行[(l_区切り + 1)..];
                }
            }

            if (l_項目.GetValueOrDefault("signature") != p_署名)
            {
                return null;
            }

            var l_ユニティグ = Get_パス(l_項目, "unitigs");
            var l_コンティグ = Get_パス(l_項目, "contigs");
            if (l_ユニティグ is null || l_コンティグ is null)
            {
                return null;
            }
            var l_スキャフォールド = Get_パス(l_項目, "scaffolds");
            if (l_項目.GetValueOrDefault("scaffolds", 無し) != 無し && l_スキャフォールド is null)
            {
                // 作ったはずのものが消えている。信用できないので作り直す。
                return null;
            }

            if (!int.TryParse(l_項目.GetValueOrDefault("k"), out var l_k長)
                || !ulong.TryParse(l_項目.GetValueOrDefault("cutoff"), out var l_カットオフ)
                || !double.TryParse(
                    l_項目.GetValueOrDefault("single_copy"),
                    NumberStyles.Float, CultureInfo.InvariantCulture, out var l_基準値))
            {
                return null;
            }

            if (p_次への引き継ぎ is not null)
            {
                p_次への引き継ぎ.Clear();
                p_次への引き継ぎ.AddRange(Get_引き継ぎ(p_作業ディレクトリ));
            }

            return new アセンブリ実行結果(
                l_k長, l_ユニティグ, l_コンティグ, l_スキャフォールド,
                l_カットオフ, l_基準値, Get_パス(l_項目, "gfa"),
                Get_整合性検査(l_項目));
        }

        /// <summary>この k を作り終えたことを記録する。</summary>
        public static void V_保存(
            string p_作業ディレクトリ, string p_署名, アセンブリ実行結果 p_結果,
            IReadOnlyList<引き継ぎ配列>? p_次への引き継ぎ)
        {
            var l_文 = new StringBuilder();
            V_追加(l_文, "signature", p_署名);
            V_追加(l_文, "k", p_結果.A_k長.ToString());
            V_追加(l_文, "cutoff", p_結果.A_kmerカットオフ.ToString());
            V_追加(
                l_文, "single_copy",
                p_結果.A_単一コピー基準値.ToString("R", CultureInfo.InvariantCulture));
            V_追加(l_文, "unitigs", p_結果.A_ユニティグパス);
            V_追加(l_文, "contigs", p_結果.A_コンティグパス);
            V_追加(l_文, "scaffolds", p_結果.A_スキャフォールドパス ?? 無し);
            V_追加(l_文, "gfa", p_結果.A_GFAパス ?? 無し);
            if (p_結果.A_整合性検査 is { } l_検査)
            {
                V_追加(l_文, "trusted_kmers", l_検査.A_信頼kmer数.ToString());
                V_追加(l_文, "assembly_kmers", l_検査.A_アセンブリ内の延べ数.ToString());
                V_追加(l_文, "assembly_distinct", l_検査.A_アセンブリ内の種類数.ToString());
                V_追加(l_文, "missing_kmers", l_検査.A_取りこぼし数.ToString());
                V_追加(l_文, "excess_kinds", l_検査.A_出しすぎkmer種類数.ToString());
                V_追加(l_文, "excess_kmers", l_検査.A_余分な延べ数.ToString());
            }

            V_保存_引き継ぎ(p_作業ディレクトリ, p_次への引き継ぎ);

            // 引き継ぎを書き終えてから署名を置く。途中で落ちた場合に
            // 揃っていない状態を再開に使わせない。
            File.WriteAllText(Path.Combine(p_作業ディレクトリ, チェックポイントファイル名), l_文.ToString());
        }

        private static void V_保存_引き継ぎ(
            string p_作業ディレクトリ, IReadOnlyList<引き継ぎ配列>? p_次への引き継ぎ)
        {
            var l_パス = Path.Combine(p_作業ディレクトリ, 引き継ぎファイル名);
            if (p_次への引き継ぎ is null)
            {
                if (File.Exists(l_パス))
                {
                    File.Delete(l_パス);
                }
                return;
            }

            using var l_書き込み = new StreamWriter(l_パス);
            foreach (var l_配列 in p_次への引き継ぎ)
            {
                l_書き込み.Write(l_配列.A_配列);
                l_書き込み.Write('\t');
                l_書き込み.Write(l_配列.A_k長);
                l_書き込み.Write('\t');
                l_書き込み.WriteLine(string.Join(',', l_配列.A_カバレッジ));
            }
        }

        private static IEnumerable<引き継ぎ配列> Get_引き継ぎ(string p_作業ディレクトリ)
        {
            var l_パス = Path.Combine(p_作業ディレクトリ, 引き継ぎファイル名);
            if (!File.Exists(l_パス))
            {
                yield break;
            }
            foreach (var l_行 in File.ReadLines(l_パス))
            {
                var l_列 = l_行.Split('\t');
                if (l_列.Length != 3 || !int.TryParse(l_列[1], out var l_k長))
                {
                    continue;
                }
                var l_カバレッジ = l_列[2].Length == 0
                    ? []
                    : l_列[2].Split(',').Select(int.Parse).ToArray();
                yield return new 引き継ぎ配列(l_列[0], l_カバレッジ, l_k長);
            }
        }

        private static 整合性検査結果? Get_整合性検査(IReadOnlyDictionary<string, string> p_項目)
        {
            return long.TryParse(p_項目.GetValueOrDefault("trusted_kmers"), out var l_信頼)
                && long.TryParse(p_項目.GetValueOrDefault("assembly_kmers"), out var l_延べ)
                && long.TryParse(p_項目.GetValueOrDefault("assembly_distinct"), out var l_種類)
                && long.TryParse(p_項目.GetValueOrDefault("missing_kmers"), out var l_取りこぼし)
                && long.TryParse(p_項目.GetValueOrDefault("excess_kinds"), out var l_出しすぎ種類)
                && long.TryParse(p_項目.GetValueOrDefault("excess_kmers"), out var l_余分)
                ? new 整合性検査結果(l_信頼, l_延べ, l_種類, l_取りこぼし, l_出しすぎ種類, l_余分)
                : null;
        }

        private static string? Get_パス(IReadOnlyDictionary<string, string> p_項目, string p_キー)
        {
            var l_値 = p_項目.GetValueOrDefault(p_キー, 無し);
            return l_値 != 無し && File.Exists(l_値) ? l_値 : null;
        }

        private static void V_追加(StringBuilder p_文, string p_キー, string p_値)
        {
            _ = p_文.Append(p_キー).Append('\t').AppendLine(p_値);
        }

        /// <summary>
        /// 入力ファイルの同一性。中身のハッシュは規模的に取れないため、
        /// パスと長さで代用する。
        /// </summary>
        private static string Get_ファイル指紋(string? p_パス)
        {
            if (string.IsNullOrWhiteSpace(p_パス) || !File.Exists(p_パス))
            {
                return 無し;
            }
            var l_情報 = new FileInfo(p_パス);
            return $"{l_情報.Name}:{l_情報.Length}";
        }
    }
}
