using System.Text;
using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Output;
using Tsumiki.Cores.Polishing;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Polishing;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// 最終成果物の仕上げと検査を管理する
    /// </summary>
    internal static class FinalAssemblyPipeline
    {
        #region 定数

        /// <summary>
        /// 見出し agp version 2 1
        /// </summary>
        private const string C_見出し_agp_version_2_1 = "##agp-version	2.1";

        /// <summary>
        /// 項目 U
        /// </summary>
        private const string C_項目_U = "U";

        /// <summary>
        /// 長さ指定ギャップ種別
        /// </summary>
        private const string C_長さ指定ギャップ種別 = "N";

        /// <summary>
        /// 余ったリードで端を繋ぐときの、跨ぐリードの必要本数<br/>
        /// 余ったリードを含む全リードは間引いたリードより深く、読み誤りや異常なリードも多く跨ぐため、要るリード数を 1 つ上げる
        /// </summary>
        private const int C_余りで要るリード数 = 3;

        /// <summary>
        /// 端として見る窓の長さ（先頭と末尾それぞれ）
        /// </summary>
        private const int C_端の窓長 = 300;

        /// <summary>
        /// 端に触れるかを判定するキーの長さ
        /// </summary>
        private const int C_端の種長 = 25;

        /// <summary>
        /// 端を繋ぐ対象にする配列の最短長
        /// </summary>
        private const int C_端つなぎの最短長 = 150;

        /// <summary>
        /// 項目 contigs
        /// </summary>
        private const string C_項目_contigs = "contigs";

        /// <summary>
        /// 項目 scaffolds
        /// </summary>
        private const string C_項目_scaffolds = "scaffolds";

        /// <summary>
        /// 項目 assembly
        /// </summary>
        private const string C_項目_assembly = "assembly";

        /// <summary>
        /// ファイル名 assembly provenance json
        /// </summary>
        private const string C_ファイル名_assembly_provenance_json = "assembly.provenance.json";

        /// <summary>
        /// 項目 schema version
        /// </summary>
        private const string C_項目_schema_version = "schema_version";

        /// <summary>
        /// 項目 build id
        /// </summary>
        private const string C_項目_build_id = "build_id";

        /// <summary>
        /// 項目 validation source
        /// </summary>
        private const string C_項目_validation_source = "validation_source";

        /// <summary>
        /// 項目 uncorrected reads from t ent holdout data
        /// </summary>
        private const string C_項目_uncorrected_reads_from_t_ent_holdout_data = "uncorrected reads from the assembly library; not independent holdout data";

        /// <summary>
        /// 項目 assembly sha256
        /// </summary>
        private const string C_項目_assembly_sha256 = "assembly_sha256";

        /// <summary>
        /// 項目 thread count
        /// </summary>
        private const string C_項目_thread_count = "thread_count";

        /// <summary>
        /// 項目 settings
        /// </summary>
        private const string C_項目_settings = "settings";

        /// <summary>
        /// 項目 assembly settings
        /// </summary>
        private const string C_項目_assembly_settings = "assembly_settings";

        /// <summary>
        /// 項目 copy number baseline requested
        /// </summary>
        private const string C_項目_copy_number_baseline_requested = "copy_number_baseline_requested";

        /// <summary>
        /// 項目 copy number baseline actual
        /// </summary>
        private const string C_項目_copy_number_baseline_actual = "copy_number_baseline_actual";

        /// <summary>
        /// 項目 trim low coverage ends
        /// </summary>
        private const string C_項目_trim_low_coverage_ends = "trim_low_coverage_ends";

        /// <summary>
        /// 項目 inputs
        /// </summary>
        private const string C_項目_inputs = "inputs";

        /// <summary>
        /// 項目 path
        /// </summary>
        private const string C_項目_path = "path";

        /// <summary>
        /// 項目 sha256
        /// </summary>
        private const string C_項目_sha256 = "sha256";

        /// <summary>
        /// 項目 corrected read hashes
        /// </summary>
        private const string C_項目_corrected_read_hashes = "corrected_read_hashes";

        /// <summary>
        /// 項目 pipeline fingerprint
        /// </summary>
        private const string C_項目_pipeline_fingerprint = "pipeline_fingerprint";

        /// <summary>
        /// 支持のない箇所のファイル名
        /// </summary>
        private const string C_支持のない箇所ファイル名 = "assembly.unsupported.tsv";

        /// <summary>
        /// リード支持検査の r 長
        /// </summary>
        private const int C_支持検査のr長 = 31;

        /// <summary>
        /// ポリッシュ結果の一時ファイル名
        /// </summary>
        private const string C_ポリッシュ済みファイル名 = "polished.fasta";

        /// <summary>
        /// 長さが分からない繋ぎ目に入れる N の数 (NCBI の慣例)
        /// </summary>
        public const int C_長さ不明のギャップ長 = 100;

        /// <summary>
        /// 長さ不明のギャップに置き換えるとき、印の左右の片の端を削る長さ (削った分は同じ N 100 個に含める)<br/>
        /// 左右どちらの片も 2 × この長さ + 1 以上のときだけ削る
        /// </summary>
        private const int C_長さ不明のギャップの両側で削る長さ = 30;

        /// <summary>
        /// 長さ不明のギャップに置き換えるとき、両側の片どうしの写しを探す種の長さ (A・C・G・T のみ、同じ向きの完全一致だけを数える)
        /// </summary>
        private const int C_写しを探す種の長さ = 25;

        /// <summary>
        /// 長さ不明のギャップに置き換えるとき、写しを探す片の端の長さ (左の片は末尾、右の片は先頭)
        /// </summary>
        private const int C_写しを探す長さ = 300;

        /// <summary>
        /// 写しを数えるとき、一致の途中で許す違いの数 (違いに当たったとき、それまでに許した数がこれ未満のときだけ許す)
        /// </summary>
        private const int C_写しで許す違いの数 = 2;

        /// <summary>
        /// 写しを数えるとき、許した違いの直後に一致が要る塩基の数 (片の端までの残りがこれより短ければ、残りが全部一致すること)
        /// </summary>
        private const int C_違いの後に要る一致の長さ = 10;

        /// <summary>
        /// 最終成果物の構成を書き出す AGP ファイルの拡張子
        /// </summary>
        private const string C_AGP拡張子 = ".agp";

        /// <summary>
        /// 最終成果物の統計表のファイル名
        /// </summary>
        private const string C_統計表ファイル名 = "stats.md";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 最終成果物から、リード長より短い配列を落とす
        /// </summary>
        /// <param name="p_パス">対象の FASTA パス</param>
        /// <param name="p_リード長">下限として使うリード長</param>
        public static void V_除外_短い配列(string p_パス, int? p_リード長)
        {
            if (p_リード長 is not { } l_下限 || l_下限 <= 0 || !File.Exists(p_パス))
            {
                return;
            }

            var l_全件 = FastaReader.Get_全エントリ(p_パス);
            var l_残す = l_全件.Where(x => x.A_配列.Length >= l_下限).ToList();
            if (l_残す.Count == l_全件.Count)
            {
                return;
            }

            var l_落とした延長 = l_全件.Sum(x => (long)x.A_配列.Length) - l_残す.Sum(x => (long)x.A_配列.Length);
            using (var l_書き込み = new FastaWriter(p_パス))
            {
                foreach (var (l_ID, l_配列) in l_残す)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }

            Logger.V_出力(メッセージID.短い配列を除外, l_全件.Count - l_残す.Count, l_下限, l_落とした延長);
        }

        /// <summary>
        /// 未確認の繋ぎ目のうち、重なりをリードで確かめられるものを畳む
        /// </summary>
        /// <param name="p_パス群">対象の FASTA パス (最初のものの件数をログに出す)</param>
        /// <param name="p_k長">採用した k 長</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_リードパス群">生リードのパス</param>
        public static void V_畳む_リードで確かめた繋ぎ目(string[] p_パス群, int p_k長, int p_リード長, IEnumerable<string> p_リードパス群)
        {
            var l_対象 = p_パス群.Where(File.Exists).Select(x => (A_パス: x, A_全件: FastaReader.Get_全エントリ(x))).Where(x => x.A_全件.Any(y => y.A_配列.Contains(Consts.未確認の繋ぎ目))).ToList();
            if (l_対象.Count == 0 || RepeatRMerVerifier.Get_リード索引(p_リードパス群) is not { } l_索引)
            {
                return;
            }

            var l_一意の出現数 = Scaffolder.Get_一意の出現数(l_索引, l_対象[0].A_全件.Select(x => x.A_配列));
            for (var i = 0; i < l_対象.Count; i++)
            {
                var l_総数 = 0;
                var l_畳んだ数 = 0;
                var l_埋めた数 = 0;
                var l_一致で畳んだ数 = 0;
                using (var l_書き込み = new FastaWriter(l_対象[i].A_パス))
                {
                    foreach (var (l_ID, l_配列) in l_対象[i].A_全件)
                    {
                        l_書き込み.V_書き込み(l_ID, Get_確かめた繋ぎ目を畳んだ配列(l_配列, l_索引, p_k長, p_リード長, l_一意の出現数, ref l_総数, ref l_畳んだ数, ref l_埋めた数, ref l_一致で畳んだ数));
                    }
                }

                if (l_対象[i].A_パス == p_パス群[0])
                {
                    Logger.V_出力(メッセージID.未確認の繋ぎ目をリードで畳んだ, l_総数, l_畳んだ数, l_埋めた数, l_一致で畳んだ数);
                }
            }
        }

        /// <summary>
        /// 配列の未確認の繋ぎ目ごとに、次の片との重なりをリードで確かめ、ただ 1 通りに決まれば畳む<br/>
        /// 畳めなければ、両側を跨ぐリードの続きで埋められるか試す<br/>
        /// それもできなければ、両側の配列が長く完全に一致する重なりでただ 1 通りに決まるか試す
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_索引">リードの索引</param>
        /// <param name="p_k長">採用した k 長</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安 (Scaffolder.Get_一意の出現数) 、分からなければ 0</param>
        /// <param name="p_総数">未確認の繋ぎ目の数 (加算する)</param>
        /// <param name="p_畳んだ数">リードで確かめた重なりで畳んだ数 (加算する)</param>
        /// <param name="p_埋めた数">リードの続きで埋めた数 (加算する)</param>
        /// <param name="p_一致で畳んだ数">両側の配列の長い一致で畳んだ数 (加算する)</param>
        /// <returns>畳んだ後の配列</returns>
        public static string Get_確かめた繋ぎ目を畳んだ配列(string p_配列, ReadMinimizerIndex p_索引, int p_k長, int p_リード長, int p_一意の出現数, ref int p_総数, ref int p_畳んだ数, ref int p_埋めた数, ref int p_一致で畳んだ数)
        {
            var l_出力 = new StringBuilder(p_配列.Length);
            var l_インデックス = 0;
            while (l_インデックス < p_配列.Length)
            {
                if (p_配列[l_インデックス] != Consts.未確認の繋ぎ目)
                {
                    _ = l_出力.Append(p_配列[l_インデックス]);
                    l_インデックス++;
                    continue;
                }

                p_総数++;
                var l_次の終わり = l_インデックス + 1;
                while (l_次の終わり < p_配列.Length && !Util.Isギャップ文字(p_配列[l_次の終わり]))
                {
                    l_次の終わり++;
                }

                var l_次の配列 = p_配列[(l_インデックス + 1)..l_次の終わり];
                if (Scaffolder.Get_リードで確かめた重なり長(p_索引, l_出力, l_次の配列, p_k長, p_リード長, p_一意の出現数) is { } l_重なり長)
                {
                    p_畳んだ数++;
                    l_インデックス += 1 + l_重なり長;
                    continue;
                }

                if (Scaffolder.Get_リードで埋めた繋ぎ目(p_索引, l_出力, l_次の配列, p_リード長, p_一意の出現数) is { } l_埋め方)
                {
                    p_埋めた数++;
                    l_出力.Length -= l_埋め方.A_左から削る長さ;
                    _ = l_出力.Append(l_埋め方.A_埋める配列);
                    l_インデックス += 1 + l_埋め方.A_右から削る長さ;
                    continue;
                }

                if (Scaffolder.Get_配列の一致で決めた重なり長(l_出力, l_次の配列, p_k長) is { } l_一致の重なり長)
                {
                    p_一致で畳んだ数++;
                    l_インデックス += 1 + l_一致の重なり長;
                    continue;
                }

                _ = l_出力.Append(Consts.未確認の繋ぎ目);
                l_インデックス++;
            }

            return l_出力.ToString();
        }

        /// <summary>
        /// 未確認の繋ぎ目の印を、長さ不明のギャップ (N を 100 個) に置き換える<br/>
        /// 印の左右の片の端 (左は末尾、右は先頭) を削り、その分は N 100 個に含める
        /// </summary>
        /// <param name="p_パス">対象の FASTA パス</param>
        /// <param name="p_Isログ出力">置き換えた数をログに出すか</param>
        /// <returns>配列 ID ごとの、長さ不明のギャップになった N の連続の番号 (0 始まり)</returns>
        public static Dictionary<string, HashSet<int>> V_置換_未確認の繋ぎ目(string p_パス, bool p_Isログ出力)
        {
            Dictionary<string, HashSet<int>> l_長さ不明の番号 = [];
            if (!File.Exists(p_パス))
            {
                return l_長さ不明の番号;
            }

            var l_全件 = FastaReader.Get_全エントリ(p_パス);
            if (!l_全件.Any(x => x.A_配列.Contains(Consts.未確認の繋ぎ目)))
            {
                return l_長さ不明の番号;
            }

            var l_置換数 = 0;
            using (var l_書き込み = new FastaWriter(p_パス))
            {
                foreach (var (l_ID, l_配列) in l_全件)
                {
                    var l_番号群 = new HashSet<int>();
                    l_書き込み.V_書き込み(l_ID, Get_長さ不明のギャップへ置換(l_配列, l_番号群));
                    if (l_番号群.Count > 0)
                    {
                        l_長さ不明の番号[Get_配列名(l_ID)] = l_番号群;
                        l_置換数 += l_番号群.Count;
                    }
                }
            }

            if (p_Isログ出力)
            {
                Logger.V_出力(メッセージID.長さ不明のギャップへ置換, l_置換数, C_長さ不明のギャップ長);
            }

            return l_長さ不明の番号;
        }

        /// <summary>
        /// 未確認の繋ぎ目の印を長さ不明のギャップに置き換え、それが何番目の N の連続になったかを集める<br/>
        /// 左の片は出力済みの配列のうち最後のギャップ文字より後ろ、右の片は印の次から次のギャップ文字の手前まで
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_長さ不明の番号">長さ不明のギャップになった N の連続の番号 (0 始まり) を足す先</param>
        /// <returns>置き換えた配列</returns>
        public static string Get_長さ不明のギャップへ置換(string p_配列, HashSet<int> p_長さ不明の番号)
        {
            var l_出力 = new StringBuilder(p_配列.Length);
            var l_連続の番号 = -1;
            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_文字 = p_配列[i];
                if (Util.Isギャップ文字(l_文字) && (i == 0 || !Util.Isギャップ文字(p_配列[i - 1])))
                {
                    l_連続の番号++;
                }

                if (l_文字 == Consts.未確認の繋ぎ目)
                {
                    var l_要る片の長さ = (2 * C_長さ不明のギャップの両側で削る長さ) + 1;
                    var l_左の片の長さ = Get_末尾の片の長さ(l_出力);
                    var l_右の片の長さ = Get_先頭の片の長さ(p_配列, i + 1);
                    var l_削る長さ = l_左の片の長さ >= l_要る片の長さ && l_右の片の長さ >= l_要る片の長さ ? C_長さ不明のギャップの両側で削る長さ : 0;
                    l_出力.Length -= l_削る長さ;

                    var l_右の始め = i + 1 + l_削る長さ;
                    var l_読み飛ばす長さ = Get_右の片の先頭の写しの長さ(l_出力, p_配列, l_右の始め, out var l_畳む);
                    if (l_読み飛ばす長さ == 0)
                    {
                        l_出力.Length -= Get_左の片の末尾の写しの長さ(l_出力, p_配列, l_右の始め);
                    }

                    i += l_削る長さ + l_読み飛ばす長さ;
                    if (l_畳む)
                    {
                        // 畳んだ印は N の連続にならないので、後ろの番号を 1 つ戻す
                        l_連続の番号--;
                    }
                    else
                    {
                        _ = l_出力.Append('N', C_長さ不明のギャップ長);
                        _ = p_長さ不明の番号.Add(l_連続の番号);
                    }
                }
                else
                {
                    _ = l_出力.Append(l_文字);
                }
            }

            return l_出力.ToString();
        }

        /// <summary>
        /// 最終成果物の構成を AGP (v2.1) で書き出す
        /// </summary>
        /// <param name="p_FASTAパス">最終成果物</param>
        /// <param name="p_長さ不明の番号">配列 ID ごとの、長さ不明のギャップの N の連続の番号</param>
        /// <param name="p_AGPパス">書き出す先</param>
        public static void V_書き出し_AGP(string p_FASTAパス, Dictionary<string, HashSet<int>> p_長さ不明の番号, string p_AGPパス)
        {
            if (!File.Exists(p_FASTAパス))
            {
                return;
            }

            using var l_書き込み = new StreamWriter(p_AGPパス);
            l_書き込み.WriteLine(C_見出し_agp_version_2_1);
            foreach (var (l_ID, l_配列) in FastaReader.Get_全エントリ(p_FASTAパス))
            {
                var l_名前 = Get_配列名(l_ID);
                foreach (var l_行 in Get_AGP行(l_名前, l_配列, p_長さ不明の番号.GetValueOrDefault(l_名前)))
                {
                    l_書き込み.WriteLine(l_行);
                }
            }
        }

        /// <summary>
        /// 1 本の配列の AGP の行
        /// </summary>
        /// <param name="p_名前">配列名</param>
        /// <param name="p_配列">配列</param>
        /// <param name="p_長さ不明の番号">長さ不明のギャップの N の連続の番号、無ければ null</param>
        /// <returns>AGP の行</returns>
        public static IEnumerable<string> Get_AGP行(string p_名前, string p_配列, IReadOnlySet<int>? p_長さ不明の番号)
        {
            var l_部品番号 = 1;
            var l_片番号 = 1;
            var l_連続の番号 = 0;
            var l_位置 = 0;
            while (l_位置 < p_配列.Length)
            {
                var l_終わり = l_位置;
                var l_Isギャップ = p_配列[l_位置] == 'N';
                while (l_終わり < p_配列.Length && (p_配列[l_終わり] == 'N') == l_Isギャップ)
                {
                    l_終わり++;
                }

                var l_長さ = l_終わり - l_位置;
                if (l_Isギャップ)
                {
                    var l_種類 = p_長さ不明の番号?.Contains(l_連続の番号) == true ? C_項目_U : C_長さ指定ギャップ種別;
                    yield return FormattableString.Invariant($"{p_名前}	{l_位置 + 1}	{l_終わり}	{l_部品番号}	{l_種類}	{l_長さ}	scaffold	yes	paired-ends");
                    l_連続の番号++;
                }
                else
                {
                    yield return FormattableString.Invariant($"{p_名前}	{l_位置 + 1}	{l_終わり}	{l_部品番号}	W	{p_名前}_{l_片番号}	1	{l_長さ}	+");
                    l_片番号++;
                }

                l_部品番号++;
                l_位置 = l_終わり;
            }
        }

        /// <summary>
        /// 最終配列を整形して元リードと照合し、レポートを出力する
        /// </summary>
        /// <param name="p_結果">採用したアセンブリ</param>
        /// <param name="p_原入力">加工前のリードを保持した設定</param>
        /// <param name="p_一時ディレクトリ">成果物の出力先</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_処理済み設定">前処理・エラー訂正後のパスを保持した設定、無ければ null</param>
        public static void V_実行(アセンブリ実行結果 p_結果, Parameters p_原入力, string p_一時ディレクトリ, int? p_リード長, Parameters? p_処理済み設定 = null)
        {
            var l_最終パス = AssemblyPipeline.V_複製_最終成果物(p_結果, p_一時ディレクトリ);
            V_除外_短い配列(l_最終パス, p_リード長);
            V_繋ぐ_ペアで橋渡し(l_最終パス, p_原入力);
            var l_scaffoldパス = Path.Combine(p_一時ディレクトリ, Consts.Scaffoldファイル名);
            V_畳む_リードで確かめた繋ぎ目([l_最終パス, l_scaffoldパス], p_結果.A_k長, p_リード長 ?? 0, AssemblyPipeline.Get_全リードパス(p_原入力));
            if (p_リード長 > 0)
            {
                var l_一意の出現数 = V_繋ぐ_リードで跨げる端(l_最終パス, p_リード長 ?? 0, AssemblyPipeline.Get_全リードパス(p_原入力));
                V_繋ぐ_余ったリードで端(l_最終パス, p_リード長 ?? 0, p_原入力, l_一意の出現数);
                V_伸ばす_端(l_最終パス, p_リード長 ?? 0, AssemblyPipeline.Get_全リードパス(p_原入力), l_一意の出現数);
            }

            RepeatRMerVerifier.V_解放_共有索引();
            var l_長さ不明の番号 = V_置換_未確認の繋ぎ目(l_最終パス, p_Isログ出力: true);
            _ = V_置換_未確認の繋ぎ目(l_scaffoldパス, p_Isログ出力: false);
            var l_ポリッシュ統計 = V_磨く(p_原入力, p_一時ディレクトリ, l_最終パス);
            V_書き出し_AGP(l_最終パス, l_長さ不明の番号, Path.ChangeExtension(l_最終パス, C_AGP拡張子));
            V_評価_継ぎ目(p_原入力, p_一時ディレクトリ, l_最終パス);
            var l_閉鎖検証 = V_検証_環状閉鎖(p_原入力, l_最終パス);
            var l_支持検査 = V_検査_リード支持(p_原入力, l_最終パス);
            p_結果 = p_結果 with
            {
                A_整合性検査 = Get_最終整合性(p_原入力, p_結果.A_k長, l_最終パス, p_一時ディレクトリ)
            };
            V_記録_出所(p_原入力, p_結果, l_最終パス, p_一時ディレクトリ, p_処理済み設定);
            V_出力_最終統計(p_一時ディレクトリ, l_最終パス);
            V_出力_完全性レポート(p_結果, p_原入力, l_最終パス, l_ポリッシュ統計, l_閉鎖検証, l_支持検査, p_一時ディレクトリ);
            Logger.V_出力(メッセージID.最終成果物, l_最終パス);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 出力済みの配列で、最後のギャップ文字より後ろの片の長さ
        /// </summary>
        /// <param name="p_出力">出力済みの配列</param>
        /// <returns>末尾の片の長さ (ギャップ文字が無ければ全体の長さ)</returns>
        private static int Get_末尾の片の長さ(StringBuilder p_出力)
        {
            var l_長さ = 0;
            while (l_長さ < p_出力.Length && !Util.Isギャップ文字(p_出力[p_出力.Length - 1 - l_長さ]))
            {
                l_長さ++;
            }

            return l_長さ;
        }

        /// <summary>
        /// 指定した位置から始まる片の長さ (次のギャップ文字の手前まで)
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_始め">片の始まりの位置</param>
        /// <returns>先頭の片の長さ (ギャップ文字が無ければ末尾までの長さ)</returns>
        private static int Get_先頭の片の長さ(string p_配列, int p_始め)
        {
            var l_終わり = p_始め;
            while (l_終わり < p_配列.Length && !Util.Isギャップ文字(p_配列[l_終わり]))
            {
                l_終わり++;
            }

            return l_終わり - p_始め;
        }

        /// <summary>
        /// 出力済みの配列の末尾の片 (最後のギャップ文字より後ろ) の配列
        /// </summary>
        /// <param name="p_出力">出力済みの配列</param>
        /// <returns>末尾の片の配列</returns>
        private static string Get_末尾の片(StringBuilder p_出力)
        {
            var l_長さ = Get_末尾の片の長さ(p_出力);
            return p_出力.ToString(p_出力.Length - l_長さ, l_長さ);
        }

        /// <summary>
        /// 右の片の先頭が、左の片の末尾の配列の写しなら、右の片の先頭から読み飛ばす長さを返す<br/>
        /// 左の片の末尾 C_写しを探す長さ 塩基の中に種 (右の片の先頭 C_写しを探す種の長さ 塩基) がちょうど 1 か所だけあり、そこから違いを許して同じ間の長さ ℓ を数えて、右の片の残りが種の長さ以上のときだけ ℓ を返す<br/>
        /// ℓ が左の片の末尾まで届いたときは、N を書かずに畳む (p_畳む を true にする)
        /// </summary>
        /// <param name="p_出力">出力済みの配列 (左の片は末尾の片)</param>
        /// <param name="p_配列">配列</param>
        /// <param name="p_右の始め">右の片の始まりの位置</param>
        /// <param name="p_畳む">ℓ が左の片の末尾まで届き、N を書かずに畳むなら true</param>
        /// <returns>読み飛ばす長さ、写しと判らなければ 0</returns>
        private static int Get_右の片の先頭の写しの長さ(StringBuilder p_出力, string p_配列, int p_右の始め, out bool p_畳む)
        {
            p_畳む = false;
            var l_左 = Get_末尾の片(p_出力);
            var l_右 = p_配列.Substring(p_右の始め, Get_先頭の片の長さ(p_配列, p_右の始め));
            if (l_左.Length < C_写しを探す種の長さ || l_右.Length < C_写しを探す種の長さ)
            {
                return 0;
            }

            var l_種 = l_右[..C_写しを探す種の長さ];
            if (!Is_ACGTのみ(l_種))
            {
                return 0;
            }

            var l_窓の始め = Math.Max(0, l_左.Length - C_写しを探す長さ);
            var l_位置 = Get_窓で一意の位置(l_左[l_窓の始め..], l_種);
            if (l_位置 < 0)
            {
                return 0;
            }

            var l_左の始め = l_窓の始め + l_位置;
            var l_一致の長さ = Get_違いを許す一致の長さ(l_左, l_左の始め, l_右, 0, 1);
            if (l_右.Length - l_一致の長さ < C_写しを探す種の長さ)
            {
                return 0;
            }

            p_畳む = l_左の始め + l_一致の長さ == l_左.Length;
            return l_一致の長さ;
        }

        /// <summary>
        /// 左の片の末尾が、右の片の先頭の配列の写しなら、左の片の末尾から除く長さを返す<br/>
        /// 左の片の末尾の種 (C_写しを探す種の長さ 塩基) が右の片の先頭 C_写しを探す長さ 塩基の中にちょうど 1 か所だけあり、そこから後ろ向きに違いを許して同じ間の長さ ℓ を数えて、左の片の残りが種の長さ以上のときだけ ℓ を返す
        /// </summary>
        /// <param name="p_出力">出力済みの配列 (左の片は末尾の片)</param>
        /// <param name="p_配列">配列</param>
        /// <param name="p_右の始め">右の片の始まりの位置</param>
        /// <returns>除く長さ、写しと判らなければ 0</returns>
        private static int Get_左の片の末尾の写しの長さ(StringBuilder p_出力, string p_配列, int p_右の始め)
        {
            var l_左 = Get_末尾の片(p_出力);
            var l_右 = p_配列.Substring(p_右の始め, Get_先頭の片の長さ(p_配列, p_右の始め));
            if (l_左.Length < C_写しを探す種の長さ || l_右.Length < C_写しを探す種の長さ)
            {
                return 0;
            }

            var l_種 = l_左[^C_写しを探す種の長さ..];
            if (!Is_ACGTのみ(l_種))
            {
                return 0;
            }

            var l_位置 = Get_窓で一意の位置(l_右[..Math.Min(C_写しを探す長さ, l_右.Length)], l_種);
            if (l_位置 < 0)
            {
                return 0;
            }

            var l_右の終わり = l_位置 + C_写しを探す種の長さ - 1;
            var l_一致の長さ = Get_違いを許す一致の長さ(l_左, l_左.Length - 1, l_右, l_右の終わり, -1);
            return l_左.Length - l_一致の長さ >= C_写しを探す種の長さ ? l_一致の長さ : 0;
        }

        /// <summary>
        /// 2 つの配列を、指定した位置から同じ向きに 1 塩基ずつ比べ、一致の長さ ℓ を返す<br/>
        /// 違う塩基に当たったとき、それまでに許した違いが C_写しで許す違いの数 未満で、その直後の C_違いの後に要る一致の長さ 塩基がすべて一致するなら、その違いを許して数え続ける (片の端までの残りがそれより短ければ、残りがすべて一致すること。残りが 0 なら許さない)<br/>
        /// 満たさなければ、その違いの手前で止める (ℓ は最後に一致した塩基までの長さ)
        /// </summary>
        /// <param name="p_左">左の配列</param>
        /// <param name="p_左の始め">左の配列で比べ始める位置</param>
        /// <param name="p_右">右の配列</param>
        /// <param name="p_右の始め">右の配列で比べ始める位置</param>
        /// <param name="p_向き">比べる向き (1 は後ろ、-1 は前)</param>
        /// <returns>一致の長さ ℓ、比べ始めで一致しなければ 0</returns>
        private static int Get_違いを許す一致の長さ(string p_左, int p_左の始め, string p_右, int p_右の始め, int p_向き)
        {
            var l_一致の長さ = 0;
            var l_違いの数 = 0;
            var l_左 = p_左の始め;
            var l_右 = p_右の始め;
            while (Is_範囲内(p_左, l_左) && Is_範囲内(p_右, l_右))
            {
                if (p_左[l_左] == p_右[l_右])
                {
                    l_一致の長さ = ((l_左 - p_左の始め) * p_向き) + 1;
                }
                else if (l_違いの数 < C_写しで許す違いの数 && Is_違いの後が一致(p_左, l_左 + p_向き, p_右, l_右 + p_向き, p_向き))
                {
                    l_違いの数++;
                }
                else
                {
                    break;
                }

                l_左 += p_向き;
                l_右 += p_向き;
            }

            return l_一致の長さ;
        }

        /// <summary>
        /// 違いの直後から、C_違いの後に要る一致の長さ 塩基 (片の端で尽きれば残り全部) が一致するか
        /// </summary>
        /// <param name="p_左">左の配列</param>
        /// <param name="p_左の始め">左の配列で、違いの直後の位置</param>
        /// <param name="p_右">右の配列</param>
        /// <param name="p_右の始め">右の配列で、違いの直後の位置</param>
        /// <param name="p_向き">比べる向き (1 は後ろ、-1 は前)</param>
        /// <returns>一致するなら true (比べられる塩基が 0 個なら false)</returns>
        private static bool Is_違いの後が一致(string p_左, int p_左の始め, string p_右, int p_右の始め, int p_向き)
        {
            var l_数 = 0;
            while (l_数 < C_違いの後に要る一致の長さ && Is_範囲内(p_左, p_左の始め + (p_向き * l_数)) && Is_範囲内(p_右, p_右の始め + (p_向き * l_数)))
            {
                if (p_左[p_左の始め + (p_向き * l_数)] != p_右[p_右の始め + (p_向き * l_数)])
                {
                    return false;
                }

                l_数++;
            }

            return l_数 > 0;
        }

        /// <summary>
        /// 位置が配列の中にあるか
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <param name="p_位置">位置</param>
        /// <returns>中にあれば true</returns>
        private static bool Is_範囲内(string p_配列, int p_位置)
        {
            return p_位置 >= 0 && p_位置 < p_配列.Length;
        }

        /// <summary>
        /// 窓の中に種がちょうど 1 か所だけ出てくるとき、その位置を返す (重なる出現も数える)
        /// </summary>
        /// <param name="p_窓">探す範囲</param>
        /// <param name="p_種">探す配列</param>
        /// <returns>窓の中の位置、出現が 1 か所でなければ -1</returns>
        private static int Get_窓で一意の位置(string p_窓, string p_種)
        {
            var l_位置 = -1;
            var l_数 = 0;
            for (var i = 0; i + p_種.Length <= p_窓.Length; i++)
            {
                if (string.CompareOrdinal(p_窓, i, p_種, 0, p_種.Length) == 0)
                {
                    l_位置 = i;
                    l_数++;
                }
            }

            return l_数 == 1 ? l_位置 : -1;
        }

        /// <summary>
        /// 配列が A・C・G・T だけからなるか
        /// </summary>
        /// <param name="p_配列">配列</param>
        /// <returns>A・C・G・T だけなら true</returns>
        private static bool Is_ACGTのみ(string p_配列)
        {
            return p_配列.All(x => x is 'A' or 'C' or 'G' or 'T');
        }

        /// <summary>
        /// 配列の端どうしを、跨ぐリード 2 本以上で繋ぐ<br/>
        /// 繋いだ数が 1 以上なら、同じパスに書き戻す
        /// </summary>
        /// <param name="p_最終パス">最終アセンブリの FASTA パス</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_リードパス群">生リードのパス</param>
        /// <returns>算出した一意の出現数、索引が作れなければ 0</returns>
        private static int V_繋ぐ_リードで跨げる端(string p_最終パス, int p_リード長, IEnumerable<string> p_リードパス群)
        {
            if (RepeatRMerVerifier.Get_リード索引(p_リードパス群) is not { } l_索引)
            {
                return 0;
            }

            var l_配列群 = FastaReader.Get_全エントリ(p_最終パス);
            var l_一意の出現数 = Scaffolder.Get_一意の出現数(l_索引, l_配列群.Select(x => x.A_配列));
            var l_繋いだ配列群 = ReadEndJoiner.Get_繋いだ配列群(l_配列群, l_索引, p_リード長, l_一意の出現数, out var l_候補数, out var l_繋いだ数);
            if (l_繋いだ数 >= 1)
            {
                using var l_書き込み = new FastaWriter(p_最終パス);
                foreach (var (l_ID, l_配列) in l_繋いだ配列群)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }

            Logger.V_出力(メッセージID.リードで端を繋いだ, l_候補数, l_繋いだ数, l_配列群.Count, l_繋いだ配列群.Count);
            return l_一意の出現数;
        }

        /// <summary>
        /// 間引く前の全リードのうち、配列の端に触れるものだけで、もう一度端をつなぐ<br/>
        /// 間引いたときだけ行う。繋いだ数が 1 以上なら、同じパスに書き戻す
        /// </summary>
        /// <param name="p_最終パス">最終アセンブリの FASTA パス</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_原入力">間引く前のライブラリ群と採用割合を保持した設定</param>
        /// <param name="p_一意の出現数">間引いたリードから算出した一意の出現数</param>
        private static void V_繋ぐ_余ったリードで端(string p_最終パス, int p_リード長, Parameters p_原入力, int p_一意の出現数)
        {
            var l_採用割合 = p_原入力.A_間引きの採用割合;
            if (p_原入力.A_間引き前ライブラリ群.Count == 0 || l_採用割合 <= 0 || l_採用割合 >= 1 || p_一意の出現数 <= 0)
            {
                return;
            }

            using var l_計測 = new StageTimer("read-join-surplus");
            var l_配列群 = FastaReader.Get_全エントリ(p_最終パス);
            var l_種集合 = Get_端の種集合(l_配列群.Where(x => x.A_配列.Length >= C_端つなぎの最短長 && !x.A_ID.Contains(Consts.環状の目印)).Select(x => x.A_配列));
            if (l_種集合.Count == 0)
            {
                return;
            }

            var l_端リード群 = new List<string>();
            foreach (var (l_順リード, l_逆リード) in p_原入力.A_間引き前ライブラリ群)
            {
                l_端リード群.AddRange(Get_端に触れるリード(l_種集合, FastqReader.Get_生リード列(l_順リード, l_逆リード)));
            }

            if (l_端リード群.Count == 0)
            {
                return;
            }

            var l_索引 = ReadMinimizerIndex.V_構築(() => l_端リード群);
            var l_目安 = (int)Math.Round(p_一意の出現数 / l_採用割合);
            var l_繋いだ配列群 = ReadEndJoiner.Get_繋いだ配列群(l_配列群, l_索引, p_リード長, l_目安, out var l_候補数, out var l_繋いだ数, p_要るリード数: C_余りで要るリード数);
            if (l_繋いだ数 >= 1)
            {
                using var l_書き込み = new FastaWriter(p_最終パス);
                foreach (var (l_ID, l_配列) in l_繋いだ配列群)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }

            Logger.V_出力(メッセージID.余ったリードで端を繋いだ, l_端リード群.Count, l_候補数, l_繋いだ数, l_配列群.Count, l_繋いだ配列群.Count);
        }

        /// <summary>
        /// 配列の端を、端の手前の錨を含むリードの続きの合意で伸ばす<br/>
        /// 伸ばした端が 1 以上なら、同じパスに書き戻す
        /// </summary>
        /// <param name="p_最終パス">最終アセンブリの FASTA パス</param>
        /// <param name="p_リード長">代表リード長</param>
        /// <param name="p_リードパス群">生リードのパス</param>
        /// <param name="p_一意の出現数">一意な配列がリードに出てくる数の目安</param>
        private static void V_伸ばす_端(string p_最終パス, int p_リード長, IEnumerable<string> p_リードパス群, int p_一意の出現数)
        {
            using var l_計測 = new StageTimer("end-extension");
            if (RepeatRMerVerifier.Get_リード索引(p_リードパス群) is not { } l_索引)
            {
                return;
            }

            var l_配列群 = FastaReader.Get_全エントリ(p_最終パス);
            var l_伸ばした配列群 = ReadEndExtender.Get_伸ばした配列群(l_配列群, l_索引, p_リード長, p_一意の出現数, out var l_伸ばした端の数, out var l_伸ばした塩基数);
            if (l_伸ばした端の数 >= 1)
            {
                using var l_書き込み = new FastaWriter(p_最終パス);
                foreach (var (l_ID, l_配列) in l_伸ばした配列群)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }

                Logger.V_出力(メッセージID.端をリードで伸ばした, l_伸ばした端の数, l_伸ばした塩基数);
            }
        }

        /// <summary>
        /// ペアで端どうしの相手を決め、間を局所組み立てで埋めて繋ぐ。埋まらなくても票が多ければ未確認の繋ぎ目で繋ぐ<br/>
        /// ライブラリは間引く前の群があればそれ、無ければ間引かないライブラリ群を使う<br/>
        /// 配列の数が減ったときだけ、同じパスに書き戻す
        /// </summary>
        /// <param name="p_最終パス">最終アセンブリの FASTA パス</param>
        /// <param name="p_原入力">ライブラリ群を保持した設定</param>
        private static void V_繋ぐ_ペアで橋渡し(string p_最終パス, Parameters p_原入力)
        {
            var l_ライブラリ群 = p_原入力.A_間引き前ライブラリ群.Count > 0 ? p_原入力.A_間引き前ライブラリ群 : p_原入力.A_ライブラリ群;
            if (l_ライブラリ群.Count == 0)
            {
                return;
            }

            using var l_計測 = new StageTimer("pair-bridge");
            var l_配列群 = FastaReader.Get_全エントリ(p_最終パス);
            var l_繋いだ配列群 = PairBridgeFiller.Get_繋いだ配列群(l_配列群, l_ライブラリ群, out var l_統計);
            if (l_繋いだ配列群.Count < l_配列群.Count)
            {
                using var l_書き込み = new FastaWriter(p_最終パス);
                foreach (var (l_ID, l_配列) in l_繋いだ配列群)
                {
                    l_書き込み.V_書き込み(l_ID, l_配列);
                }
            }

            Logger.V_出力(メッセージID.ペアで橋渡しした, l_統計.A_挿入長, l_統計.A_候補の辺の数, l_統計.A_埋めて繋いだ辺の数, l_統計.A_一意でなかった数, l_統計.A_届かなかった数, l_配列群.Count, l_繋いだ配列群.Count, l_統計.A_印で繋いだ辺の数, l_統計.A_反復の端で除いた数, l_統計.A_重なりで畳んだ辺の数, l_統計.A_多コピーの端で除いた数);
        }

        /// <summary>
        /// 配列の先頭と末尾（各 300 bp、配列が 600 bp 未満なら全体）から、長さ 25 の正準キーの集合を作る
        /// </summary>
        /// <param name="p_配列群">端を見る配列の列</param>
        /// <returns>正準キーの集合</returns>
        internal static HashSet<ulong> Get_端の種集合(IEnumerable<string> p_配列群)
        {
            HashSet<ulong> l_種集合 = [];
            foreach (var l_配列 in p_配列群)
            {
                string[] l_区間群 = l_配列.Length < 2 * C_端の窓長 ? [l_配列] : [l_配列[..C_端の窓長], l_配列[^C_端の窓長..]];
                foreach (var l_区間 in l_区間群)
                {
                    var l_窓 = new RollingKmer(C_端の種長);
                    foreach (var l_塩基 in l_区間)
                    {
                        if (l_窓.Is成功_追加(l_塩基, out var l_キー))
                        {
                            _ = l_種集合.Add((ulong)l_キー.A_下位);
                        }
                    }
                }
            }

            return l_種集合;
        }

        /// <summary>
        /// リードのうち、端の種集合と 25 塩基を共有するものを読み順のまま集める
        /// </summary>
        /// <param name="p_種集合">端の正準キーの集合</param>
        /// <param name="p_リード群">調べるリードの列</param>
        /// <returns>端に触れるリードの列</returns>
        internal static List<string> Get_端に触れるリード(HashSet<ulong> p_種集合, IEnumerable<string> p_リード群)
        {
            List<string> l_結果 = [];
            foreach (var l_リード in p_リード群)
            {
                if (LocalAssembler.Has種一致(l_リード, p_種集合, C_端の種長))
                {
                    l_結果.Add(l_リード);
                }
            }

            return l_結果;
        }

        /// <summary>
        /// 作業ディレクトリ直下に出した最終成果物の統計をログに出す
        /// </summary>
        /// <param name="p_一時ディレクトリ">成果物の出力先</param>
        /// <param name="p_最終パス">除外・ポリッシュを終えた assembly.fasta</param>
        private static void V_出力_最終統計(string p_一時ディレクトリ, string p_最終パス)
        {
            Logger.V_出力_空行();
            foreach (var (l_ラベル, l_パス) in Get_最終成果物群(p_一時ディレクトリ, p_最終パス))
            {
                if (File.Exists(l_パス))
                {
                    AssemblyStatsReporter.V_出力_統計($"final {l_ラベル}", l_パス);
                }
            }
        }

        /// <summary>
        /// 統計を出す最終成果物のラベルとパス
        /// </summary>
        /// <param name="p_一時ディレクトリ">成果物の出力先</param>
        /// <param name="p_最終パス">除外・ポリッシュを終えた assembly.fasta</param>
        /// <returns></returns>
        private static (string A_ラベル, string A_FASTAパス)[] Get_最終成果物群(string p_一時ディレクトリ, string p_最終パス)
        {
            return [
                (C_項目_contigs, Path.Combine(p_一時ディレクトリ, AssemblyPipeline.C_Contigファイル名)),
                (C_項目_scaffolds, Path.Combine(p_一時ディレクトリ, Consts.Scaffoldファイル名)),
                (C_項目_assembly, p_最終パス),
            ];
        }

        /// <summary>
        /// 加工前のリードから最終配列の整合性を再計算する
        /// </summary>
        /// <param name="p_原入力">加工前の設定</param>
        /// <param name="p_k長">採用された k 長</param>
        /// <param name="p_最終パス">最終配列</param>
        /// <param name="p_作業パス">検査用ファイルの親ディレクトリ</param>
        /// <returns>基準値を推定できた場合の整合性</returns>
        private static 整合性検査結果? Get_最終整合性(Parameters p_原入力, int p_k長, string p_最終パス, string p_作業パス)
        {
            var l_元設定 = ConfigurationManager.A_実行時引数;
            var l_元モデル = ConfigurationManager.A_スペクトルモデル;
            var l_設定 = p_原入力.Get_複製();
            l_設定.Set_推定k長(p_k長);
            l_設定.A_Is曖昧塩基許容 = false;
            try
            {
                ConfigurationManager.A_実行時引数 = l_設定;
                var l_検査パス = Path.Combine(p_作業パス, AssemblyWorkspace.C_検査ディレクトリ名);
                _ = Directory.CreateDirectory(l_検査パス);
                using var l_索引 = new TrustedKmerIndex(l_検査パス);
                KmerCounting.V_読込_リードペア(l_設定, l_索引, p_Is進行状況出力: false);
                var l_分布 = l_索引.Get_出現回数ヒストグラム();
                var l_基準 = KmerSpectrumMixtureModel.Get_解析結果(l_分布)?.A_単一コピー平均 ?? KmerHistogram.Get_解析結果(l_分布)?.A_単一コピー基準値 ?? 0D;
                KmerCutoffSelector.V_解決_kmerカットオフ(l_設定, l_索引);
                l_索引.V_適用_カットオフ(l_設定.A_kmerカットオフ);
                return l_基準 > 0D ? AssemblyValidator.Get_検査結果(p_最終パス, l_索引, p_k長, l_基準) : null;
            }
            finally
            {
                ConfigurationManager.A_実行時引数 = l_元設定;
                ConfigurationManager.A_スペクトルモデル = l_元モデル;
            }
        }

        /// <summary>
        /// 検査に使った原入力と最終配列の出所を記録する
        /// </summary>
        /// <param name="p_原入力">加工前の設定</param>
        /// <param name="p_結果"></param>
        /// <param name="p_最終パス">最終配列</param>
        /// <param name="p_作業パス">出力先</param>
        /// <param name="p_処理済み設定">前処理・エラー訂正後のパスを保持した設定、無ければ null</param>
        private static void V_記録_出所(Parameters p_原入力, アセンブリ実行結果 p_結果, string p_最終パス, string p_作業パス, Parameters? p_処理済み設定)
        {
            using var l_出力 = File.Create(Path.Combine(p_作業パス, C_ファイル名_assembly_provenance_json));
            using var l_JSON = new System.Text.Json.Utf8JsonWriter(l_出力, new System.Text.Json.JsonWriterOptions { Indented = true });
            l_JSON.WriteStartObject();
            l_JSON.WriteNumber(C_項目_schema_version, 1);
            l_JSON.WriteString(C_項目_build_id, typeof(FinalAssemblyPipeline).Assembly.ManifestModule.ModuleVersionId);
            l_JSON.WriteString(C_項目_validation_source, C_項目_uncorrected_reads_from_t_ent_holdout_data);
            l_JSON.WriteString(C_項目_assembly_sha256, StageCheckpoint.Get_ハッシュ(p_最終パス));
            l_JSON.WriteNumber(C_項目_thread_count, p_原入力.A_スレッド数);
            l_JSON.WriteString(C_項目_settings, p_原入力.ToString());
            l_JSON.WriteStartObject(C_項目_assembly_settings);
            l_JSON.WriteString(C_項目_copy_number_baseline_requested, p_原入力.A_コピー数基準の出所.ToString());
            l_JSON.WriteString(C_項目_copy_number_baseline_actual, p_結果.A_実際のコピー数基準.ToString());
            l_JSON.WriteBoolean(C_項目_trim_low_coverage_ends, p_原入力.A_Is低カバレッジ端トリミング);
            l_JSON.WriteEndObject();
            l_JSON.WriteStartArray(C_項目_inputs);
            foreach (var l_入力 in p_原入力.A_ライブラリ群.SelectMany(x => new[] { x.A_順リード, x.A_逆リード }))
            {
                if (string.IsNullOrWhiteSpace(l_入力))
                {
                    continue;
                }

                l_JSON.WriteStartObject();
                l_JSON.WriteString(C_項目_path, Path.GetFullPath(l_入力));
                l_JSON.WriteString(C_項目_sha256, StageCheckpoint.Get_ハッシュ(l_入力));
                l_JSON.WriteEndObject();
            }

            l_JSON.WriteEndArray();
            if (p_処理済み設定 is { } l_処理済み設定)
            {
                l_JSON.WriteStartArray(C_項目_corrected_read_hashes);
                foreach (var l_入力 in l_処理済み設定.A_ライブラリ群.SelectMany(x => new[] { x.A_順リード, x.A_逆リード }))
                {
                    if (!中間データ置き場.Is存在(l_入力))
                    {
                        continue;
                    }

                    l_JSON.WriteStartObject();
                    l_JSON.WriteString(C_項目_path, Path.GetFullPath(l_入力));
                    l_JSON.WriteString(C_項目_sha256, StageCheckpoint.Get_ハッシュ(l_入力));
                    l_JSON.WriteEndObject();
                }

                l_JSON.WriteEndArray();
                l_JSON.WriteString(C_項目_pipeline_fingerprint, StageCheckpoint.Get_入力署名(l_処理済み設定));
            }

            l_JSON.WriteEndObject();
        }

        /// <summary>
        /// 最終成果物の各位置がリードに裏付けられているかを調べる
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_最終パス">検査対象の最終成果物パス</param>
        /// <returns>支持検査の結果</returns>
        private static 支持検査結果? V_検査_リード支持(Parameters p_引数, string p_最終パス)
        {
            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.支持検査の開始, C_支持検査のr長);
            var l_結果 = ReadSupportChecker.Get_検査結果(p_最終パス, p_引数.A_ライブラリ群, C_支持検査のr長);
            ReadSupportChecker.V_出力_検査結果(l_結果);
            return l_結果;
        }

        /// <summary>
        /// 最終成果物にリードを貼り直して磨く
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_一時ディレクトリ">一時ディレクトリ</param>
        /// <param name="p_最終パス">磨く対象の最終成果物パス</param>
        /// <returns>ポリッシュ統計、実行しなければ null</returns>
        private static ポリッシュ統計? V_磨く(Parameters p_引数, string p_一時ディレクトリ, string p_最終パス)
        {
            if (!p_引数.A_Isポリッシュ)
            {
                return null;
            }

            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.ポリッシュ開始);
            var l_出力先 = Path.Combine(p_一時ディレクトリ, C_ポリッシュ済みファイル名);
            var l_統計 = Polisher.Get_磨いた結果(p_最終パス, p_引数.A_ライブラリ群, l_出力先);
            Polisher.V_出力_統計(l_統計);
            if (l_統計 is not null)
            {
                File.Copy(l_出力先, p_最終パス, overwrite: true);
            }

            Logger.V_出力_タイムスタンプ();
            return l_統計;
        }

        /// <summary>
        /// 最終成果物の継ぎ目ごとに、繋ぐ相手を誤っている確率を出し、確率の高い継ぎ目で切った安全版も書き出す
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_一時ディレクトリ">成果物の出力先</param>
        /// <param name="p_最終パス">評価する最終成果物パス</param>
        private static void V_評価_継ぎ目(Parameters p_引数, string p_一時ディレクトリ, string p_最終パス)
        {
            Logger.V_出力_空行();
            if (JunctionRiskEvaluator.Get_評価結果(p_最終パス, p_引数.A_ライブラリ群, Path.Combine(p_一時ディレクトリ, Consts.継ぎ目ファイル名)) is { } l_評価群)
            {
                JunctionRiskEvaluator.V_書き出し_切った配列(p_最終パス, l_評価群, Path.Combine(p_一時ディレクトリ, Consts.安全版ファイル名));
            }

            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// 環状の閉じ目を元リードで確かめる
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_最終パス">検証対象の最終成果物パス</param>
        /// <returns>検証結果、実行しなければ null</returns>
        private static IReadOnlyList<環状閉鎖検証結果>? V_検証_環状閉鎖(Parameters p_引数, string p_最終パス)
        {
            if (!p_引数.A_Is環状閉鎖検証)
            {
                return null;
            }

            Logger.V_出力_空行();
            var l_検証 = CircularClosureVerifier.Get_検証結果(p_最終パス, p_引数.A_ライブラリ群);
            CircularClosureVerifier.V_出力_検証結果(l_検証);
            Logger.V_出力_タイムスタンプ();
            return l_検証;
        }

        /// <summary>
        /// 完全長かどうかを判定し、根拠ごとレポートへ残す
        /// </summary>
        /// <param name="p_結果">アセンブリ実行結果</param>
        /// <param name="p_原入力"></param>
        /// <param name="p_最終パス">最終成果物パス</param>
        /// <param name="p_ポリッシュ統計">ポリッシュ統計、実行していなければ null</param>
        /// <param name="p_閉鎖検証">環状閉鎖の検証結果、実行していなければ null</param>
        /// <param name="p_支持検査">リード支持の検査結果</param>
        /// <param name="p_出力ディレクトリ">レポートの出力先ディレクトリ</param>
        private static void V_出力_完全性レポート(アセンブリ実行結果 p_結果, Parameters p_原入力, string p_最終パス, ポリッシュ統計? p_ポリッシュ統計, IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証, 支持検査結果? p_支持検査, string p_出力ディレクトリ)
        {
            var l_曖昧箇所 = AmbiguityRecorder.Get_記録(p_結果.A_k長);
            var l_未解決ギャップ数 = CompletenessValidator.Get_未解決ギャップ数(p_最終パス);
            var l_環状本数 = CompletenessValidator.Get_環状本数(p_最終パス);
            var l_判定 = CompletenessValidator.Get_判定結果(l_未解決ギャップ数, p_結果.A_整合性検査, p_閉鎖検証, p_ポリッシュ統計, l_曖昧箇所, p_支持検査);
            CompletenessValidator.V_出力_判定結果(l_判定);
            var l_レポートパス = Path.Combine(p_出力ディレクトリ, Consts.レポートファイル名);
            var l_曖昧箇所パス = Path.Combine(p_出力ディレクトリ, Consts.曖昧箇所ファイル名);
            var l_配列群 = FastaReader.Get_全エントリ(p_最終パス).Select(x => x.A_配列).ToList();
            const int l_統計の最小長 = 500;
            ReportWriter.V_書き出し_レポート(l_レポートパス, p_結果.A_k長, AssemblyStatsReporter.Get_統計(l_配列群), l_未解決ギャップ数, l_環状本数, l_判定, p_結果.A_整合性検査, p_閉鎖検証, p_ポリッシュ統計, l_曖昧箇所, AssemblyStatsReporter.Get_N分割統計(l_配列群, l_統計の最小長), l_統計の最小長, p_原入力.A_コピー数基準の出所.ToString(), p_結果.A_実際のコピー数基準.ToString(), p_原入力.A_Is低カバレッジ端トリミング, AssemblyStatsReporter.Get_統計(l_配列群.Where(x => x.Length >= l_統計の最小長)), p_結果.A_固定アンカー評価, PhaseTimingRecorder.Get_記録());
            Logger.V_出力(メッセージID.レポートを書き出した, l_レポートパス);
            var l_Markdownパス = Path.Combine(p_出力ディレクトリ, C_統計表ファイル名);
            ReportWriter.V_書き出し_Markdownレポート(l_Markdownパス, p_結果.A_k長, AssemblyStatsReporter.Get_統計表(Get_最終成果物群(p_出力ディレクトリ, p_最終パス)), l_判定, l_未解決ギャップ数, l_環状本数, l_曖昧箇所.Count, p_原入力.A_コピー数基準の出所.ToString(), p_結果.A_実際のコピー数基準.ToString(), p_結果.A_整合性検査, p_結果.A_固定アンカー評価, p_ポリッシュ統計, p_支持検査, p_閉鎖検証, PhaseTimingRecorder.Get_記録());
            Logger.V_出力(メッセージID.統計表を書き出した, l_Markdownパス);
            ReportWriter.V_書き出し_曖昧箇所(l_曖昧箇所パス, l_曖昧箇所);
            Logger.V_出力(メッセージID.曖昧箇所を書き出した, l_曖昧箇所.Count, l_曖昧箇所パス);
            AmbiguityRecorder.V_保存_履歴(p_出力ディレクトリ);
            if (p_支持検査 is { } l_支持検査)
            {
                var l_支持パス = Path.Combine(p_出力ディレクトリ, C_支持のない箇所ファイル名);
                ReportWriter.V_書き出し_未支持箇所(l_支持パス, l_支持検査.A_区間, l_支持検査.A_r長);
                Logger.V_出力(メッセージID.支持のない箇所を書き出した, l_支持検査.A_区間.Count, l_支持パス);
            }
        }

        /// <summary>
        /// FASTA の ID 行から配列名 (最初の空白まで) を取り出す
        /// </summary>
        /// <param name="p_ID">ID 行</param>
        /// <returns>配列名</returns>
        private static string Get_配列名(string p_ID)
        {
            return p_ID.TrimStart('>').Split(' ', 2)[0];
        }

        #endregion
    }
}
