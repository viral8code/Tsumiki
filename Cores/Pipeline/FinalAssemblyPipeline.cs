using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Output;
using Tsumiki.Cores.Polishing;
using Tsumiki.Cores.Preprocessing;
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

            var l_ポリッシュ統計 = V_磨く(p_原入力, p_一時ディレクトリ, l_最終パス);
            var l_閉鎖検証 = V_検証_環状閉鎖(p_原入力, l_最終パス);
            var l_支持検査 = V_検査_リード支持(p_原入力, l_最終パス);

            p_結果 = p_結果 with { A_整合性検査 = Get_最終整合性(p_原入力, p_結果.A_k長, l_最終パス, p_一時ディレクトリ) };
            V_記録_出所(p_原入力, p_結果, l_最終パス, p_一時ディレクトリ, p_処理済み設定);
            V_出力_最終統計(p_一時ディレクトリ, l_最終パス);
            V_出力_完全性レポート(p_結果, p_原入力, l_最終パス, l_ポリッシュ統計, l_閉鎖検証, l_支持検査, p_一時ディレクトリ);

            Logger.V_出力(メッセージID.最終成果物, l_最終パス);
        }

        #endregion

        #region 内部メソッド

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
            return
            [
                ("contigs", Path.Combine(p_一時ディレクトリ, AssemblyPipeline.C_Contigファイル名)),
                ("scaffolds", Path.Combine(p_一時ディレクトリ, Consts.Scaffoldファイル名)),
                ("assembly", p_最終パス),
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
                var l_検査パス = Path.Combine(p_作業パス, "validation");
                _ = Directory.CreateDirectory(l_検査パス);
                using var l_索引 = new TrustedKmerIndex(l_検査パス);
                KmerCounting.V_読込_リードペア(l_設定, l_索引, p_Is進行状況出力: false);
                var l_分布 = l_索引.Get_出現回数ヒストグラム();
                var l_基準 = KmerSpectrumMixtureModel.Get_解析結果(l_分布)?.A_単一コピー平均
                    ?? KmerHistogram.Get_解析結果(l_分布)?.A_単一コピー基準値 ?? 0D;
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
        /// <param name="p_最終パス">最終配列</param>
        /// <param name="p_作業パス">出力先</param>
        /// <param name="p_処理済み設定">前処理・エラー訂正後のパスを保持した設定、無ければ null</param>
        private static void V_記録_出所(Parameters p_原入力, アセンブリ実行結果 p_結果, string p_最終パス, string p_作業パス, Parameters? p_処理済み設定)
        {
            using var l_出力 = File.Create(Path.Combine(p_作業パス, "assembly.provenance.json"));
            using var l_JSON = new System.Text.Json.Utf8JsonWriter(l_出力, new System.Text.Json.JsonWriterOptions { Indented = true });
            l_JSON.WriteStartObject();
            l_JSON.WriteNumber("schema_version", 1);
            l_JSON.WriteString("build_id", typeof(FinalAssemblyPipeline).Assembly.ManifestModule.ModuleVersionId);
            l_JSON.WriteString("validation_source", "uncorrected reads from the assembly library; not independent holdout data");
            l_JSON.WriteString("assembly_sha256", StageCheckpoint.Get_ハッシュ(p_最終パス));
            l_JSON.WriteNumber("thread_count", p_原入力.A_スレッド数);
            l_JSON.WriteString("settings", p_原入力.ToString());
            l_JSON.WriteStartObject("assembly_settings");
            l_JSON.WriteString("copy_number_baseline_requested", p_原入力.A_コピー数基準の出所.ToString());
            l_JSON.WriteString("copy_number_baseline_actual", p_結果.A_実際のコピー数基準.ToString());
            l_JSON.WriteBoolean("trim_low_coverage_ends", p_原入力.A_Is低カバレッジ端トリミング);
            l_JSON.WriteEndObject();
            l_JSON.WriteStartArray("inputs");
            foreach (var l_入力 in p_原入力.A_ライブラリ群.SelectMany(x => new[] { x.A_リード1, x.A_リード2 }))
            {
                if (string.IsNullOrWhiteSpace(l_入力))
                {
                    continue;
                }

                l_JSON.WriteStartObject();
                l_JSON.WriteString("path", Path.GetFullPath(l_入力));
                l_JSON.WriteString("sha256", StageCheckpoint.Get_ハッシュ(l_入力));
                l_JSON.WriteEndObject();
            }

            l_JSON.WriteEndArray();

            if (p_処理済み設定 is { } l_処理済み設定)
            {
                l_JSON.WriteStartArray("corrected_read_hashes");
                foreach (var l_入力 in l_処理済み設定.A_ライブラリ群.SelectMany(x => new[] { x.A_リード1, x.A_リード2 }))
                {
                    if (!中間データ置き場.Is存在(l_入力))
                    {
                        continue;
                    }

                    l_JSON.WriteStartObject();
                    l_JSON.WriteString("path", Path.GetFullPath(l_入力));
                    l_JSON.WriteString("sha256", StageCheckpoint.Get_ハッシュ(l_入力));
                    l_JSON.WriteEndObject();
                }

                l_JSON.WriteEndArray();

                l_JSON.WriteString("pipeline_fingerprint", StageCheckpoint.Get_入力署名(l_処理済み設定));
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
            ReportWriter.V_書き出し_レポート(l_レポートパス, p_結果.A_k長, AssemblyStatsReporter.Get_統計(l_配列群), l_未解決ギャップ数, l_環状本数, l_判定, p_結果.A_整合性検査, p_閉鎖検証, p_ポリッシュ統計, l_曖昧箇所,
                AssemblyStatsReporter.Get_N分割統計(l_配列群, l_統計の最小長), l_統計の最小長, p_原入力.A_コピー数基準の出所.ToString(), p_結果.A_実際のコピー数基準.ToString(), p_原入力.A_Is低カバレッジ端トリミング,
                AssemblyStatsReporter.Get_統計(l_配列群.Where(x => x.Length >= l_統計の最小長)), p_結果.A_固定アンカー評価, PhaseTimingRecorder.Get_記録());
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

        #endregion
    }
}
