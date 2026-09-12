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
        #region 公開メソッド

        /// <summary>
        /// 最終成果物から、リード長より短い配列を落とす
        /// </summary>
        /// <param name="p_パス">対象の FASTA パス</param>
        /// <param name="p_リード長">下限として使うリード長</param>
        /// <remarks>
        /// リード 1 本に収まる長さの配列は、リードそのものが既に持っている以上の情報を運ばない<br/>
        /// 加えてこの帯にはタンデムリピートのコピー数を誤って繋いだ断片が集まりやすく、下流の注釈ツールも同種の閾値で捨てる<br/>
        /// 落とした分は一時ディレクトリの k ごとの成果物にそのまま残る
        /// </remarks>
        internal static void V_除外_短い配列(string p_パス, int? p_リード長)
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
        public static void V_実行(アセンブリ実行結果 p_結果, Parameters p_原入力, string p_一時ディレクトリ, int? p_リード長)
        {
            // 採用した 1 組だけを作業ディレクトリの直下へ出す (k ごとの成果物は
            // k のサブディレクトリに残る)
            // この後の除外・ポリッシュは
            // assembly.fasta にだけ効かせ、各段階の出力はそのまま残す
            var l_最終パス = AssemblyPipeline.V_複製_最終成果物(p_結果, p_一時ディレクトリ);

            V_除外_短い配列(l_最終パス, p_リード長);

            var l_ポリッシュ統計 = V_磨く(p_原入力, p_一時ディレクトリ, l_最終パス);
            var l_閉鎖検証 = V_検証_環状閉鎖(p_原入力, l_最終パス);
            var l_支持検査 = V_検査_リードの支持(p_原入力, l_最終パス);

            p_結果 = p_結果 with { A_整合性検査 = Get_最終整合性(p_原入力, p_結果.A_k長, l_最終パス, p_一時ディレクトリ) };
            V_記録_出所(p_原入力, l_最終パス, p_一時ディレクトリ);
            V_出力_完全性レポート(p_結果, l_最終パス, l_ポリッシュ統計, l_閉鎖検証, l_支持検査, p_一時ディレクトリ);

            Logger.V_出力(メッセージID.最終成果物, l_最終パス);

        }

        #endregion

        #region 内部メソッド

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
            l_設定.A_曖昧塩基を許容するか = false;
            try
            {
                ConfigurationManager.A_実行時引数 = l_設定;
                var l_検査パス = Path.Combine(p_作業パス, "validation");
                Directory.CreateDirectory(l_検査パス);
                using var l_索引 = new TrustedKmerIndex(l_検査パス);
                KmerCounting.V_読込_リードペア(l_設定, l_索引, p_進行状況を出力するか: false);
                var l_分布 = l_索引.Get_出現回数ヒストグラム();
                var l_基準 = KmerSpectrumMixtureModel.Get_解析結果(l_分布)?.A_単一コピー平均
                    ?? KmerHistogram.Get_解析結果(l_分布)?.A_ピーク出現回数 ?? 0D;
                KmerCutoffSelector.V_解決_kmerカットオフ(l_設定, l_索引);
                _ = l_索引.V_カットオフ(l_設定.A_kmerカットオフ);
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
        private static void V_記録_出所(Parameters p_原入力, string p_最終パス, string p_作業パス)
        {
            using var l_出力 = File.Create(Path.Combine(p_作業パス, "assembly.provenance.json"));
            using var l_JSON = new System.Text.Json.Utf8JsonWriter(l_出力, new System.Text.Json.JsonWriterOptions { Indented = true });
            l_JSON.WriteStartObject();
            l_JSON.WriteNumber("schema_version", 1);
            l_JSON.WriteString("build_id", typeof(FinalAssemblyPipeline).Assembly.ManifestModule.ModuleVersionId);
            l_JSON.WriteString("validation_source", "uncorrected reads from the assembly library; not independent holdout data");
            l_JSON.WriteString("assembly_sha256", StageCheckpoint.Get_ハッシュ(p_最終パス));
            l_JSON.WriteString("settings", p_原入力.ToString());
            l_JSON.WriteStartArray("inputs");
            foreach (var l_入力 in new[] { p_原入力.A_リード1のパス, p_原入力.A_リード2のパス })
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
            l_JSON.WriteEndObject();
        }

        /// <summary>
        /// 最終成果物の各位置がリードに裏付けられているかを調べる
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <param name="p_最終パス">検査対象の最終成果物パス</param>
        /// <returns>支持検査の結果</returns>
        /// <remarks>
        /// ポリッシュの後に行う<br/>
        /// ギャップ充填・局所アセンブリ・ポリッシュはどれも組み立て後に配列を書き換えるので、それより前に調べても最後に手が入った箇所を見ないことになる<br/>
        /// 加工前のリードを使うが、組み立てと同じライブラリなので独立した検証データではない
        /// </remarks>
        private static 支持検査結果? V_検査_リードの支持(Parameters p_引数, string p_最終パス)
        {
            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.支持検査の開始, Consts.支持検査のr長);
            var l_結果 = ReadSupportChecker.Get_検査結果(p_最終パス, p_引数.A_リード1のパス, p_引数.A_リード2のパス, Consts.支持検査のr長);
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
        /// <remarks>
        /// -po が無ければ何もしない<br/>
        /// 磨いた結果は同じファイル名へ被せ、利用者が受け取るものを 1 つに保つ
        /// </remarks>
        private static ポリッシュ統計? V_磨く(Parameters p_引数, string p_一時ディレクトリ, string p_最終パス)
        {
            if (!p_引数.A_ポリッシュするか)
            {
                return null;
            }

            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.ポリッシュ開始);
            var l_出力先 = Path.Combine(p_一時ディレクトリ, Consts.ポリッシュ済みファイル名);
            var l_統計 = Polisher.Get_磨いた結果(p_最終パス, p_引数.A_リード1のパス, p_引数.A_リード2のパス, l_出力先);
            Polisher.V_出力_統計(l_統計);
            if (l_統計 is not null)
            {
                File.Copy(l_出力先, p_最終パス, overwrite: true);
                var l_最終深度 = Polisher.Get_磨いた結果(p_最終パス, p_引数.A_リード1のパス, p_引数.A_リード2のパス, l_出力先, p_訂正するか: false);
                l_統計 = l_最終深度 is { } l_測定 ? l_測定 with { A_訂正した塩基数 = l_統計.Value.A_訂正した塩基数 } : null;
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
        /// <remarks>
        /// -cc が無ければ何もしない<br/>
        /// 検証していないことと、検証して支持が無かったことは別なので、前者は null を返して判定不能として扱わせる
        /// </remarks>
        private static IReadOnlyList<環状閉鎖検証結果>? V_検証_環状閉鎖(Parameters p_引数, string p_最終パス)
        {
            if (!p_引数.A_環状閉鎖を検証するか)
            {
                return null;
            }

            Logger.V_出力_空行();
            var l_検証 = CircularClosureVerifier.Get_検証結果(p_最終パス, p_引数.A_リード1のパス, p_引数.A_リード2のパス);
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
        private static void V_出力_完全性レポート(アセンブリ実行結果 p_結果, string p_最終パス, ポリッシュ統計? p_ポリッシュ統計, IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証, 支持検査結果? p_支持検査, string p_出力ディレクトリ)
        {
            var l_曖昧箇所 = AmbiguityRecorder.Get_記録(p_結果.A_k長);
            var l_未解決ギャップ数 = CompletenessValidator.Get_未解決ギャップ数(p_最終パス);
            var l_環状本数 = CompletenessValidator.Get_環状本数(p_最終パス);

            var l_判定 = CompletenessValidator.Get_判定結果(l_未解決ギャップ数, p_結果.A_整合性検査, p_閉鎖検証, p_ポリッシュ統計, l_曖昧箇所, p_支持検査);
            CompletenessValidator.V_出力_判定結果(l_判定);

            var l_レポートパス = Path.Combine(p_出力ディレクトリ, Consts.レポートファイル名);
            var l_曖昧箇所パス = Path.Combine(p_出力ディレクトリ, Consts.曖昧箇所ファイル名);

            ReportWriter.V_書き出し_レポート(l_レポートパス, p_結果.A_k長, AssemblyStatsReporter.Get_統計_FASTA(p_最終パス), l_未解決ギャップ数, l_環状本数, l_判定, p_結果.A_整合性検査, p_閉鎖検証, p_ポリッシュ統計, l_曖昧箇所);
            Logger.V_出力(メッセージID.レポートを書き出した, l_レポートパス);

            ReportWriter.V_書き出し_曖昧箇所(l_曖昧箇所パス, l_曖昧箇所);
            Logger.V_出力(メッセージID.曖昧箇所を書き出した, l_曖昧箇所.Count, l_曖昧箇所パス);

            if (p_支持検査 is { } l_支持検査)
            {
                var l_支持パス = Path.Combine(p_出力ディレクトリ, Consts.支持のない箇所ファイル名);
                ReportWriter.V_書き出し_支持のない箇所(l_支持パス, l_支持検査.A_区間, l_支持検査.A_r長);
                Logger.V_出力(メッセージID.支持のない箇所を書き出した, l_支持検査.A_区間.Count, l_支持パス);
            }
        }

        #endregion
    }
}
