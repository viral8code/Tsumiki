using System.Text;
using Tsumiki.Common;
using Tsumiki.Core.Evaluation;
using Tsumiki.Core.Output;
using Tsumiki.Core.Pipeline;
using Tsumiki.Core.Polishing;
using Tsumiki.Core.Preprocessing;
using Tsumiki.IO;
using Tsumiki.Model.Evaluation;
using Tsumiki.Model.Foundation;
using Tsumiki.Model.Polishing;
using Tsumiki.Utility;

namespace Tsumiki
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            try
            {
                V_実行(args);
            }
            catch (Exception ex)
            {
                Logger.V_出力_エラー("Unhandled Tsumiki's method", ex);
            }
            finally
            {
                // 一時ディレクトリには k ごとの成果物が入っており、後から
                // 見比べたくなることが多い。消すかどうかは利用者に決めさせる。
                var l_引数 = ConfigurationManager.A_実行時引数;
                var l_一時ディレクトリ = Path.Combine(
                    Environment.CurrentDirectory, l_引数.A_一時ディレクトリ);
                if (Directory.Exists(l_一時ディレクトリ))
                {
                    if (l_引数.A_一時ディレクトリを削除するか)
                    {
                        V_削除_中間ファイル(l_一時ディレクトリ);
                        Logger.V_出力(メッセージID.中間ファイルを削除した, l_引数.A_一時ディレクトリ);
                    }
                    else
                    {
                        Logger.V_出力(メッセージID.一時ディレクトリを残した, l_引数.A_一時ディレクトリ, Consts.引数キー.一時ディレクトリ削除);
                    }
                    Logger.V_出力(
                        メッセージID.ログの保存先,
                        Path.Combine(l_引数.A_一時ディレクトリ, Consts.ログファイル名));
                }
                Logger.V_終了_ファイル出力();
            }
        }

        /// <summary>
        /// 作業ディレクトリから中間ファイルだけを消す。最終成果物とログは残す。
        ///
        /// 作業ディレクトリは利用者が受け取る成果物の置き場でもあるので、
        /// 消してよいのは k ごとの途中経過や訂正済みリードのほうだけになる。
        /// ログを残すのは、消す指定をした実行こそ後から確かめる手段が
        /// それしかなくなるため。
        /// </summary>
        internal static void V_削除_中間ファイル(string p_作業ディレクトリ)
        {
            foreach (var l_ディレクトリ in Directory.EnumerateDirectories(p_作業ディレクトリ))
            {
                Directory.Delete(l_ディレクトリ, recursive: true);
            }
            foreach (var l_ファイル in Directory.EnumerateFiles(p_作業ディレクトリ))
            {
                if (!Consts.最終成果物のファイル名.Contains(Path.GetFileName(l_ファイル)))
                {
                    File.Delete(l_ファイル);
                }
            }
        }

        private static void V_実行(string[] p_引数列)
        {
            // 日本語・中国語の文言をそのまま出せるようにする。
            Console.OutputEncoding = Encoding.UTF8;

            if (p_引数列.Length == 0)
            {
                Console.WriteLine(HelpText.Get_概要());
                Environment.Exit(0);
            }

            var l_引数 = ArgumentsReader.Get_実行時引数(p_引数列);
            ConfigurationManager.A_実行時引数 = l_引数;
            Messages.A_言語 = l_引数.A_言語;
            Logger.A_水準 = l_引数.A_ログ水準;

            if (l_引数.A_バージョンモードか)
            {
                Console.WriteLine(HelpText.Get_概要());
                Environment.Exit(0);
            }

            if (l_引数.A_ヘルプモードか)
            {
                Console.WriteLine(HelpText.Get_ヘルプ());
                Environment.Exit(0);
            }

            // データを見なければ決まらないパラメータは、パラメータ一覧を
            // 表示する前に確定させる。表示された値と実際に使う値が食い違うと、
            // 後からログを読んだときに何が起きたのか分からなくなる。
            PhredSniffer.V_解決_Phredオフセット(l_引数, l_引数.A_リード1のパス, l_引数.A_リード2のパス);

            var l_リード長 = ReadLengthSniffer.Get_代表リード長(l_引数.A_リード1のパス, l_引数.A_リード2のパス);
            if (l_リード長 is { } l_観測リード長)
            {
                Logger.V_出力(メッセージID.リード長の観測値, l_観測リード長);
            }
            KmerLengthSelector.V_解決_k長(l_引数, l_リード長);

            Logger.V_出力_そのまま(l_引数.ToString());

            Logger.V_出力_タイムスタンプ();

            var l_一時ディレクトリ = Path.Combine(Environment.CurrentDirectory, l_引数.A_一時ディレクトリ);

            if (Path.Exists(l_一時ディレクトリ))
            {
                // 中身を上書きすると、前回の成果と今回の成果が混ざった状態になる。
                // 再開を明示されたときだけ、残っているものを使うことを許す。
                if (!l_引数.A_再開するか)
                {
                    Logger.V_出力(メッセージID.一時ディレクトリが既にある, l_引数.A_一時ディレクトリ);
                    Logger.V_出力(メッセージID.パスの確認);
                    Environment.Exit(0);
                }
                Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_引数.A_一時ディレクトリ);
            }

            _ = Directory.CreateDirectory(l_一時ディレクトリ);

            // ここまでに出た行(Phred の推定・パラメータ一覧など)も
            // 控えから書き出される。
            Logger.V_開始_ファイル出力(l_一時ディレクトリ);

            if (l_引数.A_前処理するか)
            {
                if (string.IsNullOrWhiteSpace(l_引数.A_リード2のパス))
                {
                    Logger.V_出力(メッセージID.前処理省略_ペアなし);
                }
                else
                {
                    Logger.V_出力(メッセージID.前処理開始);

                    var l_前処理済み1 = Path.Combine(l_一時ディレクトリ, "preprocessed.1.fq");
                    var l_前処理済み2 = Path.Combine(l_一時ディレクトリ, "preprocessed.2.fq");

                    if (Get_再利用できるか(l_引数, l_前処理済み1, l_前処理済み2))
                    {
                        Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_前処理済み1);
                    }
                    else
                    {
                        var l_前処理統計 = Preprocessor.V_前処理_リードファイル(
                            l_引数.A_リード1のパス, l_引数.A_リード2のパス, l_前処理済み1, l_前処理済み2);
                        Preprocessor.V_出力_前処理統計(l_前処理統計);
                    }

                    // 以降の全処理(エラー訂正・k-merカウント・グラフ構築)は
                    // 前処理済みファイルを見るようにする。
                    l_引数.A_リード1のパス = l_前処理済み1;
                    l_引数.A_リード2のパス = l_前処理済み2;

                    Logger.V_出力_タイムスタンプ();
                }
            }

            if (l_引数.A_エラー訂正するか)
            {
                Logger.V_出力(メッセージID.エラー訂正開始);

                var l_訂正済み1 = Path.Combine(l_一時ディレクトリ, "corrected.1.fq");
                var l_リード2があるか = !string.IsNullOrWhiteSpace(l_引数.A_リード2のパス);
                var l_訂正済み2 = l_リード2があるか ? Path.Combine(l_一時ディレクトリ, "corrected.2.fq") : null;

                if (Get_再利用できるか(l_引数, l_訂正済み1, l_訂正済み2))
                {
                    Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_訂正済み1);
                }
                else
                {
                    ErrorCorrector.V_訂正_リードファイル(
                        l_引数.A_リード1のパス,
                        l_リード2があるか ? l_引数.A_リード2のパス : null,
                        l_一時ディレクトリ,
                        l_訂正済み1,
                        l_訂正済み2);
                }

                // 以降の全処理(k-merカウント・グラフ構築・リードの再マッピング)は
                // 訂正済みファイルを見るようにする。
                l_引数.A_リード1のパス = l_訂正済み1;
                if (l_リード2があるか)
                {
                    l_引数.A_リード2のパス = l_訂正済み2!;
                }

                Logger.V_出力_タイムスタンプ();
            }

            // -k に複数指定するのは「これらを試して選べ」という意味なので、
            // -mk を別途書かせない。
            var l_結果 = l_引数.A_マルチkか || l_引数.A_k長一覧.Count > 1
                ? MultiKAssembler.Get_実行結果(l_引数, l_一時ディレクトリ, l_リード長)
                : AssemblyPipeline.Get_実行結果(
                    l_引数, l_引数.A_k長, l_一時ディレクトリ, l_リード長);

            if (l_結果 is null)
            {
                Logger.V_出力(メッセージID.アセンブリ不能);
                return;
            }

            // 採用した1組だけを作業ディレクトリの直下へ出す(k ごとの成果物は
            // k のサブディレクトリに残る)。
            AssemblyPipeline.V_複製_最終成果物(l_結果, l_一時ディレクトリ);

            var l_最終パス = Path.Combine(
                l_一時ディレクトリ,
                l_結果.A_スキャフォールドパス is null
                    ? Consts.コンティグファイル名
                    : Consts.スキャフォールドファイル名);

            V_除外_短い配列(l_最終パス, l_リード長);

            var l_ポリッシュ統計 = V_磨く(l_引数, l_一時ディレクトリ, l_最終パス);
            var l_閉鎖検証 = V_検証_環状閉鎖(l_引数, l_最終パス);

            V_出力_完全性レポート(l_結果, l_最終パス, l_ポリッシュ統計, l_閉鎖検証, l_一時ディレクトリ);

            Logger.V_出力(メッセージID.開発中);

            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// 最終成果物から、リード長より短い配列を落とす。
        ///
        /// リード1本に収まる長さの配列は、リードそのものが既に持っている以上の
        /// 情報を運ばない。加えてこの帯にはタンデムリピートのコピー数を誤って
        /// 繋いだ断片が集まりやすく、下流の注釈ツールも同種の閾値で捨てる。
        /// 落とした分は一時ディレクトリの k ごとの成果物にそのまま残る。
        /// </summary>
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
        /// 再開が指定されていて、その工程の出力が既に揃っているか。
        /// 揃っていれば作り直さずそのまま使う。
        /// </summary>
        private static bool Get_再利用できるか(Parameters p_引数, string p_出力1, string? p_出力2)
        {
            return p_引数.A_再開するか
                && File.Exists(p_出力1)
                && (p_出力2 is null || File.Exists(p_出力2));
        }

        /// <summary>
        /// 最終成果物にリードを貼り直して磨く。-po が無ければ何もしない。
        /// 磨いた結果は同じファイル名へ被せ、利用者が受け取るものを1つに保つ。
        /// </summary>
        private static ポリッシュ統計? V_磨く(
            Parameters p_引数, string p_一時ディレクトリ, string p_最終パス)
        {
            if (!p_引数.A_ポリッシュするか)
            {
                return null;
            }

            Logger.V_出力_空行();
            Logger.V_出力(メッセージID.ポリッシュ開始);
            var l_出力先 = Path.Combine(p_一時ディレクトリ, Consts.ポリッシュ済みファイル名);
            var l_統計 = Polisher.Get_磨いた結果(
                p_最終パス, p_引数.A_リード1のパス, p_引数.A_リード2のパス, l_出力先);
            Polisher.V_出力_統計(l_統計);
            if (l_統計 is not null)
            {
                File.Copy(l_出力先, p_最終パス, overwrite: true);
            }
            Logger.V_出力_タイムスタンプ();
            return l_統計;
        }

        /// <summary>
        /// 環状の閉じ目を元リードで確かめる。-cc が無ければ何もしない。
        /// 検証していないことと、検証して支持が無かったことは別なので、
        /// 前者は null を返して判定不能として扱わせる。
        /// </summary>
        private static IReadOnlyList<環状閉鎖検証結果>? V_検証_環状閉鎖(
            Parameters p_引数, string p_最終パス)
        {
            if (!p_引数.A_環状閉鎖を検証するか)
            {
                return null;
            }

            Logger.V_出力_空行();
            var l_検証 = CircularClosureVerifier.Get_検証結果(
                p_最終パス, p_引数.A_リード1のパス, p_引数.A_リード2のパス);
            CircularClosureVerifier.V_出力_検証結果(l_検証);
            Logger.V_出力_タイムスタンプ();
            return l_検証;
        }

        /// <summary>
        /// 完全長かどうかを判定し、根拠ごとレポートへ残す。
        /// </summary>
        private static void V_出力_完全性レポート(
            アセンブリ実行結果 p_結果,
            string p_最終パス,
            ポリッシュ統計? p_ポリッシュ統計,
            IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証,
            string p_出力ディレクトリ)
        {
            var l_曖昧箇所 = AmbiguityRecorder.Get_記録(p_結果.A_k長);
            var l_未解決ギャップ数 = CompletenessValidator.Get_未解決ギャップ数(p_最終パス);
            var l_環状本数 = CompletenessValidator.Get_環状本数(p_最終パス);

            var l_判定 = CompletenessValidator.Get_判定結果(
                l_未解決ギャップ数, p_結果.A_整合性検査, p_閉鎖検証, p_ポリッシュ統計, l_曖昧箇所);
            CompletenessValidator.V_出力_判定結果(l_判定);

            var l_レポートパス = Path.Combine(p_出力ディレクトリ, Consts.レポートファイル名);
            var l_曖昧箇所パス = Path.Combine(p_出力ディレクトリ, Consts.曖昧箇所ファイル名);

            ReportWriter.V_書き出し_レポート(
                l_レポートパス,
                p_結果.A_k長,
                AssemblyStatsReporter.Get_統計_FASTA(p_最終パス),
                l_未解決ギャップ数,
                l_環状本数,
                l_判定,
                p_結果.A_整合性検査,
                p_閉鎖検証,
                p_ポリッシュ統計,
                l_曖昧箇所);
            Logger.V_出力(メッセージID.レポートを書き出した, l_レポートパス);

            ReportWriter.V_書き出し_曖昧箇所(l_曖昧箇所パス, l_曖昧箇所);
            Logger.V_出力(メッセージID.曖昧箇所を書き出した, l_曖昧箇所.Count, l_曖昧箇所パス);
        }
    }
}
