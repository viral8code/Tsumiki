using System.Text;
using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Model;
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
                        Directory.Delete(l_一時ディレクトリ, true);
                    }
                    else
                    {
                        Logger.V_出力(メッセージID.一時ディレクトリを残した, l_引数.A_一時ディレクトリ, Consts.引数キー.一時ディレクトリ削除);
                    }
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

            Console.WriteLine(l_引数);

            Logger.V_出力_タイムスタンプ();

            var l_一時ディレクトリ = Path.Combine(Environment.CurrentDirectory, l_引数.A_一時ディレクトリ);

            if (Path.Exists(l_一時ディレクトリ))
            {
                Logger.V_出力(メッセージID.一時ディレクトリが既にある, l_引数.A_一時ディレクトリ);
                Logger.V_出力(メッセージID.パスの確認);
                Environment.Exit(0);
            }

            _ = Directory.CreateDirectory(l_一時ディレクトリ);

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

                    var l_前処理統計 = Preprocessor.V_前処理_リードファイル(
                        l_引数.A_リード1のパス, l_引数.A_リード2のパス, l_前処理済み1, l_前処理済み2);
                    Preprocessor.V_出力_前処理統計(l_前処理統計);

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

                ErrorCorrector.V_訂正_リードファイル(
                    l_引数.A_リード1のパス,
                    l_リード2があるか ? l_引数.A_リード2のパス : null,
                    l_一時ディレクトリ,
                    l_訂正済み1,
                    l_訂正済み2);

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

            // 採用した1組だけを実行ディレクトリへ出す(k ごとの成果物は
            // 一時ディレクトリに残る)。
            AssemblyPipeline.V_複製_最終成果物(l_結果);

            Logger.V_出力(メッセージID.開発中);

            Logger.V_出力_タイムスタンプ();
        }
    }
}
