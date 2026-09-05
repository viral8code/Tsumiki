using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
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
                var l_一時ディレクトリ = Path.Combine(
                    Environment.CurrentDirectory, ConfigurationManager.A_実行時引数.A_一時ディレクトリ);
                if (Directory.Exists(l_一時ディレクトリ))
                {
                    Directory.Delete(l_一時ディレクトリ, true);
                }
            }
        }

        private static void V_実行(string[] p_引数列)
        {
            if (p_引数列.Length == 0)
            {
                Console.WriteLine(Consts.概要テキスト);
                Environment.Exit(0);
            }

            var l_引数 = ArgumentsReader.Get_実行時引数(p_引数列);
            ConfigurationManager.A_実行時引数 = l_引数;

            if (l_引数.A_ヘルプモードか)
            {
                Console.WriteLine(Consts.ヘルプテキスト);
                Environment.Exit(0);
            }

            // データを見なければ決まらないパラメータは、パラメータ一覧を
            // 表示する前に確定させる。表示された値と実際に使う値が食い違うと、
            // 後からログを読んだときに何が起きたのか分からなくなる。
            PhredSniffer.V_解決_Phredオフセット(l_引数, l_引数.A_リード1のパス, l_引数.A_リード2のパス);

            var l_リード長 = ReadLengthSniffer.Get_代表リード長(l_引数.A_リード1のパス, l_引数.A_リード2のパス);
            if (l_リード長 is { } l_観測リード長)
            {
                Console.WriteLine($"[Info] Read length (median of sampled reads): {l_観測リード長} bp");
            }
            KmerLengthSelector.V_解決_k長(l_引数, l_リード長);

            Console.WriteLine(l_引数);

            Logger.V_出力_タイムスタンプ();

            var l_一時ディレクトリ = Path.Combine(Environment.CurrentDirectory, l_引数.A_一時ディレクトリ);

            if (Path.Exists(l_一時ディレクトリ))
            {
                Console.WriteLine($"{l_引数.A_一時ディレクトリ} already exists");
                Console.WriteLine("Please check path!");
                Environment.Exit(0);
            }

            _ = Directory.CreateDirectory(l_一時ディレクトリ);

            if (l_引数.A_エラー訂正するか)
            {
                Console.WriteLine("Correcting reads before assembly");

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

            var l_結果 = l_引数.A_マルチkか
                ? MultiKAssembler.Get_実行結果(l_引数, l_一時ディレクトリ, l_リード長)
                : AssemblyPipeline.Get_実行結果(
                    l_引数, l_引数.A_k長, l_一時ディレクトリ, p_出力接頭辞: string.Empty, l_リード長);

            if (l_結果 is null)
            {
                Console.WriteLine("""

                    This genome is too complex to assembly...
                    Please adjust the parameters!

                    """);
                return;
            }

            Console.WriteLine("開発中！");

            Logger.V_出力_タイムスタンプ();
        }
    }
}
