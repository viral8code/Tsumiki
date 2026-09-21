using System.Text;
using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// コマンドラインの実行とパイプラインを調停する
    /// </summary>
    internal static class AssemblyApplication
    {
        #region 公開メソッド

        /// <summary>
        /// エントリポイント
        /// </summary>
        /// <param name="p_引数列">コマンドライン引数</param>
        /// <returns>成功なら 0、失敗なら 1</returns>
        public static int Get_終了コード(string[] p_引数列)
        {
            var l_元設定 = ConfigurationManager.A_実行時引数;
            try
            {
                V_実行(p_引数列);
                return 0;
            }
            catch (Exception l_例外)
            {
                Logger.V_出力_エラー("Unhandled Tsumiki's method", l_例外);
                return 1;
            }
            finally
            {
                Logger.V_終了_ファイル出力();
                ConfigurationManager.A_実行時引数 = l_元設定;
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 引数を解釈してアセンブリを最後まで走らせる
        /// </summary>
        /// <param name="p_引数列">コマンドライン引数</param>
        private static void V_実行(string[] p_引数列)
        {
            // 日本語・中国語の文言をそのまま出せるようにする
            Console.OutputEncoding = Encoding.UTF8;

            if (p_引数列.Length == 0)
            {
                Console.WriteLine(HelpText.Get_概要());
                return;
            }

            var l_引数 = ArgumentsReader.Get_実行時引数(p_引数列);
            ConfigurationManager.A_実行時引数 = l_引数;
            Messages.A_言語 = l_引数.A_言語;
            Logger.A_水準 = l_引数.A_ログ水準;

            if (l_引数.A_Isバージョンモード)
            {
                Console.WriteLine(HelpText.Get_概要());
                return;
            }

            if (l_引数.A_Isヘルプモード)
            {
                Console.WriteLine(HelpText.Get_ヘルプ());
                return;
            }

            // データを見なければ決まらないパラメータは、パラメータ一覧を
            // 表示する前に確定させる
            // 表示された値と実際に使う値が食い違うと、
            // 後からログを読んだときに何が起きたのか分からなくなる
            PhredSniffer.V_解決_Phredオフセット(l_引数, l_引数.A_リード1のパス, l_引数.A_リード2のパス);

            // 長さの違うライブラリが混ざっていても高い k へ届くよう、
            // 中央値ではなく「その長さ以上のリードが塩基の 1 割以上を出している」最長の長さで k を決める
            // 長さが 1 種類なら従来と同じ値になる
            var l_リード長 = ReadLengthSniffer.Get_梯子上限のリード長(l_引数.A_リード1のパス, l_引数.A_リード2のパス, out var l_リード長分布);
            if (l_リード長 is { } l_観測リード長)
            {
                Logger.V_出力(メッセージID.リード長の観測値, l_観測リード長);
            }
            if (l_リード長分布.Count > 1)
            {
                Logger.V_出力_そのまま(FormattableString.Invariant(
                    $"[Info] リード長の分布: {string.Join(", ", l_リード長分布.OrderByDescending(x => x.Key).Select(x => $"{x.Key}bp x {x.Value:N0}"))}"));
            }
            KmerLengthSelector.V_解決_k長(l_引数, l_リード長);

            Logger.V_出力_そのまま(l_引数.ToString());

            Logger.V_出力_タイムスタンプ();

            var l_一時ディレクトリ = Path.Combine(Environment.CurrentDirectory, l_引数.A_一時ディレクトリ);

            var l_絶対作業パス = Path.GetFullPath(l_一時ディレクトリ).TrimEnd(Path.DirectorySeparatorChar) + Path.DirectorySeparatorChar;
            foreach (var l_入力 in new[] { l_引数.A_リード1のパス, l_引数.A_リード2のパス })
            {
                if (!string.IsNullOrWhiteSpace(l_入力) && Path.GetFullPath(l_入力).StartsWith(l_絶対作業パス, StringComparison.OrdinalIgnoreCase))
                {
                    throw new ArgumentException("Input reads must be outside the output directory");
                }
            }

            if (Path.Exists(l_一時ディレクトリ))
            {
                // 中身を上書きすると、前回の成果と今回の成果が混ざった状態になる
                // 再開を明示されたときだけ、残っているものを使うことを許す
                if (!l_引数.A_Is再開)
                {
                    Logger.V_出力(メッセージID.一時ディレクトリが既にある, l_引数.A_一時ディレクトリ);
                    Logger.V_出力(メッセージID.パスの確認);
                    throw new IOException("Output directory already exists; choose a new -t directory or specify -rs");
                }
                Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_引数.A_一時ディレクトリ);
            }

            _ = Directory.CreateDirectory(l_一時ディレクトリ);

            // ここまでに出た行 (Phred の推定・パラメータ一覧など) も
            // 控えから書き出される
            Logger.V_開始_ファイル出力(l_一時ディレクトリ);

            var l_原入力 = l_引数.Get_複製();
            ReadPreparationPipeline.V_実行(l_引数, l_一時ディレクトリ);

            // アンカー評価が推定k長を書き換える前に、前処理・エラー訂正後のパスを持つ
            // 時点の設定を控えておく (出所記録で使う)
            var l_処理済み設定 = l_引数.Get_複製();

            // -k に複数指定するのは「これらを試して選べ」という意味なので、
            // -mk を別途書かせない
            アセンブリ実行結果 l_結果;
            if (l_引数.A_Isマルチk || l_引数.A_k長一覧.Count > 1)
            {
                l_結果 = MultiKAssembler.Get_実行結果(l_引数, l_一時ディレクトリ, l_リード長, l_原入力) ?? throw new InvalidOperationException("Assembly could not produce a result");
            }
            else
            {
                // マルチk比較が無くても、採用した k やコピー数基準に依らない固定の物差しでの
                // 評価をレポートへ残す (以前は単一 k 指定時にこの評価が一切残らなかった)
                var l_単一結果 = AssemblyPipeline.Get_実行結果(l_引数, l_引数.A_k長, l_一時ディレクトリ, l_リード長, p_原入力: l_原入力) ?? throw new InvalidOperationException("Assembly could not produce a result");
                l_結果 = MultiKAssembler.Get_固定アンカー評価を付与(l_単一結果, l_引数, l_一時ディレクトリ, l_リード長);
            }
            FinalAssemblyPipeline.V_実行(l_結果, l_原入力, l_一時ディレクトリ, l_リード長, l_処理済み設定);

            if (l_引数.A_Is一時ディレクトリ削除)
            {
                AssemblyWorkspace.V_削除_中間ファイル(l_一時ディレクトリ);
            }

            Logger.V_出力(メッセージID.開発中);

            Logger.V_出力_タイムスタンプ();
        }

        #endregion
    }
}
