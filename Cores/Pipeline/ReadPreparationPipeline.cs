using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// 前処理とエラー訂正を管理する
    /// </summary>
    internal static class ReadPreparationPipeline
    {
        #region 公開メソッド

        /// <summary>
        /// 作業用の入力を前処理し、処理済みパスへ切り替える
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
        /// <param name="p_一時ディレクトリ">処理済みリードの出力先</param>
        public static void V_実行(Parameters p_引数, string p_一時ディレクトリ)
        {
            if (p_引数.A_前処理するか)
            {
                if (string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
                {
                    Logger.V_出力(メッセージID.前処理省略_ペアなし);
                }
                else
                {
                    Logger.V_出力(メッセージID.前処理開始);

                    var l_前処理済み1 = Path.Combine(p_一時ディレクトリ, "preprocessed.1.fq");
                    var l_前処理済み2 = Path.Combine(p_一時ディレクトリ, "preprocessed.2.fq");

                    var l_署名 = StageCheckpoint.Get_入力署名(p_引数);
                    if (p_引数.A_再開するか && StageCheckpoint.Get_再利用可能(l_署名, l_前処理済み1, l_前処理済み2))
                    {
                        Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_前処理済み1);
                    }
                    else
                    {
                        var l_前処理統計 = Preprocessor.V_前処理_リードファイル(p_引数.A_リード1のパス, p_引数.A_リード2のパス, l_前処理済み1, l_前処理済み2);
                        Preprocessor.V_出力_前処理統計(l_前処理統計);
                        StageCheckpoint.V_保存(l_署名, l_前処理済み1, l_前処理済み2);
                    }

                    // 以降の全処理 (エラー訂正・ k-mer カウント・グラフ構築) は
                    // 前処理済みファイルを見るようにする
                    p_引数.A_リード1のパス = l_前処理済み1;
                    p_引数.A_リード2のパス = l_前処理済み2;

                    Logger.V_出力_タイムスタンプ();
                }
            }

            if (p_引数.A_エラー訂正するか)
            {
                Logger.V_出力(メッセージID.エラー訂正開始);

                var l_訂正済み1 = Path.Combine(p_一時ディレクトリ, "corrected.1.fq");
                var l_リード2があるか = !string.IsNullOrWhiteSpace(p_引数.A_リード2のパス);
                var l_訂正済み2 = l_リード2があるか ? Path.Combine(p_一時ディレクトリ, "corrected.2.fq") : null;

                var l_署名 = StageCheckpoint.Get_入力署名(p_引数);
                if (p_引数.A_再開するか && StageCheckpoint.Get_再利用可能(l_署名, l_訂正済み1, l_訂正済み2))
                {
                    Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_訂正済み1);
                }
                else
                {
                    ErrorCorrector.V_訂正_リードファイル(p_引数.A_リード1のパス, l_リード2があるか ? p_引数.A_リード2のパス : null, p_一時ディレクトリ, l_訂正済み1, l_訂正済み2);
                    StageCheckpoint.V_保存(l_署名, l_訂正済み1, l_訂正済み2);
                }

                // 以降の全処理 (k-mer カウント・グラフ構築・リードの再マッピング) は
                // 訂正済みファイルを見るようにする
                p_引数.A_リード1のパス = l_訂正済み1;
                if (l_リード2があるか)
                {
                    p_引数.A_リード2のパス = l_訂正済み2!;
                }

                Logger.V_出力_タイムスタンプ();
            }
        }

        #endregion

    }
}
