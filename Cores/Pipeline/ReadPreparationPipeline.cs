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
            // 訂正済みリードがそのまま使えるなら前処理まで遡らない
            // 前処理済みリードは訂正の入力にしか使わないので、消してあっても再開できる
            if (p_引数.A_Is再開 && p_引数.A_Isエラー訂正 && Get_再利用できる訂正済み(p_引数, p_一時ディレクトリ) is { } l_再利用)
            {
                Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_再利用.A_訂正済み1);
                p_引数.A_リード1のパス = l_再利用.A_訂正済み1;
                if (l_再利用.A_訂正済み2 is not null)
                {
                    p_引数.A_リード2のパス = l_再利用.A_訂正済み2;
                }
                V_削除_前処理済みリード(p_一時ディレクトリ);
                return;
            }

            if (p_引数.A_Is前処理)
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
                    if (p_引数.A_Is再開 && StageCheckpoint.Is再利用可能(l_署名, l_前処理済み1, l_前処理済み2))
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

            if (p_引数.A_Isエラー訂正)
            {
                Logger.V_出力(メッセージID.エラー訂正開始);

                var l_訂正済み1 = Path.Combine(p_一時ディレクトリ, "corrected.1.fq");
                var l_Hasリード2 = !string.IsNullOrWhiteSpace(p_引数.A_リード2のパス);
                var l_訂正済み2 = l_Hasリード2 ? Path.Combine(p_一時ディレクトリ, "corrected.2.fq") : null;

                var l_署名 = StageCheckpoint.Get_入力署名(p_引数);
                if (p_引数.A_Is再開 && StageCheckpoint.Is再利用可能(l_署名, l_訂正済み1, l_訂正済み2))
                {
                    Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_訂正済み1);
                }
                else
                {
                    ErrorCorrector.V_訂正_リードファイル(p_引数.A_リード1のパス, l_Hasリード2 ? p_引数.A_リード2のパス : null, p_一時ディレクトリ, l_訂正済み1, l_訂正済み2);
                    StageCheckpoint.V_保存(l_署名, l_訂正済み1, l_訂正済み2);
                }

                // 以降の全処理 (k-mer カウント・グラフ構築・リードの再マッピング) は
                // 訂正済みファイルを見るようにする
                p_引数.A_リード1のパス = l_訂正済み1;
                if (l_Hasリード2)
                {
                    p_引数.A_リード2のパス = l_訂正済み2!;
                }

                V_削除_前処理済みリード(p_一時ディレクトリ);

                Logger.V_出力_タイムスタンプ();
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 訂正済みリードが揃った後、要らなくなった前処理済みリードを消す
        /// </summary>
        /// <param name="p_一時ディレクトリ">処理済みリードの置き場</param>
        /// <remarks>
        /// 前処理済みリードは訂正の入力にしか使わない (以降の工程は訂正済み、最終検査は元のリードを読む)<br/>
        /// 再開の照合は工程の記録 (.sha256) だけで足りるので、記録は残して実体だけ消す
        /// </remarks>
        internal static void V_削除_前処理済みリード(string p_一時ディレクトリ)
        {
            foreach (var l_名前 in new[] { "preprocessed.1.fq", "preprocessed.2.fq" })
            {
                var l_パス = Path.Combine(p_一時ディレクトリ, l_名前);
                if (File.Exists(l_パス))
                {
                    File.Delete(l_パス);
                    Logger.V_出力(メッセージID.中間リードを削除, l_パス);
                }
            }
        }

        /// <summary>
        /// 前処理を省いて使える訂正済みリードがあるか調べる
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
        /// <param name="p_一時ディレクトリ">処理済みリードの置き場</param>
        /// <remarks>
        /// 訂正工程の入力は前処理済みリードだが、その識別には前処理工程の記録 (.sha256) を使うので実体は要らない
        /// </remarks>
        /// <returns>そのまま使える訂正済みリード、無ければ null</returns>
        private static (string A_訂正済み1, string? A_訂正済み2)? Get_再利用できる訂正済み(Parameters p_引数, string p_一時ディレクトリ)
        {
            var l_訂正済み1 = Path.Combine(p_一時ディレクトリ, "corrected.1.fq");
            var l_Hasリード2 = !string.IsNullOrWhiteSpace(p_引数.A_リード2のパス);
            var l_訂正済み2 = l_Hasリード2 ? Path.Combine(p_一時ディレクトリ, "corrected.2.fq") : null;

            var l_訂正前設定 = p_引数.Get_複製();
            if (p_引数.A_Is前処理 && l_Hasリード2)
            {
                var l_前処理済み1 = Path.Combine(p_一時ディレクトリ, "preprocessed.1.fq");
                var l_前処理済み2 = Path.Combine(p_一時ディレクトリ, "preprocessed.2.fq");
                if (StageCheckpoint.Get_保存済み記録(l_前処理済み1) is null && !(File.Exists(l_前処理済み1) && File.Exists(l_前処理済み2)))
                {
                    return null;
                }
                l_訂正前設定.A_リード1のパス = l_前処理済み1;
                l_訂正前設定.A_リード2のパス = l_前処理済み2;
            }

            return StageCheckpoint.Is再利用可能(StageCheckpoint.Get_入力署名(l_訂正前設定), l_訂正済み1, l_訂正済み2) ? (l_訂正済み1, l_訂正済み2) : null;
        }

        #endregion
    }
}
