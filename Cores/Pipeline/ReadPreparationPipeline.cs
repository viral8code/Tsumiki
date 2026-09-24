using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// 前処理とエラー訂正を管理する
    /// </summary>
    internal static class ReadPreparationPipeline
    {
        #region 定数

        /// <summary>
        /// 前処理済みリードのファイル名の幹
        /// </summary>
        private const string 前処理済みの幹 = "preprocessed";

        /// <summary>
        /// 訂正済みリードのファイル名の幹
        /// </summary>
        private const string 訂正済みの幹 = "corrected";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 作業用の入力を前処理し、処理済みパスへ切り替える
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
        /// <param name="p_一時ディレクトリ">処理済みリードの出力先</param>
        /// <remarks>
        /// 前処理も訂正もライブラリごとに独立して行う (アダプタの読み抜けも誤りの出方もライブラリで違う)
        /// </remarks>
        public static void V_実行(Parameters p_引数, string p_一時ディレクトリ)
        {
            // 訂正済みリードがそのまま使えるなら前処理まで遡らない
            // 前処理済みリードは訂正の入力にしか使わないので、消してあっても再開できる
            if (p_引数.A_Is再開 && p_引数.A_Isエラー訂正 && Get_再利用できる訂正済み(p_引数, p_一時ディレクトリ) is { } l_再利用)
            {
                Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_再利用[0].A_リード1);
                p_引数.Set_ライブラリ群(l_再利用);
                V_削除_前処理済みリード(p_一時ディレクトリ);
                return;
            }

            if (p_引数.A_Is前処理)
            {
                V_前処理(p_引数, p_一時ディレクトリ);
            }

            if (p_引数.A_Isエラー訂正)
            {
                V_訂正(p_引数, p_一時ディレクトリ);
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// ライブラリごとにアダプタ除去とペア相互訂正を行う
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
        /// <param name="p_一時ディレクトリ">出力先</param>
        private static void V_前処理(Parameters p_引数, string p_一時ディレクトリ)
        {
            if (!p_引数.Hasペア)
            {
                Logger.V_出力(メッセージID.前処理省略_ペアなし);
                return;
            }

            Logger.V_出力(メッセージID.前処理開始);
            using var l_計測 = new StageTimer("preprocess");
            var l_ペアなしを飛ばした = false;

            // 署名は入力を差し替える前に採る (差し替えた後では前処理の入力を指さなくなる)
            var l_署名 = Get_署名(p_引数);
            List<(string A_リード1, string A_リード2)> l_出力群 = [];
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var (A_リード1, A_リード2) = p_引数.A_ライブラリ群[i];

                // アダプタの読み抜けもペアの相互訂正も相方が要る。シングルエンドのライブラリはそのまま通す
                if (string.IsNullOrWhiteSpace(A_リード2))
                {
                    l_出力群.Add((A_リード1, string.Empty));
                    l_ペアなしを飛ばした = true;
                    continue;
                }

                var l_出力1 = Get_中間パス(p_一時ディレクトリ, 前処理済みの幹, i, 1);
                var l_出力2 = Get_中間パス(p_一時ディレクトリ, 前処理済みの幹, i, 2);

                if (p_引数.A_Is再開 && StageCheckpoint.Is再利用可能(l_署名!, l_出力1, l_出力2))
                {
                    Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_出力1);
                }
                else
                {
                    var l_前処理統計 = Preprocessor.V_前処理_リードファイル(A_リード1, A_リード2, l_出力1, l_出力2, p_引数.Get_Phredオフセット(i));
                    Preprocessor.V_出力_前処理統計(l_前処理統計);
                    V_保存_記録(l_署名, l_出力1, l_出力2);
                }
                l_出力群.Add((l_出力1, l_出力2));
            }

            if (l_ペアなしを飛ばした)
            {
                Logger.V_出力(メッセージID.前処理省略_ペアなし);
            }

            // 以降の全処理 (エラー訂正・ k-mer カウント・グラフ構築) は
            // 前処理済みファイルを見るようにする
            p_引数.Set_ライブラリ群(l_出力群);

            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// ライブラリごとにエラー訂正を行う
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
        /// <param name="p_一時ディレクトリ">出力先</param>
        private static void V_訂正(Parameters p_引数, string p_一時ディレクトリ)
        {
            Logger.V_出力(メッセージID.エラー訂正開始);

            var l_署名 = Get_署名(p_引数);
            List<(string A_リード1, string A_リード2)> l_出力群 = [];
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var (A_リード1, A_リード2) = p_引数.A_ライブラリ群[i];
                var l_Hasリード2 = !string.IsNullOrWhiteSpace(A_リード2);
                var l_出力1 = Get_中間パス(p_一時ディレクトリ, 訂正済みの幹, i, 1);
                var l_出力2 = l_Hasリード2 ? Get_中間パス(p_一時ディレクトリ, 訂正済みの幹, i, 2) : null;

                if (p_引数.A_Is再開 && StageCheckpoint.Is再利用可能(l_署名!, l_出力1, l_出力2))
                {
                    Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_出力1);
                }
                else
                {
                    ErrorCorrector.V_訂正_リードファイル(A_リード1, l_Hasリード2 ? A_リード2 : null, p_一時ディレクトリ, l_出力1, l_出力2, p_引数.Get_Phredオフセット(i));
                    V_保存_記録(l_署名, l_出力1, l_出力2);
                }
                l_出力群.Add((l_出力1, l_出力2 ?? string.Empty));
            }

            // 以降の全処理 (k-mer カウント・グラフ構築・リードの再マッピング) は
            // 訂正済みファイルを見るようにする
            p_引数.Set_ライブラリ群(l_出力群);

            V_削除_前処理済みリード(p_一時ディレクトリ);

            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// 中間リードのパス
        /// </summary>
        /// <param name="p_一時ディレクトリ">置き場</param>
        /// <param name="p_幹">ファイル名の幹</param>
        /// <param name="p_ライブラリ番号">0 起点のライブラリ番号</param>
        /// <param name="p_side">1 か 2</param>
        /// <remarks>
        /// 先頭のライブラリだけ従来の名前にするのは、単一ライブラリの再開が過去の中間ファイルで効くようにするため
        /// </remarks>
        /// <returns></returns>
        private static string Get_中間パス(string p_一時ディレクトリ, string p_幹, int p_ライブラリ番号, int p_side)
        {
            var l_接尾 = p_ライブラリ番号 == 0 ? string.Empty : FormattableString.Invariant($".lib{p_ライブラリ番号 + 1}");
            return Path.Combine(p_一時ディレクトリ, FormattableString.Invariant($"{p_幹}{l_接尾}.{p_side}.fq"));
        }

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
            foreach (var l_パス in 中間データ置き場.Get_一覧(p_一時ディレクトリ, 前処理済みの幹, ".fq"))
            {
                中間データ置き場.V_削除(l_パス);
                Logger.V_出力(メッセージID.中間リードを削除, l_パス);
            }
        }

        /// <summary>
        /// 再開の照合に使う署名
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
        /// <remarks>
        /// 署名は入力リードを丸ごと読んでハッシュを取るので重い。中間データをメモリに置くときは再開の元が残らず使い道が無いので取らない
        /// </remarks>
        /// <returns>取らないときは null</returns>
        private static string? Get_署名(Parameters p_引数)
        {
            return 中間データ置き場.A_Is有効 ? null : StageCheckpoint.Get_入力署名(p_引数);
        }

        /// <summary>
        /// 工程の記録 (.sha256) を残す
        /// </summary>
        /// <param name="p_署名">入力の署名、取っていなければ null</param>
        /// <param name="p_出力">出力</param>
        /// <param name="p_対出力">対になる出力</param>
        private static void V_保存_記録(string? p_署名, string p_出力, string? p_対出力)
        {
            if (p_署名 is not null)
            {
                StageCheckpoint.V_保存(p_署名, p_出力, p_対出力);
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
        /// <returns>そのまま使える訂正済みリード、1 つでも欠けていれば null</returns>
        private static List<(string A_リード1, string A_リード2)>? Get_再利用できる訂正済み(Parameters p_引数, string p_一時ディレクトリ)
        {
            var l_訂正前設定 = p_引数.Get_複製();
            if (p_引数.A_Is前処理 && p_引数.Hasペア)
            {
                List<(string A_リード1, string A_リード2)> l_前処理済み群 = [];
                for (var i = 0; i < p_引数.A_ライブラリ数; i++)
                {
                    var (A_リード1, A_リード2) = p_引数.A_ライブラリ群[i];
                    if (string.IsNullOrWhiteSpace(A_リード2))
                    {
                        l_前処理済み群.Add((A_リード1, string.Empty));
                        continue;
                    }

                    var l_前処理済み1 = Get_中間パス(p_一時ディレクトリ, 前処理済みの幹, i, 1);
                    var l_前処理済み2 = Get_中間パス(p_一時ディレクトリ, 前処理済みの幹, i, 2);
                    if (StageCheckpoint.Get_保存済み記録(l_前処理済み1) is null && !(File.Exists(l_前処理済み1) && File.Exists(l_前処理済み2)))
                    {
                        return null;
                    }
                    l_前処理済み群.Add((l_前処理済み1, l_前処理済み2));
                }
                l_訂正前設定.Set_ライブラリ群(l_前処理済み群);
            }

            var l_署名 = StageCheckpoint.Get_入力署名(l_訂正前設定);
            List<(string A_リード1, string A_リード2)> l_訂正済み群 = [];
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var l_Hasリード2 = !string.IsNullOrWhiteSpace(p_引数.A_ライブラリ群[i].A_リード2);
                var l_訂正済み1 = Get_中間パス(p_一時ディレクトリ, 訂正済みの幹, i, 1);
                var l_訂正済み2 = l_Hasリード2 ? Get_中間パス(p_一時ディレクトリ, 訂正済みの幹, i, 2) : null;
                if (!StageCheckpoint.Is再利用可能(l_署名, l_訂正済み1, l_訂正済み2))
                {
                    return null;
                }
                l_訂正済み群.Add((l_訂正済み1, l_訂正済み2 ?? string.Empty));
            }
            return l_訂正済み群;
        }

        #endregion
    }
}
