using Tsumiki.Commons;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Models.Correction;
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
        /// ファイル名 fq
        /// </summary>
        private const string C_ファイル名_fq = ".fq";

        /// <summary>
        /// 項目 preprocess
        /// </summary>
        private const string C_項目_preprocess = "preprocess";

        /// <summary>
        /// 前処理済みリードのファイル名の幹
        /// </summary>
        private const string C_前処理済みの幹 = "preprocessed";

        /// <summary>
        /// 訂正済みリードのファイル名の幹
        /// </summary>
        private const string C_訂正済みの幹 = "corrected";

        /// <summary>
        /// 訂正済みリードの横に置く、誤りの無い区間の度数のファイルの接尾辞
        /// </summary>
        private const string C_度数の接尾辞 = ".segments.tsv";

        /// <summary>
        /// 項目 thinning
        /// </summary>
        private const string C_項目_thinning = "thinning";

        /// <summary>
        /// 目標の深度へ間引いたリードのファイル名の幹
        /// </summary>
        private const string C_間引き済みの幹 = "thinned";

        /// <summary>
        /// 間引いたときの元のライブラリ群と採用割合を記録するファイル名 (一時ディレクトリに置く、1 行の採用割合)
        /// </summary>
        private const string C_間引きの情報ファイル名 = "thinned.info";

        /// <summary>
        /// 圧縮して置くリードのファイル名の接尾辞 (gzip)
        /// </summary>
        private const string C_圧縮ファイル拡張子 = ".gz";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 作業用の入力を前処理し、処理済みパスへ切り替える
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
        /// <param name="p_一時ディレクトリ">処理済みリードの出力先</param>
        public static void V_実行(Parameters p_引数, string p_一時ディレクトリ)
        {
            p_引数.A_無誤り区間の度数群.Clear();
            if (p_引数.A_Is再開 && p_引数.A_Isエラー訂正 && Get_再利用できる訂正済み(p_引数, p_一時ディレクトリ) is { } l_再利用)
            {
                Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_再利用[0].A_順リード);
                l_再利用.ForEach(x => V_読込_度数(p_引数, x.A_順リード));
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

        /// <summary>
        /// 入力の推定深度が目標の深度を十分に超えていれば、リードのペアを一様に間引き、入力を間引いたリードへ切り替える
        /// </summary>
        /// <param name="p_引数">入力の設定</param>
        /// <param name="p_一時ディレクトリ">間引いたリードの出力先</param>
        public static void V_間引き(Parameters p_引数, string p_一時ディレクトリ)
        {
            if (p_引数.A_目標深度 == 0)
            {
                return;
            }

            using var l_計測 = new StageTimer(C_項目_thinning);
            var l_署名 = Get_署名(p_引数);
            var l_圧縮 = !中間データ置き場.A_Is有効;
            List<(string A_出力1, string? A_出力2)> l_出力群 = [];
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var l_Hasペア = !string.IsNullOrWhiteSpace(p_引数.A_ライブラリ群[i].A_逆リード);
                l_出力群.Add((Get_間引きパス(p_一時ディレクトリ, i, 1, l_圧縮), l_Hasペア ? Get_間引きパス(p_一時ディレクトリ, i, 2, l_圧縮) : null));
            }

            if (p_引数.A_Is再開 && l_署名 is not null && l_出力群.Count > 0 && l_出力群.All(x => StageCheckpoint.Is再利用可能(l_署名, x.A_出力1, x.A_出力2)))
            {
                Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_出力群[0].A_出力1);
                if (ReadThinner.Get_読込_採用割合(Path.Combine(p_一時ディレクトリ, C_間引きの情報ファイル名)) is { } l_再利用の採用割合)
                {
                    p_引数.A_間引き前ライブラリ群 = [.. p_引数.A_ライブラリ群];
                    p_引数.A_間引きの採用割合 = l_再利用の採用割合;
                }

                p_引数.Set_ライブラリ群(l_出力群.Select(x => (x.A_出力1, x.A_出力2 ?? string.Empty)));
                return;
            }

            var l_総塩基数 = ReadThinner.Get_総塩基数(p_引数);
            var l_推定ゲノムサイズ = ReadThinner.Get_推定ゲノムサイズ(p_引数, p_一時ディレクトリ);
            if (l_推定ゲノムサイズ is not > 0)
            {
                Logger.V_出力(メッセージID.間引き_推定できず);
                return;
            }

            var l_深度 = (double)l_総塩基数 / l_推定ゲノムサイズ.Value;
            var l_割合 = ReadThinner.Get_採用割合(l_総塩基数, l_推定ゲノムサイズ.Value, p_引数.A_目標深度);
            if (l_割合 is not { } l_採用割合)
            {
                Logger.V_出力(メッセージID.間引き_不要, l_深度, p_引数.A_目標深度, Consts.間引きを始める倍率);
                return;
            }

            long l_元の件数 = 0;
            long l_採用件数 = 0;
            List<(string A_順リード, string A_逆リード)> l_間引き後群 = [];
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var (l_順リード, l_逆リード) = p_引数.A_ライブラリ群[i];
                var (l_出力1, l_出力2) = l_出力群[i];
                var (l_元, l_採用) = ReadThinner.V_間引き_ファイル(l_順リード, l_逆リード, l_出力1, l_出力2, l_採用割合, p_Is圧縮: l_圧縮);
                l_元の件数 += l_元;
                l_採用件数 += l_採用;
                V_保存_記録(l_署名, l_出力1, l_出力2);
                l_間引き後群.Add((l_出力1, l_出力2 ?? string.Empty));
            }

            Logger.V_出力(メッセージID.間引き_結果, l_深度, l_総塩基数, l_推定ゲノムサイズ.Value, p_引数.A_目標深度, l_元の件数, l_採用件数);
            ReadThinner.V_書き出し_採用割合(Path.Combine(p_一時ディレクトリ, C_間引きの情報ファイル名), l_採用割合);
            p_引数.A_間引き前ライブラリ群 = [.. p_引数.A_ライブラリ群];
            p_引数.A_間引きの採用割合 = l_採用割合;
            p_引数.Set_ライブラリ群(l_間引き後群);
            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// 訂正済みリードが揃った後、要らなくなった前処理済みリードを消す
        /// </summary>
        /// <param name="p_一時ディレクトリ">処理済みリードの置き場</param>
        public static void V_削除_前処理済みリード(string p_一時ディレクトリ)
        {
            foreach (var l_パス in 中間データ置き場.Get_一覧(p_一時ディレクトリ, C_前処理済みの幹, C_ファイル名_fq))
            {
                中間データ置き場.V_削除(l_パス);
                Logger.V_出力(メッセージID.中間リードを削除, l_パス);
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
            if (!p_引数.A_Hasペア)
            {
                Logger.V_出力(メッセージID.前処理省略_ペアなし);
                return;
            }

            Logger.V_出力(メッセージID.前処理開始);
            using var l_計測 = new StageTimer(C_項目_preprocess);
            var l_ペアなしを飛ばした = false;
            var l_署名 = Get_署名(p_引数);
            List<(string A_順リード, string A_逆リード)> l_出力群 = [];
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var (l_順リード, l_逆リード) = p_引数.A_ライブラリ群[i];
                if (string.IsNullOrWhiteSpace(l_逆リード))
                {
                    l_出力群.Add((l_順リード, string.Empty));
                    l_ペアなしを飛ばした = true;
                    continue;
                }

                var l_出力1 = Get_中間パス(p_一時ディレクトリ, C_前処理済みの幹, i, 1);
                var l_出力2 = Get_中間パス(p_一時ディレクトリ, C_前処理済みの幹, i, 2);
                if (p_引数.A_Is再開 && StageCheckpoint.Is再利用可能(l_署名!, l_出力1, l_出力2))
                {
                    Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_出力1);
                }
                else
                {
                    var l_前処理統計 = Preprocessor.V_前処理_リードファイル(l_順リード, l_逆リード, l_出力1, l_出力2, p_引数.Get_Phredオフセット(i));
                    Preprocessor.V_出力_前処理統計(l_前処理統計);
                    V_保存_記録(l_署名, l_出力1, l_出力2);
                }

                l_出力群.Add((l_出力1, l_出力2));
            }

            if (l_ペアなしを飛ばした)
            {
                Logger.V_出力(メッセージID.前処理省略_ペアなし);
            }

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
            List<(string A_順リード, string A_逆リード)> l_出力群 = [];
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var (l_順リード, l_逆リード) = p_引数.A_ライブラリ群[i];
                var l_Has逆リード = !string.IsNullOrWhiteSpace(l_逆リード);
                var l_出力1 = Get_中間パス(p_一時ディレクトリ, C_訂正済みの幹, i, 1);
                var l_出力2 = l_Has逆リード ? Get_中間パス(p_一時ディレクトリ, C_訂正済みの幹, i, 2) : null;
                if (p_引数.A_Is再開 && StageCheckpoint.Is再利用可能(l_署名!, l_出力1, l_出力2))
                {
                    Logger.V_出力(メッセージID.再開_中間ファイルを再利用, l_出力1);
                    V_読込_度数(p_引数, l_出力1);
                }
                else
                {
                    var l_度数 = ErrorCorrector.V_訂正_リードファイル(l_順リード, l_Has逆リード ? l_逆リード : null, p_一時ディレクトリ, l_出力1, l_出力2, p_引数.Get_Phredオフセット(i));
                    if (l_度数 is not null)
                    {
                        p_引数.A_無誤り区間の度数群.Add(l_度数);
                        if (!中間データ置き場.A_Is有効)
                        {
                            l_度数.V_書き出し(l_出力1 + C_度数の接尾辞);
                        }
                    }

                    V_保存_記録(l_署名, l_出力1, l_出力2);
                }

                l_出力群.Add((l_出力1, l_出力2 ?? string.Empty));
            }

            p_引数.Set_ライブラリ群(l_出力群);
            V_削除_前処理済みリード(p_一時ディレクトリ);
            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// 訂正済みリードの横に残した誤りの無い区間の度数を読み込む
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
        /// <param name="p_訂正済み1">訂正済みのリード 1 のパス</param>
        private static void V_読込_度数(Parameters p_引数, string p_訂正済み1)
        {
            if (無誤り区間の度数.Get_読込(p_訂正済み1 + C_度数の接尾辞) is { } l_度数)
            {
                p_引数.A_無誤り区間の度数群.Add(l_度数);
            }
        }

        /// <summary>
        /// 中間リードのパス
        /// </summary>
        /// <param name="p_一時ディレクトリ">置き場</param>
        /// <param name="p_幹">ファイル名の幹</param>
        /// <param name="p_ライブラリ番号">0 起点のライブラリ番号</param>
        /// <param name="p_side">1 か 2</param>
        /// <returns></returns>
        private static string Get_中間パス(string p_一時ディレクトリ, string p_幹, int p_ライブラリ番号, int p_side)
        {
            var l_接尾 = p_ライブラリ番号 == 0 ? string.Empty : FormattableString.Invariant($".lib{p_ライブラリ番号 + 1}");
            return Path.Combine(p_一時ディレクトリ, FormattableString.Invariant($"{p_幹}{l_接尾}.{p_side}.fq"));
        }

        /// <summary>
        /// 間引いたリードのパス (圧縮するときは gzip の接尾辞を付ける)。判定と書き出しの両方がこれを使う
        /// </summary>
        /// <param name="p_一時ディレクトリ">置き場</param>
        /// <param name="p_ライブラリ番号">0 起点のライブラリ番号</param>
        /// <param name="p_side">1 か 2</param>
        /// <param name="p_Is圧縮">gzip で置くか</param>
        /// <returns>間引いたリードのパス</returns>
        private static string Get_間引きパス(string p_一時ディレクトリ, int p_ライブラリ番号, int p_side, bool p_Is圧縮)
        {
            var l_パス = Get_中間パス(p_一時ディレクトリ, C_間引き済みの幹, p_ライブラリ番号, p_side);
            return p_Is圧縮 ? l_パス + C_圧縮ファイル拡張子 : l_パス;
        }

        /// <summary>
        /// 再開の照合に使う署名
        /// </summary>
        /// <param name="p_引数">作業用設定</param>
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
        /// <returns>そのまま使える訂正済みリード、1 つでも欠けていれば null</returns>
        private static List<(string A_順リード, string A_逆リード)>? Get_再利用できる訂正済み(Parameters p_引数, string p_一時ディレクトリ)
        {
            var l_訂正前設定 = p_引数.Get_複製();
            if (p_引数.A_Is前処理 && p_引数.A_Hasペア)
            {
                List<(string A_順リード, string A_逆リード)> l_前処理済み群 = [];
                for (var i = 0; i < p_引数.A_ライブラリ数; i++)
                {
                    var (l_順リード, l_逆リード) = p_引数.A_ライブラリ群[i];
                    if (string.IsNullOrWhiteSpace(l_逆リード))
                    {
                        l_前処理済み群.Add((l_順リード, string.Empty));
                        continue;
                    }

                    var l_前処理済み1 = Get_中間パス(p_一時ディレクトリ, C_前処理済みの幹, i, 1);
                    var l_前処理済み2 = Get_中間パス(p_一時ディレクトリ, C_前処理済みの幹, i, 2);
                    if (StageCheckpoint.Get_保存済み記録(l_前処理済み1) is null && !(File.Exists(l_前処理済み1) && File.Exists(l_前処理済み2)))
                    {
                        return null;
                    }

                    l_前処理済み群.Add((l_前処理済み1, l_前処理済み2));
                }

                l_訂正前設定.Set_ライブラリ群(l_前処理済み群);
            }

            var l_署名 = StageCheckpoint.Get_入力署名(l_訂正前設定);
            List<(string A_順リード, string A_逆リード)> l_訂正済み群 = [];
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var l_Has逆リード = !string.IsNullOrWhiteSpace(p_引数.A_ライブラリ群[i].A_逆リード);
                var l_訂正済み1 = Get_中間パス(p_一時ディレクトリ, C_訂正済みの幹, i, 1);
                var l_訂正済み2 = l_Has逆リード ? Get_中間パス(p_一時ディレクトリ, C_訂正済みの幹, i, 2) : null;
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
