using System.Runtime.CompilerServices;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Commons
{
    /// <summary>
    /// 進行状況と結果を、画面とファイルへ書き出す
    /// </summary>
    internal class Logger
    {
        #region 定数

        /// <summary>
        /// CallerMemberName が置換する既定値
        /// </summary>
        private const string 呼び出し元の既定値 = "";

        /// <summary>
        /// 控えの上限
        /// </summary>
        private const int 控えの上限 = 10_000;

        #endregion

        #region 内部変数

        /// <summary>
        /// 錠
        /// </summary>
        private static readonly Lock _錠 = new();

        /// <summary>
        /// 全量を残す書き出し先、まだ開いていなければ null
        /// </summary>
        private static StreamWriter? _ファイル;

        /// <summary>
        /// 一時ディレクトリを作る前に出た行の控え
        /// </summary>
        private static readonly List<string> _書き出し待ち = [];

        #endregion

        #region プロパティ

        /// <summary>
        /// 画面へ出す量
        /// </summary>
        public static ログ水準 A_水準 { get; set; } = ログ水準.標準;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 呼び出し元のメソッド名を返す
        /// </summary>
        /// <param name="p_メソッド名">呼び出し元のメソッド名、コンパイラが埋める</param>
        /// <returns>呼び出し元のメソッド名</returns>
        public static string Get_メソッド名([CallerMemberName] string p_メソッド名 = 呼び出し元の既定値)
        {
            return p_メソッド名;
        }

        /// <summary>
        /// 以降の出力をファイルにも残す
        /// </summary>
        /// <param name="p_一時ディレクトリ"></param>
        public static void V_開始_ファイル出力(string p_一時ディレクトリ)
        {
            lock (_錠)
            {
                if (_ファイル is not null)
                {
                    return;
                }
                _ファイル = new StreamWriter(Path.Combine(p_一時ディレクトリ, Consts.ログファイル名), append: true)
                {
                    AutoFlush = true,
                };
                foreach (var l_行 in _書き出し待ち)
                {
                    _ファイル.WriteLine(l_行);
                }
                _書き出し待ち.Clear();
            }
        }

        /// <summary>
        /// 標準出力へ 1 行出す
        /// </summary>
        /// <param name="p_ID"></param>
        /// <param name="p_引数"></param>
        public static void V_出力(メッセージID p_ID, params object?[] p_引数)
        {
            V_書き出し(Messages.Get_文言(p_ID, p_引数), p_Is標準エラー: false);
        }

        /// <summary>
        /// 標準エラーへ 1 行出す
        /// </summary>
        /// <param name="p_ID"></param>
        /// <param name="p_引数"></param>
        public static void V_出力_標準エラー(メッセージID p_ID, params object?[] p_引数)
        {
            V_書き出し(Messages.Get_文言(p_ID, p_引数), p_Is標準エラー: true);
        }

        /// <summary>
        /// カタログを通さない文字列をそのまま出す
        /// </summary>
        /// <param name="p_文"></param>
        public static void V_出力_そのまま(string p_文)
        {
            V_書き出し(p_文, p_Is標準エラー: false);
        }

        /// <summary>
        /// 警告を出力する
        /// </summary>
        /// <param name="p_メソッド名">警告を出した場所</param>
        /// <param name="p_例外">警告の元になった例外</param>
        public static void V_出力_警告(string p_メソッド名, Exception p_例外)
        {
            V_出力_標準エラー(メッセージID.例外を無視_見出し);
            V_出力_標準エラー(メッセージID.例外を無視_メソッド, p_メソッド名);

            V_書き出し(p_例外.ToString(), p_Is標準エラー: true);
        }

        /// <summary>
        /// エラーを出力する
        /// </summary>
        /// <param name="p_メソッド名">エラーを出した場所</param>
        /// <param name="p_例外">エラーの元になった例外</param>
        public static void V_出力_エラー(string p_メソッド名, Exception p_例外)
        {
            V_出力_標準エラー(メッセージID.停止_見出し);
            V_出力_標準エラー(メッセージID.停止_メソッド, p_メソッド名);
            V_書き出し(p_例外.ToString(), p_Is標準エラー: true);
        }

        /// <summary>
        /// 現在時刻を出力する
        /// </summary>
        public static void V_出力_タイムスタンプ()
        {
            V_出力(メッセージID.タイムスタンプ, DateTime.Now);
        }

        /// <summary>
        /// 区切りの空行
        /// </summary>
        public static void V_出力_空行()
        {
            V_書き出し(string.Empty, p_Is標準エラー: false);
        }

        /// <summary>
        /// 記録を閉じる
        /// </summary>
        public static void V_終了_ファイル出力()
        {
            lock (_錠)
            {
                _ファイル?.Dispose();
                _ファイル = null;
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 行を、必要ならば画面へ出し、常にファイルへ残す
        /// </summary>
        /// <param name="p_行"></param>
        /// <param name="p_Is標準エラー"></param>
        private static void V_書き出し(string p_行, bool p_Is標準エラー)
        {
            lock (_錠)
            {
                if (_ファイル is { } l_ファイル)
                {
                    l_ファイル.WriteLine(p_行);
                }
                else if (_書き出し待ち.Count < 控えの上限)
                {
                    _書き出し待ち.Add(p_行);
                }

                if (p_Is標準エラー)
                {
                    Console.Error.WriteLine(p_行);
                    return;
                }

                if (Get_水準(p_行) <= A_水準)
                {
                    Console.WriteLine(p_行);
                }
            }
        }

        /// <summary>
        /// その行を出すのに必要な水準
        /// </summary>
        /// <param name="p_行"></param>
        /// <returns></returns>
        private static ログ水準 Get_水準(string p_行)
        {
            return p_行.StartsWith(Consts.ログ目印.詳細, StringComparison.Ordinal)
                ? ログ水準.詳細
                : p_行.StartsWith(Consts.ログ目印.完全性, StringComparison.Ordinal)
                || p_行.StartsWith(Consts.ログ目印.レポート, StringComparison.Ordinal)
                ? ログ水準.最小
                : ログ水準.標準;
        }

        #endregion
    }
}
