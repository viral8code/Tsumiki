using System.Runtime.CompilerServices;
using Tsumiki.Model.Foundation;

namespace Tsumiki.Common
{
    /// <summary>
    /// 進行状況と結果を、画面とファイルへ書き出す
    /// </summary>
    /// <remarks>
    /// 画面へ出す量は水準で絞れるが、ファイルへは常に全量を残す
    /// </remarks>
    internal class Logger
    {
        /// <summary>
        /// 錠
        /// </summary>
        protected internal static readonly object _錠 = new();

        /// <summary>
        /// 画面へ出す量
        /// </summary>
        /// <remarks>
        /// ファイルへの記録はこれに関わらず全量を残す
        /// </remarks>
        public static ログ水準 A_水準 { get; set; } = ログ水準.標準;

        /// <summary>
        /// 全量を残す書き出し先、まだ開いていなければ null
        /// </summary>
        private static StreamWriter? _ファイル;

        /// <summary>
        /// 0 より大きい間は何も出さない
        /// </summary>
        /// <remarks>
        /// 入れ子にできるよう数で持つ
        /// </remarks>
        protected internal static int _休止の深さ;

        /// <summary>
        /// 一時ディレクトリを作る前に出た行の控え
        /// </summary>
        /// <remarks>
        /// Phred の推定やパラメータ一覧はディレクトリの用意より前に出るため、
        /// そのままでは記録から漏れる
        /// </remarks>
        private static readonly List<string> _書き出し待ち = [];

        /// <summary>
        /// 控えの上限
        /// </summary>
        /// <remarks>
        /// 通常は一時ディレクトリを作るまでの数十行しか溜まらないが、
        /// ファイルを開かないまま使われ続けても際限なく積まないようにする
        /// </remarks>
        private const int 控えの上限 = 10_000;

        /// <summary>
        /// 呼び出し元のメソッド名を返す
        /// </summary>
        /// <param name="p_メソッド名">呼び出し元のメソッド名、コンパイラが埋める</param>
        /// <returns>呼び出し元のメソッド名</returns>
        public static string Get_メソッド名([CallerMemberName] string p_メソッド名 = "")
        {
            return p_メソッド名;
        }

        /// <summary>
        /// 以降の出力をファイルにも残す
        /// </summary>
        /// <remarks>
        /// 既にあれば追記する
        /// (再開したときに前回までの経過が消えないようにする)
        /// </remarks>
        public static void V_開始_ファイル出力(string p_一時ディレクトリ)
        {
            lock (_錠)
            {
                if (_ファイル is not null)
                {
                    return;
                }
                _ファイル = new StreamWriter(
                    Path.Combine(p_一時ディレクトリ, Consts.ログファイル名), append: true)
                {
                    // 長時間走るので、途中で落ちても直前までが残るようにする
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
        /// <remarks>
        /// 文言は言語ごとのカタログから引く
        /// </remarks>
        public static void V_出力(メッセージID p_ID, params object?[] p_引数)
        {
            V_書き出し(Messages.Get_文言(p_ID, p_引数), p_標準エラーか: false);
        }

        /// <summary>
        /// 標準エラーへ 1 行出す
        /// </summary>
        public static void V_出力_標準エラー(メッセージID p_ID, params object?[] p_引数)
        {
            V_書き出し(Messages.Get_文言(p_ID, p_引数), p_標準エラーか: true);
        }

        /// <summary>
        /// カタログを通さない文字列をそのまま出す
        /// </summary>
        /// <remarks>
        /// パラメータ一覧のように、訳す対象ではないが記録には残したいもの向け
        /// </remarks>
        public static void V_出力_そのまま(string p_文)
        {
            V_書き出し(p_文, p_標準エラーか: false);
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

            // 例外の内容そのものは訳す対象ではない
            V_書き出し(p_例外.ToString(), p_標準エラーか: true);
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
            V_書き出し(p_例外.ToString(), p_標準エラーか: true);
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
        /// <remarks>
        /// 文言を持たないのでカタログには載せない
        /// </remarks>
        public static void V_出力_空行()
        {
            V_書き出し(string.Empty, p_標準エラーか: false);
        }

        /// <summary>
        /// 1 行を、必要ならば画面へ出し、常にファイルへ残す
        /// </summary>
        /// <remarks>
        /// 標準エラーへ出すもの (警告・エラー) は水準によらず必ず画面にも出す
        /// </remarks>
        private static void V_書き出し(string p_行, bool p_標準エラーか)
        {
            lock (_錠)
            {
                if (_休止の深さ > 0)
                {
                    return;
                }
                if (_ファイル is { } l_ファイル)
                {
                    l_ファイル.WriteLine(p_行);
                }
                else if (_書き出し待ち.Count < 控えの上限)
                {
                    _書き出し待ち.Add(p_行);
                }

                if (p_標準エラーか)
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
        /// <remarks>
        /// 行頭の目印で決まる<br/>
        /// 目印を持たない行 (進行状況の見出しなど) は標準扱いとする
        /// </remarks>
        private static ログ水準 Get_水準(string p_行)
        {
            return p_行.StartsWith(Consts.ログ目印.詳細, StringComparison.Ordinal)
                ? ログ水準.詳細
                : p_行.StartsWith(Consts.ログ目印.完全性, StringComparison.Ordinal)
                || p_行.StartsWith(Consts.ログ目印.レポート, StringComparison.Ordinal)
                ? ログ水準.最小
                : ログ水準.標準;
        }

        /// <summary>
        /// この場を抜けるまで、画面にもファイルにも何も出さない
        /// </summary>
        /// <remarks>
        /// 局所アセンブリのように、小さな使い捨ての処理を数百回繰り返す
        /// 区間で使う<br/>
        /// 1 回あたりの索引の統計は、集めても読む意味が無い割に
        /// 本来のログを埋め尽くす (実データでは k=21 だけで千行を超えた)
        /// </remarks>
        public static IDisposable V_止める_記録()
        {
            return new 記録の休止();
        }

        /// <summary>
        /// 記録を閉じる
        /// </summary>
        /// <remarks>
        /// ここまでに書いたものは失われない
        /// </remarks>
        public static void V_終了_ファイル出力()
        {
            lock (_錠)
            {
                _ファイル?.Dispose();
                _ファイル = null;
            }
        }
    }

    /// <summary>
    /// Logger の記録を一時的に止めるための解放用ハンドル
    /// </summary>
    internal sealed class 記録の休止 : IDisposable
    {
        public 記録の休止()
        {
            lock (Logger._錠)
            {
                Logger._休止の深さ++;
            }
        }

        /// <summary>
        /// 保持している資源を解放する
        /// </summary>
        public void Dispose()
        {
            lock (Logger._錠)
            {
                Logger._休止の深さ--;
            }
        }
    }
}
