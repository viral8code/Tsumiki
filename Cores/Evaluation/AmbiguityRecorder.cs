using System.Globalization;
using Tsumiki.Cores.Evidence;
using Tsumiki.Models.Reporting;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 決めきれずに打ち切った箇所を、その場で書き留めておくための収集器
    /// </summary>
    internal static class AmbiguityRecorder
    {
        #region 定数

        /// <summary>
        /// k ごとの作業ディレクトリに残す控え
        /// </summary>
        /// <remarks>
        /// 再開時に読み直す
        /// </remarks>
        private const string 保存ファイル名 = "ambiguous.tsv";

        /// <summary>
        /// 場所名の既定接頭辞
        /// </summary>
        private const string 既定接頭辞 = "unitig";

        #endregion

        #region 内部変数

        /// <summary>
        /// 錠
        /// </summary>
        private static readonly Lock _錠 = new();

        /// <summary>
        /// k 長ごとの記録
        /// </summary>
        private static readonly Dictionary<int, List<曖昧箇所>> _k長ごとの記録 = [];

        /// <summary>
        /// いま記録している k の長さ
        /// </summary>
        private static int _現在のk長;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// この k の記録を集め直す
        /// </summary>
        /// <param name="p_k長"></param>
        /// <remarks>
        /// 同じ k を再実行した場合は上書きする
        /// </remarks>
        public static void V_開始(int p_k長)
        {
            lock (_錠)
            {
                _現在のk長 = p_k長;
                _k長ごとの記録[p_k長] = [];
            }
        }

        /// <summary>
        /// 1 箇所ぶん書き留める
        /// </summary>
        /// <param name="p_種別"></param>
        /// <param name="p_場所"></param>
        /// <param name="p_首位の支持"></param>
        /// <param name="p_次点の支持"></param>
        /// <param name="p_首位の生支持数"></param>
        /// <remarks>
        /// 確信度は独立な支持本数を飽和関数に通した値で、同じ種類の証拠がいくら積み上がっても 1 に近づくだけになる
        /// </remarks>
        public static void V_記録(曖昧箇所の種別 p_種別, string p_場所, double p_首位の支持 = 0D, double p_次点の支持 = 0D, long p_首位の生支持数 = 0L)
        {
            lock (_錠)
            {
                if (!_k長ごとの記録.TryGetValue(_現在のk長, out var l_一覧))
                {
                    return;
                }
                l_一覧.Add(new 曖昧箇所(_現在のk長, p_種別, p_場所, p_首位の支持, p_次点の支持, p_首位の生支持数, 証拠較正器.Get_飽和支持(p_首位の生支持数)));
            }
        }

        /// <summary>
        /// 符号付き頂点番号を、記録に残す読める名前にする
        /// </summary>
        /// <param name="p_頂点"></param>
        /// <param name="p_接頭辞"></param>
        /// <returns></returns>
        public static string Get_場所名(int p_頂点, string p_接頭辞 = 既定接頭辞)
        {
            return $"{p_接頭辞}{p_頂点 >> 1}{((p_頂点 & 1) == 0 ? '+' : '-')}";
        }

        /// <summary>
        /// その k の記録を作業ディレクトリへ残す
        /// </summary>
        /// <param name="p_作業ディレクトリ"></param>
        /// <param name="p_k長"></param>
        /// <remarks>
        /// 再開でこの k を飛ばしたときに、決めきれなかった箇所だけが失われてレポートが実態より綺麗に見えるのを防ぐ
        /// </remarks>
        public static void V_保存(string p_作業ディレクトリ, int p_k長)
        {
            var l_文 = new System.Text.StringBuilder();
            foreach (var l_箇所 in Get_記録(p_k長))
            {
                _ = l_文.AppendLine(string.Join('	', (int)l_箇所.A_種別, l_箇所.A_場所, l_箇所.A_首位の支持.ToString("R", CultureInfo.InvariantCulture), l_箇所.A_次点の支持.ToString("R", CultureInfo.InvariantCulture), l_箇所.A_首位の生支持数, l_箇所.A_確信度.ToString("R", CultureInfo.InvariantCulture)));
            }
            File.WriteAllText(Path.Combine(p_作業ディレクトリ, 保存ファイル名), l_文.ToString());
        }

        /// <summary>
        /// その k で書き留めた箇所の一覧
        /// </summary>
        /// <param name="p_k長"></param>
        /// <remarks>
        /// 記録が無ければ空
        /// </remarks>
        /// <returns></returns>
        public static IReadOnlyList<曖昧箇所> Get_記録(int p_k長)
        {
            lock (_錠)
            {
                return _k長ごとの記録.TryGetValue(p_k長, out var l_一覧) ? [.. l_一覧] : [];
            }
        }

        #endregion
    }
}
