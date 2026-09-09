using System.Globalization;
using Tsumiki.Core.Evidence;
using Tsumiki.Model.Reporting;

namespace Tsumiki.Core.Evaluation
{
    /// <summary>
    /// 決めきれずに打ち切った箇所を、その場で書き留めておくための収集器<br/>
    /// 打ち切った事実は「N を入れた」「繋がなかった」という結果にしか残らず、
    /// なぜそこで止めたのかは失われる<br/>
    /// 後から人や別のツールが再解析できるよう、
    /// 判断に使った数値ごと残す<br/>
    /// k ごとに分けて持つのは、multi-k では採用しなかった k の記録まで
    /// 混ざるため<br/>
    /// レポートには採用した k のぶんだけを出す
    /// </summary>
    internal static class AmbiguityRecorder
    {
        private static readonly object _錠 = new();

        private static readonly Dictionary<int, List<曖昧箇所>> _k長ごとの記録 = [];

        private static int _現在のk長;

        /// <summary>
        /// k ごとの作業ディレクトリに残す控え<br/>
        /// 再開時に読み直す
        /// </summary>
        private const string 保存ファイル名 = "ambiguous.tsv";

        /// <summary>
        /// この k の記録を集め直す<br/>
        /// 同じ k を再実行した場合は上書きする
        /// </summary>
        public static void V_開始(int p_k長)
        {
            lock (_錠)
            {
                _現在のk長 = p_k長;
                _k長ごとの記録[p_k長] = [];
            }
        }

        /// <summary>
        /// 1 箇所ぶん書き留める<br/>
        /// 確信度は独立な支持本数を飽和関数に通した値で、
        /// 同じ種類の証拠がいくら積み上がっても 1 に近づくだけになる
        /// </summary>
        public static void V_記録(
            曖昧箇所の種別 p_種別, string p_場所,
            double p_首位の支持 = 0D, double p_次点の支持 = 0D, long p_首位の生支持数 = 0L)
        {
            lock (_錠)
            {
                if (!_k長ごとの記録.TryGetValue(_現在のk長, out var l_一覧))
                {
                    return;
                }
                l_一覧.Add(new 曖昧箇所(
                    _現在のk長, p_種別, p_場所,
                    p_首位の支持, p_次点の支持, p_首位の生支持数,
                    証拠較正器.Get_飽和支持(p_首位の生支持数)));
            }
        }

        /// <summary>
        /// 符号付き頂点番号を、記録に残す読める名前にする
        /// </summary>
        public static string Get_場所名(int p_頂点, string p_接頭辞 = "unitig")
        {
            return $"{p_接頭辞}{p_頂点 >> 1}{((p_頂点 & 1) == 0 ? '+' : '-')}";
        }

        /// <summary>
        /// その k の記録を作業ディレクトリへ残す<br/>
        /// 再開でこの k を飛ばしたときに、
        /// 決めきれなかった箇所だけが失われてレポートが実態より綺麗に見えるのを防ぐ
        /// </summary>
        public static void V_保存(string p_作業ディレクトリ, int p_k長)
        {
            var l_文 = new System.Text.StringBuilder();
            foreach (var l_箇所 in Get_記録(p_k長))
            {
                _ = l_文.AppendLine(string.Join(
                    '	',
                    (int)l_箇所.A_種別,
                    l_箇所.A_場所,
                    l_箇所.A_首位の支持.ToString("R", CultureInfo.InvariantCulture),
                    l_箇所.A_次点の支持.ToString("R", CultureInfo.InvariantCulture),
                    l_箇所.A_首位の生支持数,
                    l_箇所.A_確信度.ToString("R", CultureInfo.InvariantCulture)));
            }
            File.WriteAllText(Path.Combine(p_作業ディレクトリ, 保存ファイル名), l_文.ToString());
        }

        /// <summary>
        /// その k で書き留めた箇所の一覧<br/>
        /// 記録が無ければ空
        /// </summary>
        public static IReadOnlyList<曖昧箇所> Get_記録(int p_k長)
        {
            lock (_錠)
            {
                return _k長ごとの記録.TryGetValue(p_k長, out var l_一覧) ? [.. l_一覧] : [];
            }
        }
    }
}
