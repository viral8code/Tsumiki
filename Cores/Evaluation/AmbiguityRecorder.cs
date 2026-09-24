using System.Globalization;
using System.Security.Cryptography;
using System.Text;
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
        private const string 保存ファイル名 = "ambiguous.tsv";

        /// <summary>
        /// 全 k を通じた累積履歴のファイル名
        /// </summary>
        private const string 履歴ファイル名 = "ambiguous.history.tsv";

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
        /// 再実行で上書きされる前の判定も含めた、全 k を通じた累積履歴
        /// </summary>
        private static readonly List<曖昧箇所> _履歴 = [];

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
        /// <param name="p_安定ID">再実行や ID 採番替えをまたいで同じ箇所を追跡するための ID、無ければ空</param>
        public static void V_記録(曖昧箇所の種別 p_種別, string p_場所, double p_首位の支持 = 0D, double p_次点の支持 = 0D, long p_首位の生支持数 = 0L, string p_安定ID = "")
        {
            lock (_錠)
            {
                var l_記録 = new 曖昧箇所(_現在のk長, p_種別, p_場所, p_首位の支持, p_次点の支持, p_首位の生支持数, 証拠較正器.Get_飽和支持(p_首位の生支持数), p_安定ID);
                _履歴.Add(l_記録);
                if (_k長ごとの記録.TryGetValue(_現在のk長, out var l_一覧))
                {
                    l_一覧.Add(l_記録);
                }
            }
        }

        /// <summary>
        /// 左右アンカー配列から、座標に依らない安定 ID を作る
        /// </summary>
        /// <param name="p_左アンカー">左側の足場配列</param>
        /// <param name="p_右アンカー">右側の足場配列</param>
        /// <returns>16 文字の16進ID</returns>
        public static string Get_安定ID(string p_左アンカー, string p_右アンカー)
        {
            var l_ハッシュ = SHA256.HashData(Encoding.ASCII.GetBytes(p_左アンカー + "|" + p_右アンカー));
            return Convert.ToHexStringLower(l_ハッシュ)[..16];
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
        public static void V_保存(string p_作業ディレクトリ, int p_k長)
        {
            File.WriteAllText(Path.Combine(p_作業ディレクトリ, 保存ファイル名), Get_行群(Get_記録(p_k長)));
        }

        /// <summary>
        /// 再実行で上書きされた分も含む累積履歴を作業ディレクトリへ残す
        /// </summary>
        /// <param name="p_作業ディレクトリ"></param>
        public static void V_保存_履歴(string p_作業ディレクトリ)
        {
            lock (_錠)
            {
                File.WriteAllText(Path.Combine(p_作業ディレクトリ, 履歴ファイル名), Get_行群(_履歴));
            }
        }

        /// <summary>
        /// 記録を TSV の本文へ整形する
        /// </summary>
        /// <param name="p_記録"></param>
        /// <returns></returns>
        private static string Get_行群(IReadOnlyList<曖昧箇所> p_記録)
        {
            var l_文 = new System.Text.StringBuilder();
            foreach (var l_箇所 in p_記録)
            {
                _ = l_文.AppendLine(string.Join('	', (int)l_箇所.A_種別, l_箇所.A_場所, l_箇所.A_安定ID, l_箇所.A_首位の支持.ToString("R", CultureInfo.InvariantCulture), l_箇所.A_次点の支持.ToString("R", CultureInfo.InvariantCulture), l_箇所.A_首位の生支持数, l_箇所.A_確信度.ToString("R", CultureInfo.InvariantCulture)));
            }
            return l_文.ToString();
        }

        /// <summary>
        /// その k で書き留めた箇所の一覧
        /// </summary>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        public static IReadOnlyList<曖昧箇所> Get_記録(int p_k長)
        {
            lock (_錠)
            {
                return _k長ごとの記録.TryGetValue(p_k長, out var l_一覧) ? [.. l_一覧] : [];
            }
        }

        #endregion

        #region テストメソッド

        /// <summary>
        /// 再実行で上書きされた分も含む、全 k を通じた累積履歴
        /// </summary>
        /// <returns></returns>
        public static IReadOnlyList<曖昧箇所> Get_履歴()
        {
            lock (_錠)
            {
                return [.. _履歴];
            }
        }

        #endregion
    }
}
