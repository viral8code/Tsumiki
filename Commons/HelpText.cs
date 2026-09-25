using System.Text;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Commons
{
    /// <summary>
    /// 概要表示 (引数なし・-v) とヘルプ (-h) の本文を組み立てる
    /// </summary>
    internal static class HelpText
    {
        #region 定数

        /// <summary>
        /// 説明を書き始める桁
        /// </summary>
        private const int 説明の開始桁 = 20;

        /// <summary>
        /// 入力形式 (パス)
        /// </summary>
        private const string 形式_パス = "path";

        /// <summary>
        /// 入力形式 (整数)
        /// </summary>
        private const string 形式_整数 = "int";

        /// <summary>
        /// 入力形式 (小数)
        /// </summary>
        private const string 形式_小数 = "decimal";

        /// <summary>
        /// 入力形式 (容量)
        /// </summary>
        private const string 形式_容量 = "size";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 名前とバージョン
        /// </summary>
        /// <returns></returns>
        public static string Get_概要()
        {
            return new StringBuilder()
                .AppendLine()
                .AppendLine(Messages.Get_文言(メッセージID.概要_説明))
                .AppendLine(Messages.Get_文言(メッセージID.概要_バージョン, Consts.バージョン))
                .ToString();
        }

        /// <summary>
        /// ヘルプの全文を組み立てて返す
        /// </summary>
        /// <returns>ヘルプの全文</returns>
        public static string Get_ヘルプ()
        {
            var l_文 = new StringBuilder(Get_概要())
                .AppendLine(Messages.Get_文言(メッセージID.ヘルプ_使い方, Consts.引数キー.リード1のパス, Consts.引数キー.リード2のパス));

            V_追加_節(l_文, メッセージID.ヘルプ節_入力);
            V_追加_行(l_文, $"{Consts.引数キー.リード1のパス} <{形式_パス}>", メッセージID.ヘルプ_リード1);
            V_追加_行(l_文, $"{Consts.引数キー.リード2のパス} <{形式_パス}>", メッセージID.ヘルプ_リード2);
            V_追加_行(l_文, $"{Consts.引数キー.シングルのパス} <{形式_パス}>", メッセージID.ヘルプ_シングル);
            V_追加_行(l_文, Consts.引数キー.曖昧塩基を許容, メッセージID.ヘルプ_曖昧塩基);

            V_追加_節(l_文, メッセージID.ヘルプ節_kmerと品質);
            V_追加_行(l_文, $"{Consts.引数キー.k長} <{形式_整数}[,{形式_整数}...]>", メッセージID.ヘルプ_k長, Consts.自動k長の上限);
            V_追加_行(l_文, $"{Consts.引数キー.kmerカットオフ} <{形式_整数}>", メッセージID.ヘルプ_kmerカットオフ);
            V_追加_行(l_文, $"{Consts.引数キー.Phredオフセット} <{形式_整数}>", メッセージID.ヘルプ_Phred, string.Join(" or ", Consts.許容Phredオフセット), Consts.Phredオフセットの既定値);
            V_追加_行(l_文, $"{Consts.引数キー.クオリティカットオフ} <{形式_整数}>", メッセージID.ヘルプ_クオリティカットオフ, Consts.クオリティカットオフの既定値);
            V_追加_行(l_文, $"{Consts.引数キー.品質トリム閾値} <{形式_整数}>", メッセージID.ヘルプ_品質トリム閾値, Consts.品質トリム閾値の既定値);
            V_追加_行(l_文, $"{Consts.引数キー.メモリ予算} <{形式_容量}>", メッセージID.ヘルプ_メモリ予算, Util.Get_表示用メモリサイズ(Consts.メモリ予算の既定値));
            V_追加_行(l_文, Consts.引数キー.救済kmerなし, メッセージID.ヘルプ_救済kmer);
            V_追加_行(l_文, $"{Consts.引数キー.コピー数基準} <{Consts.コピー数基準の出所.スペクトラム}|{Consts.コピー数基準の出所.重みづけ}>", メッセージID.ヘルプ_コピー数基準);
            V_追加_行(l_文, Consts.引数キー.低カバレッジ端トリミングなし, メッセージID.ヘルプ_低カバレッジ端トリミングなし);

            V_追加_節(l_文, メッセージID.ヘルプ節_ペアエンド);
            V_追加_行(l_文, $"{Consts.引数キー.インサートサイズ} <{形式_整数}>", メッセージID.ヘルプ_インサートサイズ);
            V_追加_行(l_文, $"{Consts.引数キー.ペア結合閾値} <{形式_小数}>", メッセージID.ヘルプ_ペア結合閾値, Consts.ペア結合閾値の既定値);
            V_追加_行(l_文, $"{Consts.引数キー.ペア支持数閾値} <{形式_整数}>", メッセージID.ヘルプ_ペア支持数閾値, Consts.ペア支持数閾値の既定値);
            V_追加_行(l_文, $"{Consts.引数キー.積極性モード} <{Consts.積極性モード名.保守的}|{Consts.積極性モード名.標準}|{Consts.積極性モード名.積極的}>", メッセージID.ヘルプ_積極性モード, Consts.引数キー.ペア結合閾値, Consts.引数キー.ペア支持数閾値, Consts.積極性モード名.標準);

            V_追加_節(l_文, メッセージID.ヘルプ節_前処理);
            V_追加_行(l_文, Consts.引数キー.前処理なし, メッセージID.ヘルプ_前処理);
            V_追加_行(l_文, Consts.引数キー.エラー訂正なし, メッセージID.ヘルプ_エラー訂正);

            V_追加_節(l_文, メッセージID.ヘルプ節_マルチk);
            V_追加_行(l_文, Consts.引数キー.マルチkなし, メッセージID.ヘルプ_マルチk, Consts.マルチkで試す個数 + 1);
            V_追加_行(l_文, Consts.引数キー.引き継ぎなし, メッセージID.ヘルプ_引き継ぎなし);
            V_追加_行(l_文, Consts.引数キー.SuperReadなし, メッセージID.ヘルプ_SuperRead);
            V_追加_行(l_文, Consts.引数キー.マージ, メッセージID.ヘルプ_マージ);

            V_追加_節(l_文, メッセージID.ヘルプ節_反復の安全策);
            V_追加_行(l_文, Consts.引数キー.反復r_mer検証なし, メッセージID.ヘルプ_反復rMer検証);
            V_追加_行(l_文, Consts.引数キー.局所アセンブリなし, メッセージID.ヘルプ_局所アセンブリ, Consts.引数キー.マージ);

            V_追加_節(l_文, メッセージID.ヘルプ節_完全性の検証);
            V_追加_行(l_文, Consts.引数キー.ポリッシュなし, メッセージID.ヘルプ_ポリッシュ);
            V_追加_行(l_文, Consts.引数キー.環状閉鎖検証なし, メッセージID.ヘルプ_環状閉鎖検証);
            V_追加_説明(l_文, メッセージID.ヘルプ_レポートの説明, Consts.レポートファイル名, Consts.曖昧箇所ファイル名);

            V_追加_節(l_文, メッセージID.ヘルプ節_出力とその他);
            V_追加_行(l_文, Consts.引数キー.GFA出力なし, メッセージID.ヘルプ_GFA出力, Consts.GFAファイル名);
            V_追加_行(l_文, $"{Consts.引数キー.一時ディレクトリ} <{形式_パス}>", メッセージID.ヘルプ_一時ディレクトリ, Consts.一時ディレクトリの既定値);
            V_追加_行(l_文, Consts.引数キー.一時ディレクトリ削除, メッセージID.ヘルプ_一時ディレクトリ削除);
            V_追加_行(l_文, Consts.引数キー.再開, メッセージID.ヘルプ_再開);
            V_追加_行(l_文, Consts.引数キー.オンメモリ, メッセージID.ヘルプ_オンメモリ);
            V_追加_行(l_文, $"{Consts.引数キー.スレッド数} <{形式_整数}>", メッセージID.ヘルプ_スレッド数);
            V_追加_行(l_文, $"{Consts.引数キー.言語} <{Consts.言語名.日本語}|{Consts.言語名.英語}|{Consts.言語名.中国語}>", メッセージID.ヘルプ_言語, Consts.言語名.日本語);
            V_追加_行(l_文, $"{Consts.引数キー.ログ水準} <{Consts.ログ水準名.最小}|{Consts.ログ水準名.標準}|{Consts.ログ水準名.詳細}>", メッセージID.ヘルプ_ログ水準, Consts.ログ水準名.標準, Consts.ログファイル名);
            V_追加_行(l_文, Consts.引数キー.バージョン, メッセージID.ヘルプ_バージョン);
            V_追加_行(l_文, Consts.引数キー.ヘルプ, メッセージID.ヘルプ_ヘルプ);

            return l_文.ToString();
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 節の見出しを 1 行足す
        /// </summary>
        /// <param name="p_文">組み立て中のヘルプ</param>
        /// <param name="p_見出し">見出しの文言</param>
        private static void V_追加_節(StringBuilder p_文, メッセージID p_見出し)
        {
            _ = p_文.AppendLine().AppendLine(Messages.Get_文言(p_見出し));
        }

        /// <summary>
        /// オプションに紐づかない補足を、説明と同じ桁に 1 行足す
        /// </summary>
        /// <param name="p_文"></param>
        /// <param name="p_説明"></param>
        /// <param name="p_引数"></param>
        private static void V_追加_説明(StringBuilder p_文, メッセージID p_説明, params object?[] p_引数)
        {
            _ = p_文.AppendLine(new string(' ', 説明の開始桁) + Messages.Get_文言(p_説明, p_引数));
        }

        /// <summary>
        /// オプション 1 つぶんの説明を足す
        /// </summary>
        /// <param name="p_文">組み立て中のヘルプ</param>
        /// <param name="p_オプション">オプションの綴り</param>
        /// <param name="p_説明">説明の文言</param>
        /// <param name="p_引数">説明へ埋め込む値</param>
        private static void V_追加_行(StringBuilder p_文, string p_オプション, メッセージID p_説明, params object?[] p_引数)
        {
            var l_説明 = Messages.Get_文言(p_説明, p_引数);
            _ = p_オプション.Length < 説明の開始桁
                ? p_文.AppendLine(p_オプション.PadRight(説明の開始桁) + l_説明)
                : p_文.AppendLine(p_オプション).AppendLine(new string(' ', 説明の開始桁) + l_説明);
        }

        #endregion
    }
}
