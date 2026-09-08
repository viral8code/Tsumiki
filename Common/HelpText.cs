using System.Text;
using Tsumiki.Model;

namespace Tsumiki.Common
{
    /// <summary>
    /// 概要表示(引数なし・-v)とヘルプ(-h)の本文を組み立てる。
    ///
    /// 1つの大きな文字列として持たせると言語ごとに全文を複製することになり、
    /// 既定値やオプション名の変更が全言語に波及する。ここでは「オプション名と
    /// 既定値はコード側、説明文だけカタログ側」に分けて、行ごとに組み立てる。
    /// </summary>
    internal static class HelpText
    {
        /// <summary>説明を書き始める桁。オプション名がこれを超える行は説明を次行へ送る。</summary>
        private const int 説明の開始桁 = 20;

        /// <summary>名前・作者・バージョン。引数なしの起動と -v で出す。</summary>
        public static string Get_概要()
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine();
            _ = l_文.AppendLine(Messages.Get_文言(メッセージID.概要_説明));
            _ = l_文.AppendLine(Messages.Get_文言(メッセージID.概要_作者, string.Join(", ", Consts.作者一覧)));
            _ = l_文.AppendLine(Messages.Get_文言(メッセージID.概要_バージョン, Consts.バージョン));
            return l_文.ToString();
        }

        public static string Get_ヘルプ()
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine(Get_概要());
            _ = l_文.AppendLine(Messages.Get_文言(
                メッセージID.ヘルプ_使い方, Consts.引数キー.リード1のパス, Consts.引数キー.リード2のパス));

            V_追加_節(l_文, メッセージID.ヘルプ節_入力);
            V_追加_行(l_文, $"{Consts.引数キー.リード1のパス} <path>", メッセージID.ヘルプ_リード1);
            V_追加_行(l_文, $"{Consts.引数キー.リード2のパス} <path>", メッセージID.ヘルプ_リード2);
            V_追加_行(l_文, Consts.引数キー.曖昧塩基を許容, メッセージID.ヘルプ_曖昧塩基);

            V_追加_節(l_文, メッセージID.ヘルプ節_kmerと品質);
            V_追加_行(
                l_文, $"{Consts.引数キー.k長} <int[,int...]>", メッセージID.ヘルプ_k長,
                Consts.引数キー.マルチk, Consts.自動k長の上限);
            V_追加_行(l_文, $"{Consts.引数キー.kmerカットオフ} <int>", メッセージID.ヘルプ_kmerカットオフ);
            V_追加_行(
                l_文, $"{Consts.引数キー.Phredオフセット} <int>", メッセージID.ヘルプ_Phred,
                string.Join(" or ", Consts.許容Phredオフセット), Consts.Phredオフセットの既定値);
            V_追加_行(
                l_文, $"{Consts.引数キー.クオリティカットオフ} <int>", メッセージID.ヘルプ_クオリティカットオフ,
                Consts.クオリティカットオフの既定値);
            V_追加_行(
                l_文, $"{Consts.引数キー.メモリ予算} <size>", メッセージID.ヘルプ_メモリ予算,
                Util.Get_表示用メモリサイズ(Consts.メモリ予算の既定値));

            V_追加_節(l_文, メッセージID.ヘルプ節_ペアエンド);
            V_追加_行(l_文, $"{Consts.引数キー.インサートサイズ} <int>", メッセージID.ヘルプ_インサートサイズ);
            V_追加_行(
                l_文, $"{Consts.引数キー.ペア結合閾値} <decimal>", メッセージID.ヘルプ_ペア結合閾値,
                Consts.ペア結合閾値の既定値);
            V_追加_行(
                l_文, $"{Consts.引数キー.ペア支持数閾値} <int>", メッセージID.ヘルプ_ペア支持数閾値,
                Consts.ペア支持数閾値の既定値);
            V_追加_行(
                l_文,
                $"{Consts.引数キー.積極性モード} <{Consts.積極性モード名.保守的}|{Consts.積極性モード名.標準}|{Consts.積極性モード名.積極的}>",
                メッセージID.ヘルプ_積極性モード,
                Consts.引数キー.ペア結合閾値, Consts.引数キー.ペア支持数閾値, Consts.積極性モード名.標準);

            V_追加_節(l_文, メッセージID.ヘルプ節_前処理);
            V_追加_行(l_文, Consts.引数キー.前処理, メッセージID.ヘルプ_前処理);
            V_追加_行(l_文, Consts.引数キー.エラー訂正, メッセージID.ヘルプ_エラー訂正);

            V_追加_節(l_文, メッセージID.ヘルプ節_マルチk);
            V_追加_行(l_文, Consts.引数キー.マルチk, メッセージID.ヘルプ_マルチk, Consts.マルチkで試す個数 + 1);
            V_追加_行(l_文, Consts.引数キー.引き継ぎなし, メッセージID.ヘルプ_引き継ぎなし);
            V_追加_行(l_文, Consts.引数キー.SuperRead, メッセージID.ヘルプ_SuperRead);
            V_追加_行(l_文, Consts.引数キー.マージ, メッセージID.ヘルプ_マージ);

            V_追加_節(l_文, メッセージID.ヘルプ節_反復の安全策);
            V_追加_行(l_文, Consts.引数キー.反復r_mer検証, メッセージID.ヘルプ_反復rMer検証);
            V_追加_行(l_文, Consts.引数キー.局所アセンブリ, メッセージID.ヘルプ_局所アセンブリ, Consts.引数キー.マージ);

            V_追加_節(l_文, メッセージID.ヘルプ節_出力とその他);
            V_追加_行(l_文, Consts.引数キー.GFA出力, メッセージID.ヘルプ_GFA出力, Consts.GFAファイル名);
            V_追加_行(
                l_文, $"{Consts.引数キー.一時ディレクトリ} <path>", メッセージID.ヘルプ_一時ディレクトリ,
                Consts.一時ディレクトリの既定値);
            V_追加_行(l_文, Consts.引数キー.一時ディレクトリ削除, メッセージID.ヘルプ_一時ディレクトリ削除);
            V_追加_行(l_文, $"{Consts.引数キー.スレッド数} <int>", メッセージID.ヘルプ_スレッド数);
            V_追加_行(
                l_文,
                $"{Consts.引数キー.言語} <{Consts.言語名.日本語}|{Consts.言語名.英語}|{Consts.言語名.中国語}>",
                メッセージID.ヘルプ_言語, Consts.言語名.日本語);
            V_追加_行(l_文, Consts.引数キー.バージョン, メッセージID.ヘルプ_バージョン);
            V_追加_行(l_文, Consts.引数キー.ヘルプ, メッセージID.ヘルプ_ヘルプ);

            return l_文.ToString();
        }

        private static void V_追加_節(StringBuilder p_文, メッセージID p_見出し)
        {
            _ = p_文.AppendLine();
            _ = p_文.AppendLine(Messages.Get_文言(p_見出し));
        }

        private static void V_追加_行(
            StringBuilder p_文, string p_オプション, メッセージID p_説明, params object?[] p_引数)
        {
            var l_説明 = Messages.Get_文言(p_説明, p_引数);
            _ = p_オプション.Length < 説明の開始桁
                ? p_文.AppendLine(p_オプション.PadRight(説明の開始桁) + l_説明)
                : p_文.AppendLine(p_オプション).AppendLine(new string(' ', 説明の開始桁) + l_説明);
        }
    }
}
