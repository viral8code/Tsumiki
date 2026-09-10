namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// 画面へ出す量の段階
    /// </summary>
    /// <remarks>
    /// 段階の判定は文言の先頭に付く目印 ([Debug] など) で行う<br/>
    /// 目印は
    /// 言語によらず同じ綴りで、その行が何のための行なのかを表しているため、
    /// 別に対応表を持たずに済む<br/>
    /// ファイルへの記録はこの設定に関わらず常に全量を残す<br/>
    /// 画面を静かに
    /// したことで、後から原因を追う手掛かりまで失われないようにする
    /// </remarks>
    internal enum ログ水準
    {
        /// <summary>
        /// 結論だけ
        /// </summary>
        /// <remarks>
        /// 警告・エラーと、完全性の判定・レポートの出力先
        /// </remarks>
        最小 = 0,

        /// <summary>
        /// 進行状況と結果
        /// </summary>
        /// <remarks>
        /// 内部の判断過程は出さない
        /// </remarks>
        標準 = 1,

        /// <summary>
        /// 内部の判断過程まで含めて全部
        /// </summary>
        詳細 = 2,
    }
}
