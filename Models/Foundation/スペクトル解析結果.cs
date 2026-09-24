namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// k-mer 出現回数ヒストグラム (k-mer スペクトル) の解析結果
    /// </summary>
    /// <param name="A_谷"></param>
    /// <param name="A_ピーク出現回数"></param>
    /// <param name="A_谷の頻度"></param>
    /// <param name="A_ピークの頻度"></param>
    /// <param name="A_ゲノム由来の延べ数"></param>
    /// <param name="A_延べ数の総和"></param>
    /// <param name="A_推定ゲノムサイズ"></param>
    /// <param name="A_単一コピー基準値">単一コピーの k-mer カバレッジ</param>
    /// <param name="A_単一コピー上限">これ未満のカバレッジを単一コピーとみなす</param>
    internal record スペクトル解析結果(ulong A_谷, ulong A_ピーク出現回数, long A_谷の頻度, long A_ピークの頻度, long A_ゲノム由来の延べ数, long A_延べ数の総和, long A_推定ゲノムサイズ, double A_単一コピー基準値, double A_単一コピー上限);
}
