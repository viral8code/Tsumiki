namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// k-mer 出現回数ヒストグラム (k-mer スペクトル) の解析結果
    /// </summary>
    /// <remarks>
    /// 実データのスペクトルは二峰性になる<br/>
    /// 出現回数 1 付近にシーケンスエラー由来の巨大な山があり、そこから離れたところに真のゲノム由来の山がある<br/>
    /// 前者と後者を分ける「谷」がカットオフの適正値であり、後者の山の位置が 1 コピーあたりの k-mer カバレッジになる
    /// </remarks>
    /// <param name="A_谷"></param>
    /// <param name="A_ピーク出現回数"></param>
    /// <param name="A_谷の頻度"></param>
    /// <param name="A_ピークの頻度"></param>
    /// <param name="A_ゲノム由来の延べ数"></param>
    /// <param name="A_延べ数の総和"></param>
    /// <param name="A_推定ゲノムサイズ"></param>
    internal record スペクトル解析結果(ulong A_谷, ulong A_ピーク出現回数, long A_谷の頻度, long A_ピークの頻度, long A_ゲノム由来の延べ数, long A_延べ数の総和, long A_推定ゲノムサイズ);
}
