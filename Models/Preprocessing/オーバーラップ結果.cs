namespace Tsumiki.Models.Preprocessing
{
    /// <summary>
    /// ペアの重なり解析 (R1 と RC (R2)) が見つけた最良の位置合わせ
    /// </summary>
    /// <param name="A_offset"></param>
    /// <param name="A_重なり長"></param>
    /// <param name="A_不一致数"></param>
    internal readonly record struct オーバーラップ結果(int A_offset, int A_重なり長, int A_不一致数);
}
