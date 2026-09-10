namespace Tsumiki.Models.Preprocessing
{
    /// <summary>
    /// ペアの重なり解析 (R1 と RC (R2)) が見つけた最良の位置合わせ
    /// </summary>
    /// <remarks>
    /// A_offset は RC (R2) の先頭が R1 の座標系で来る位置 (R1 の座標系がそのまま
    /// 断片上の座標系になる)<br/>
    /// フラグメント長は A_offset + R2 の長さ で求まる
    /// </remarks>
    internal readonly record struct オーバーラップ結果(int A_offset, int A_重なり長, int A_不一致数);
}
