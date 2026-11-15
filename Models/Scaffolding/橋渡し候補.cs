namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 骨格アセンブリの 2 本の配列が隣接しているという、別の k のアセンブリからの証拠
    /// </summary>
    /// <param name="A_始点"></param>
    /// <param name="A_終点"></param>
    /// <param name="A_橋渡し配列"></param>
    /// <param name="A_由来のk長"></param>
    /// <param name="A_重なり長"></param>
    /// <param name="A_Is錨が反復だけ">片側の錨が反復の中だけに当たる橋渡しか<br/>true なら確定の候補にせず、行き先の競合の証拠にだけ使う</param>
    /// <param name="A_Is端が合わない">跨いだ配列が、骨格の端までの配列と片側では合い、もう片側では合わない橋渡しか<br/>true なら確定の候補にせず、行き先の競合の証拠にだけ使う</param>
    /// <param name="A_Is通り抜け">跨いだ配列が、間の骨格配列に短く触れて通り抜けた先への橋渡しか<br/>true なら確定の候補にせず、行き先の競合の証拠にだけ使う</param>
    /// <param name="A_届き">跨いだ配列が、前の骨格配列の端より手前と、後の骨格配列の端より先に持つ長さの短い方 (bp)<br/>確定の候補のどれもが短い辺は確定しない</param>
    internal record 橋渡し候補(int A_始点, int A_終点, string A_橋渡し配列, int A_由来のk長, int A_重なり長 = 0, bool A_Is錨が反復だけ = false, bool A_Is端が合わない = false, bool A_Is通り抜け = false, int A_届き = int.MaxValue);
}
