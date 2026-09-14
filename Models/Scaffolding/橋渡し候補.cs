namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 骨格アセンブリの 2 本の配列が隣接しているという、別の k のアセンブリからの証拠
    /// </summary>
    /// <remarks>
    /// 頂点番号は骨格配列 ID に向きを付けたもの (2 i が順鎖、2 i+1 が逆鎖、双子は v^1) <br/>
    /// A_橋渡し配列は 2 本の間に挟まる塩基<br/>
    /// 2 本が重なって繋がる場合は橋渡し配列が空になり、終点側の先頭を A_重なり長 だけ削って繋ぐ
    /// </remarks>
    /// <param name="A_始点"></param>
    /// <param name="A_終点"></param>
    /// <param name="A_橋渡し配列"></param>
    /// <param name="A_由来のk長"></param>
    /// <param name="A_重なり長"></param>
    internal record 橋渡し候補(int A_始点, int A_終点, string A_橋渡し配列, int A_由来のk長, int A_重なり長 = 0);
}
