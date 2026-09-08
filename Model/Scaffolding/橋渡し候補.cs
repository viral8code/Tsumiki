namespace Tsumiki.Model
{
    /// <summary>
    /// 骨格アセンブリの2本の配列が隣接しているという、別の k のアセンブリからの証拠。
    ///
    /// 頂点番号は骨格配列 ID に向きを付けたもの(2i が順鎖、2i+1 が逆鎖、双子は v^1)。
    /// A_橋渡し配列 は2本の間に挟まる塩基。重なりだけで繋がる場合は空になる。
    /// </summary>
    internal record 橋渡し候補(
        int A_始点,
        int A_終点,
        string A_橋渡し配列,
        int A_由来のk長);
}
