namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// ペアで橋渡し (PairBridgeFiller) の結果の統計
    /// </summary>
    /// <param name="A_挿入長">挿入長の標本の中央値、標本不足なら 0</param>
    /// <param name="A_候補の辺の数">相互に一意で、間の見込みが範囲内の端の組の数</param>
    /// <param name="A_埋めて繋いだ辺の数">局所組み立てで埋めて繋いだ辺の数 (鎖の連結に使った数)</param>
    /// <param name="A_一意でなかった数">局所組み立てで経路が一意に定まらなかった候補の数</param>
    /// <param name="A_届かなかった数">局所組み立てで埋まらなかった候補の数 (一意でない以外)</param>
    /// <param name="A_印で繋いだ辺の数">埋まらなかったが票が多いため、未確認の繋ぎ目 (N) で繋いだ辺の数</param>
    /// <param name="A_反復の端で除いた数">端の配列が反復の端と見なされ、繋がずに除いた候補の辺の数</param>
    /// <param name="A_重なりで畳んだ辺の数">局所組み立てで埋まらず、配列の一致による重なりで畳んだ辺の数 (埋めた数には含めない)</param>
    internal readonly record struct 橋渡し統計(int A_挿入長, int A_候補の辺の数, int A_埋めて繋いだ辺の数, int A_一意でなかった数, int A_届かなかった数, int A_印で繋いだ辺の数, int A_反復の端で除いた数, int A_重なりで畳んだ辺の数);
}
