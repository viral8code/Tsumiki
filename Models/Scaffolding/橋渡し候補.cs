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
    internal record 橋渡し候補(int A_始点, int A_終点, string A_橋渡し配列, int A_由来のk長, int A_重なり長 = 0);
}
