namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 1 つのギャップを埋めようとした結果の区分
    /// </summary>
    internal enum ギャップ充填判定
    {
        /// <summary>
        /// 経路が一意に定まり、実配列で埋められた
        /// </summary>
        充填済み,

        /// <summary>
        /// 経路が複数見つかり、どれが正しいか決められなかった
        /// </summary>
        一意でない,

        /// <summary>
        /// 両端を繋ぐ経路がグラフ上に存在しなかった
        /// </summary>
        到達不能,

        /// <summary>
        /// 探索状態数が上限に達し、経路が存在するかどうかも決められないまま打ち切った
        /// </summary>
        /// <remarks>
        /// 到達不能(経路が無いと確認できた)とは異なり、時間をかければ見つかったかもしれない
        /// </remarks>
        探索打切り,
    }
}
