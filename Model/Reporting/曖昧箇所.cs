namespace Tsumiki.Model.Reporting
{
    /// <summary>
    /// 決めきれずに打ち切った箇所の区分
    /// </summary>
    internal enum 曖昧箇所の種別
    {
        /// <summary>
        /// どの候補にも証拠が無く、伸ばす根拠が得られなかった
        /// </summary>
        支持なし,

        /// <summary>
        /// 首位と次点の差が小さく、一方に確定できなかった
        /// </summary>
        僅差,

        /// <summary>
        /// 両端を繋ぐ経路が複数あり、どれか 1 本に定まらなかった
        /// </summary>
        経路が一意でない,

        /// <summary>
        /// 両端を繋ぐ経路がグラフ上に見つからなかった
        /// </summary>
        到達不能,

        /// <summary>
        /// 分岐元が多コピーで、いまどのコピーの上にいるのかを区別できなかった
        /// </summary>
        /// <remarks>
        /// 反復の内側から読まれたリードはどの行き先にも支持を付けるため、
        /// 選ぶ根拠が原理的に存在しない
        /// </remarks>
        反復の内側,
    }

    /// <summary>
    /// 決めきれなかった 1 箇所の記録
    /// </summary>
    /// <remarks>
    /// N で埋めて黙って落とすと、後から
    /// 人や別のツールが再解析するための手掛かりが残らない
    /// </remarks>
    internal readonly record struct 曖昧箇所(
        int A_k長,
        曖昧箇所の種別 A_種別,
        string A_場所,
        double A_首位の支持,
        double A_次点の支持,
        long A_首位の生支持数,
        double A_確信度)
    {
        /// <summary>
        /// 首位と次点の差
        /// </summary>
        /// <remarks>
        /// これが小さいほど選ぶ根拠が薄い
        /// </remarks>
        public double A_余裕 => this.A_首位の支持 - this.A_次点の支持;
    }
}
