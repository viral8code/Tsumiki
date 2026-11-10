namespace Tsumiki.Models.Scaffolding
{
    /// <summary>
    /// 継ぎ目をリードの合意で埋めたかどうかの判定
    /// </summary>
    internal enum 継ぎ目合意判定
    {
        /// <summary>
        /// 埋めた
        /// </summary>
        埋めた,

        /// <summary>
        /// 両側の片に置けるリードが無い
        /// </summary>
        跨ぐリード無し,

        /// <summary>
        /// 最多の間の長さを示すリードが要る数に満たない
        /// </summary>
        票が足りない,

        /// <summary>
        /// 最多の間の長さが次点の優勢比倍に届かない
        /// </summary>
        票が割れた,

        /// <summary>
        /// 列ごとの多数決で決まらない列がある
        /// </summary>
        合意が弱い,

        /// <summary>
        /// 繋いだ配列が縦に並んだ反復の中
        /// </summary>
        縦列反復,

        /// <summary>
        /// 片側だけに置けるリードの多くが、繋いだ形のもう一方の片と食い違う
        /// </summary>
        片側のリードが食い違う,

        /// <summary>
        /// 重なりが短いほうの片の長さ以上で、片がまるごともう一方の片の範囲に入る
        /// </summary>
        片が重なりに収まる,
    }
}
