namespace Tsumiki.Models.ContigBuilding
{
    /// <summary>
    /// リードが通った unitig の並びを辞書のキーにするための値
    /// </summary>
    /// <param name="p_頂点列">符号付き unitig ID の並び</param>
    internal readonly struct 経路キー(int[] p_頂点列) : IEquatable<経路キー>
    {
        #region 内部変数

        /// <summary>
        /// 符号付き unitig ID の並び
        /// </summary>
        public readonly int[] A_頂点列 = p_頂点列;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 並びが要素ごとに一致するか
        /// </summary>
        /// <param name="p_他"></param>
        /// <returns></returns>
        public bool Equals(経路キー p_他)
        {
            return this.A_頂点列.AsSpan().SequenceEqual(p_他.A_頂点列);
        }

        #endregion

        #region 継承メソッド

        /// <summary>
        /// (オーバーライド) 並びが要素ごとに一致するか
        /// </summary>
        /// <param name="p_他"></param>
        /// <returns></returns>
        public override bool Equals(object? p_他)
        {
            return p_他 is 経路キー l_キー && this.Equals(l_キー);
        }

        /// <summary>
        /// (オーバーライド) 並び全体から求めたハッシュ値
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            var l_ハッシュ = new HashCode();
            foreach (var l_ID in this.A_頂点列)
            {
                l_ハッシュ.Add(l_ID);
            }

            return l_ハッシュ.ToHashCode();
        }

        #endregion
    }
}
