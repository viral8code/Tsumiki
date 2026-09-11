namespace Tsumiki.Utilities
{
    /// <summary>
    /// k-mer が信頼できる集合に含まれるかどうかを判定できる型が実装するインターフェース
    /// </summary>
    /// <remarks>
    /// ConstrainedPathFinder はゲノム全体規模の TrustedKmerIndex と、
    /// 局所アセンブリ規模の軽量な LocalKmerSet のどちらが相手でも同じ探索ロジックを
    /// 使い回せるよう、この越しに k-mer 集合へ問い合わせる
    /// </remarks>
    internal interface IKmerLookup
    {
        #region 公開メソッド

        /// <summary>
        /// kmer(順鎖・逆鎖いずれの向きでもよい) が集合に含まれるかどうかを判定する
        /// </summary>
        /// <param name="p_kmer"></param>
        bool Get_含まれるか(Span<byte> p_kmer);

        #endregion
    }
}
