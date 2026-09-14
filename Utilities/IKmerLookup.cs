namespace Tsumiki.Utilities
{
    /// <summary>
    /// k-mer が信頼できる集合に含まれるかどうかを判定できる型が実装するインターフェース
    /// </summary>
    /// <remarks>
    /// ConstrainedPathFinder はゲノム全体規模の TrustedKmerIndex と、局所アセンブリ規模の軽量な LocalKmerSet のどちらが相手でも同じ探索ロジックを使い回せるよう、この越しに k-mer 集合へ問い合わせる
    /// </remarks>
    internal interface IKmerLookup
    {
        #region 公開メソッド

        /// <summary>
        /// kmer (順鎖・逆鎖いずれの向きでもよい) が集合に含まれるかどうかを判定する
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        bool Haskmer(Span<byte> p_kmer);

        /// <summary>
        /// 正規形の右詰めパック値 (k &lt;= 128) が集合に含まれるかどうかを判定する
        /// </summary>
        /// <param name="p_上位">128 bit を超える側、k &lt;= 64 なら 0</param>
        /// <param name="p_下位"></param>
        /// <returns></returns>
        bool Haskmer_正規形(UInt128 p_上位, UInt128 p_下位);

        #endregion
    }
}
