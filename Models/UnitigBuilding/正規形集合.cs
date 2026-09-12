using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// 逆相補を同一視して k-mer を覚える集合
    /// </summary>
    /// <remarks>
    /// k で表現を切り替えるのは、k &lt;= 64 なら 2 bit パックが UInt128 に収まり、鍵 1 つあたりの大きさが半分以下になるため
    /// </remarks>
    /// <param name="p_k長">k 長</param>
    internal sealed class 正規形集合(int p_k長)
    {
        #region 内部変数

        /// <summary>
        /// 小
        /// </summary>
        private readonly HashSet<UInt128>? _小 = p_k長 <= 64 ? [] : null;

        /// <summary>
        /// 大
        /// </summary>
        private readonly HashSet<KmerKey>? _大 = p_k長 > 64 ? [] : null;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 訪問済みとして正規形を覚える
        /// </summary>
        /// <param name="p_正規形">k-mer の正規形</param>
        public void V_追加(UInt128 p_正規形)
        {
            _ = this._小!.Add(p_正規形);
        }

        /// <summary>
        /// 訪問済みとして k-mer を覚える
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        public void V_追加(ReadOnlySpan<byte> p_kmer)
        {
            if (this._小 is { } l_小)
            {
                _ = l_小.Add(KmerPacking.TryGet_正規化パック(p_kmer));
                return;
            }
            _ = this._大!.Add(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// 既に訪問済みか
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列</param>
        /// <returns>訪問済みなら true</returns>
        public bool Haskmer(ReadOnlySpan<byte> p_kmer)
        {
            return this._小 is { } l_小
                ? l_小.Contains(KmerPacking.TryGet_正規化パック(p_kmer))
                : this._大!.Contains(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// 2 つの k-mer が、逆相補を同一視して同じ座位を指すか
        /// </summary>
        /// <param name="p_左">比較する k-mer の片方</param>
        /// <param name="p_右">比較する k-mer のもう片方</param>
        /// <param name="p_k長">k 長</param>
        /// <returns>同じ座位を指すなら true</returns>
        public static bool Is同一座位(ReadOnlySpan<byte> p_左, ReadOnlySpan<byte> p_右, int p_k長)
        {
            return p_k長 <= 64
                ? KmerPacking.TryGet_正規化パック(p_左) == KmerPacking.TryGet_正規化パック(p_右)
                : new KmerKey(p_左).Get_正規形().Equals(new KmerKey(p_右).Get_正規形());
        }

        #endregion
    }
}
