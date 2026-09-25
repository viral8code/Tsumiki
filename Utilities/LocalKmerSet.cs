using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 局所アセンブリ用の、ディスクを使わないインメモリの k-mer 集合
    /// </summary>
    internal sealed class LocalKmerSet : IKmerLookup
    {
        #region 内部変数

        /// <summary>
        /// 信頼できる k-mer 集合 (k &lt;= 32)
        /// </summary>
        private readonly HashSet<ulong>? _小;

        /// <summary>
        /// 信頼できる k-mer 集合 (33 &lt;= k &lt;= 64)
        /// </summary>
        private readonly HashSet<UInt128>? _中;

        /// <summary>
        /// 信頼できる k-mer 集合 (65 &lt;= k &lt;= 128)
        /// </summary>
        private readonly HashSet<(UInt128 A_上位, UInt128 A_下位)>? _長;

        /// <summary>
        /// 信頼できる k-mer 集合 (k &gt; 128)
        /// </summary>
        private readonly HashSet<KmerKey>? _大;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// k 長に応じた表現で空の集合を用意する
        /// </summary>
        /// <param name="p_k長"></param>
        public LocalKmerSet(int p_k長)
        {
            if (p_k長 <= 32)
            {
                this._小 = [];
            }
            else if (p_k長 <= 64)
            {
                this._中 = [];
            }
            else if (p_k長 <= TrustedKmerIndex.C_パック値のk上限)
            {
                this._長 = [];
            }
            else
            {
                this._大 = [];
            }
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// k-mer を 1 件登録する
        /// </summary>
        /// <param name="p_kmer"></param>
        public void V_登録(ReadOnlySpan<byte> p_kmer)
        {
            _ = this._小 is { } l_小
                ? l_小.Add(TrustedKmerIndex.Get_正規形_小(p_kmer))
                : this._中 is { } l_中
                ? l_中.Add(TrustedKmerIndex.Get_正規形_中(p_kmer))
                : this._長 is { } l_長 ? l_長.Add(TrustedKmerIndex.Get_正規形_長(p_kmer)) : this._大!.Add(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// kmer (順鎖・逆鎖いずれの向きでもよい) が集合に含まれるかどうかを判定する
        /// </summary>
        /// <param name="p_kmer"></param>
        /// <returns></returns>
        public bool Haskmer(Span<byte> p_kmer)
        {
            return this._小 is { } l_小
                ? l_小.Contains(TrustedKmerIndex.Get_正規形_小(p_kmer))
                : this._中 is { } l_中
                ? l_中.Contains(TrustedKmerIndex.Get_正規形_中(p_kmer))
                : this._長 is { } l_長 ? l_長.Contains(TrustedKmerIndex.Get_正規形_長(p_kmer)) : this._大!.Contains(new KmerKey(p_kmer).Get_正規形());
        }

        /// <summary>
        /// 正規形の右詰めパック値 (k &lt;= 128) が集合に含まれるかどうかを判定する
        /// </summary>
        /// <param name="p_上位">128 bit を超える側、k &lt;= 64 なら 0</param>
        /// <param name="p_下位"></param>
        /// <returns></returns>
        public bool Haskmer_正規形(UInt128 p_上位, UInt128 p_下位)
        {
            return this._小 is { } l_小
                ? l_小.Contains((ulong)p_下位)
                : this._中 is { } l_中 ? l_中.Contains(p_下位) : this._長!.Contains((p_上位, p_下位));
        }

        #endregion
    }
}
