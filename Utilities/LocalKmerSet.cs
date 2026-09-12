using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 局所アセンブリ用の、ディスクを使わないインメモリの k-mer 集合
    /// </summary>
    /// <remarks>
    /// LocalAssembler は 1 ギャップあたり高々数千リード・数十万 k-mer 程度しか扱わないため、TrustedKmerIndex のシャード分割・外部ソート・ディスクマージは過剰でしかない (ギャップの数だけ一時ディレクトリの作成とファイルの生成・マージ・削除が走り、それ自体が支配的なコストになっていた) <br/>
    /// この規模ならインメモリの HashSet で完結できるため、カウントは持たず「見たことがあるか」だけを覚える
    /// </remarks>
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
        /// 信頼できる k-mer 集合 (k &gt; 64)
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
                : this._中 is { } l_中 ? l_中.Add(TrustedKmerIndex.Get_正規形_中(p_kmer)) : this._大!.Add(new KmerKey(p_kmer).Get_正規形());
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
                : this._中 is { } l_中 ? l_中.Contains(TrustedKmerIndex.Get_正規形_中(p_kmer)) : this._大!.Contains(new KmerKey(p_kmer).Get_正規形());
        }

        #endregion
    }
}
