using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Utilities;

namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// ワーカーごとに持つ走査用の状態
    /// </summary>
    /// <remarks>
    /// k &lt;= 64 なら転がし更新の実装を使い、
    /// それを超える場合だけ従来の実装へ落ちる
    /// </remarks>
    /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
    internal sealed class 走査状態(TrustedKmerIndex p_kmerインデックス)
    {
        /// <summary>
        /// 転がし
        /// </summary>
        private readonly UnitigWalk? _転がし =
            UnitigWalk.Get_扱えるか(ConfigurationManager.A_実行時引数.A_k長)
                ? new UnitigWalk(p_kmerインデックス, ConfigurationManager.A_実行時引数.A_k長)
                : null;

        /// <summary>
        /// 従来
        /// </summary>
        private readonly UnitigMaker _従来 = new(p_kmerインデックス);

        /// <summary>
        /// 訪問済み
        /// </summary>
        private readonly HashSet<UInt128> _訪問済み = [];

        /// <summary>
        /// 開始 k-mer から walk して配列を返す
        /// </summary>
        /// <param name="p_開始kmer">walk を始める k-mer</param>
        /// <returns>組み上がった配列</returns>
        public string Get_配列(byte[] p_開始kmer)
        {
            if (this._転がし is not { } l_転がし)
            {
                return this._従来.Get_ユニティグ(p_開始kmer).A_配列;
            }
            var l_塩基列 = l_転がし.Get_塩基列(p_開始kmer, this._訪問済み);
            return string.Create(l_塩基列.Count, l_塩基列,
                static (l_文字, l_元) =>
                {
                    for (var i = 0; i < l_元.Count; i++)
                    {
                        l_文字[i] = Util.Get_塩基文字(l_元[i]);
                    }
                });
        }
    }
}
