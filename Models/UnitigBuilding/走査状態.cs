using Tsumiki.Commons;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.Utilities;

namespace Tsumiki.Models.UnitigBuilding
{
    /// <summary>
    /// ワーカーごとに持つ走査用の状態
    /// </summary>
    /// <param name="p_kmerインデックス">信頼できる k-mer 集合</param>
    internal sealed class 走査状態(TrustedKmerIndex p_kmerインデックス)
    {
        #region 内部変数

        /// <summary>
        /// 固定幅キーを転がす高速 walk
        /// </summary>
        private readonly UnitigWalk? _高速walk =
            UnitigWalk.Is対応k長(ConfigurationManager.A_実行時引数.A_k長)
                ? new UnitigWalk(p_kmerインデックス, ConfigurationManager.A_実行時引数.A_k長)
                : null;

        /// <summary>
        /// 固定幅キーを使えない k 長向けの参照実装
        /// </summary>
        private readonly UnitigMaker _参照実装 = new(p_kmerインデックス);

        /// <summary>
        /// 訪問済み
        /// </summary>
        private readonly HashSet<UInt128> _訪問済み = [];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 開始 k-mer から walk して配列を返す
        /// </summary>
        /// <param name="p_開始kmer">walk を始める k-mer</param>
        /// <returns>組み上がった配列</returns>
        public string Get_配列(byte[] p_開始kmer)
        {
            if (this._高速walk is not { } l_高速walk)
            {
                return this._参照実装.Get_ユニティグ(p_開始kmer).A_配列;
            }
            var l_塩基列 = l_高速walk.Get_塩基列(p_開始kmer, this._訪問済み);
            return string.Create(l_塩基列.Count, l_塩基列,
                static (l_文字, l_元) =>
                {
                    for (var i = 0; i < l_元.Count; i++)
                    {
                        l_文字[i] = Util.Get_塩基文字(l_元[i]);
                    }
                });
        }

        #endregion
    }
}
