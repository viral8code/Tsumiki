namespace Tsumiki.Utilities
{
    /// <summary>
    /// ワーカーごとに k-mer のパック値をシャード単位で溜め、まとめてカウンタへ渡す
    /// </summary>
    /// <remarks>
    /// k-mer ごとにシャードの錠を取ると、全リードの全 k-mer の回数だけワーカー同士が錠を取り合う
    /// </remarks>
    internal sealed class KmerCountBatch
    {
        #region 定数

        /// <summary>
        /// 1 シャードぶん溜める件数
        /// </summary>
        private const int 束の件数 = 4_096;

        #endregion

        #region 内部変数

        /// <summary>
        /// 渡し先の索引
        /// </summary>
        private readonly TrustedKmerIndex _索引;

        /// <summary>
        /// シャードごとに溜めたパック値
        /// </summary>
        private readonly (UInt128 A_上位, UInt128 A_下位)[][] _束;

        /// <summary>
        /// シャードごとの溜まっている件数
        /// </summary>
        private readonly int[] _件数;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 索引のシャード数ぶんの束を用意する
        /// </summary>
        /// <param name="p_索引"></param>
        public KmerCountBatch(TrustedKmerIndex p_索引)
        {
            this._索引 = p_索引;
            var l_シャード数 = Math.Max(1, p_索引.A_シャード数);
            this._束 = new (UInt128 A_上位, UInt128 A_下位)[l_シャード数][];
            for (var i = 0; i < l_シャード数; i++)
            {
                this._束[i] = new (UInt128 A_上位, UInt128 A_下位)[束の件数];
            }
            this._件数 = new int[l_シャード数];
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 正規形のパック値を 1 件溜める
        /// </summary>
        /// <param name="p_上位"></param>
        /// <param name="p_下位"></param>
        public void V_追加(UInt128 p_上位, UInt128 p_下位)
        {
            var l_シャード = TrustedKmerIndex.Get_シャード番号((p_上位, p_下位), this._束.Length);
            this._束[l_シャード][this._件数[l_シャード]++] = (p_上位, p_下位);
            if (this._件数[l_シャード] == 束の件数)
            {
                this._索引.V_登録_値群(l_シャード, this._束[l_シャード]);
                this._件数[l_シャード] = 0;
            }
        }

        /// <summary>
        /// 溜まっている分をすべて渡す
        /// </summary>
        public void V_吐き出し()
        {
            for (var i = 0; i < this._束.Length; i++)
            {
                if (this._件数[i] > 0)
                {
                    this._索引.V_登録_値群(i, this._束[i].AsSpan(0, this._件数[i]));
                    this._件数[i] = 0;
                }
            }
        }

        #endregion
    }
}
