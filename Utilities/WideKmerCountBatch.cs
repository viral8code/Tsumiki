namespace Tsumiki.Utilities
{
    /// <summary>
    /// ワーカーごとに k &gt; 128 の k-mer のパック済みバイト列をシャード単位で溜め、まとめてカウンタへ渡す
    /// </summary>
    internal sealed class WideKmerCountBatch
    {
        #region 定数

        /// <summary>
        /// 1 シャードぶん溜める件数
        /// </summary>
        private const int C_束の件数 = 4_096;

        #endregion

        #region 内部変数

        /// <summary>
        /// 渡し先の索引
        /// </summary>
        private readonly TrustedKmerIndex _索引;

        /// <summary>
        /// シャードごとに溜めたパック済みバイト列
        /// </summary>
        private readonly List<byte[]>[] _束;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// 索引のシャード数ぶんの束を用意する
        /// </summary>
        /// <param name="p_索引"></param>
        public WideKmerCountBatch(TrustedKmerIndex p_索引)
        {
            this._索引 = p_索引;
            this._束 = [.. Enumerable.Range(0, Math.Max(1, p_索引.A_シャード数)).Select(_ => new List<byte[]>(C_束の件数))];
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// パック済みバイト列を 1 件溜める
        /// </summary>
        /// <param name="p_パック済み"></param>
        public void V_追加(byte[] p_パック済み)
        {
            var l_シャード = TrustedKmerIndex.Get_シャード番号(p_パック済み, this._束.Length);
            this._束[l_シャード].Add(p_パック済み);
            if (this._束[l_シャード].Count == C_束の件数)
            {
                this._索引.V_登録_パック済み群(l_シャード, this._束[l_シャード]);
                this._束[l_シャード].Clear();
            }
        }

        /// <summary>
        /// 溜まっている分をすべて渡す
        /// </summary>
        public void V_吐き出し()
        {
            for (var i = 0; i < this._束.Length; i++)
            {
                if (this._束[i].Count > 0)
                {
                    this._索引.V_登録_パック済み群(i, this._束[i]);
                    this._束[i].Clear();
                }
            }
        }

        #endregion
    }
}
