namespace Tsumiki.Utilities
{
    /// <summary>
    /// 1 ワーカーが見つけた r-mer を分割ごとに溜め、溜まったらまとめて登録する
    /// </summary>
    /// <param name="p_検証器">登録先</param>
    internal sealed class 登録束(RepeatRMerVerifier p_検証器)
    {
        #region 定数

        /// <summary>
        /// 1 分割に溜める件数
        /// </summary>
        private const int C_束の件数 = 1_024;

        #endregion

        #region 内部変数

        /// <summary>
        /// 分割ごとの溜め置き
        /// </summary>
        private readonly (UInt128 A_上位, UInt128 A_下位)[]?[] _束 = new (UInt128 A_上位, UInt128 A_下位)[]?[RepeatRMerVerifier.C_分割数];

        /// <summary>
        /// 分割ごとの溜めた件数
        /// </summary>
        private readonly int[] _件数 = new int[RepeatRMerVerifier.C_分割数];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// r-mer を 1 件溜める
        /// </summary>
        /// <param name="p_値">正準化した r-mer の値</param>
        public void V_追加((UInt128 A_上位, UInt128 A_下位) p_値)
        {
            var l_分割 = RepeatRMerVerifier.Get_分割番号(p_値);
            var l_束 = this._束[l_分割] ??= new (UInt128 A_上位, UInt128 A_下位)[C_束の件数];
            l_束[this._件数[l_分割]++] = p_値;
            if (this._件数[l_分割] == C_束の件数)
            {
                p_検証器.V_登録_束(l_分割, l_束);
                this._件数[l_分割] = 0;
            }
        }

        /// <summary>
        /// 溜めた分をすべて登録する
        /// </summary>
        public void V_吐き出し()
        {
            for (var l_分割 = 0; l_分割 < RepeatRMerVerifier.C_分割数; l_分割++)
            {
                if (this._件数[l_分割] > 0)
                {
                    p_検証器.V_登録_束(l_分割, this._束[l_分割].AsSpan(0, this._件数[l_分割]));
                    this._件数[l_分割] = 0;
                }
            }
        }

        #endregion
    }
}
