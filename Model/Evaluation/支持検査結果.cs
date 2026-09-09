namespace Tsumiki.Model.Evaluation
{
    /// <summary>
    /// 出した配列がリードに裏付けられているかの検査結果
    /// </summary>
    internal readonly record struct 支持検査結果(int A_r長, long A_調べた位置数, long A_支持のない位置数, IReadOnlyList<支持のない区間> A_区間)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 調べた位置のうち、どのリードにも現れなかった割合
        /// </summary>
        public double A_支持のない率 => this.A_調べた位置数 == 0L ? 0D : 100D * this.A_支持のない位置数 / this.A_調べた位置数;

        #endregion
    }
}
