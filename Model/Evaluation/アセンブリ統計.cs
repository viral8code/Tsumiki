namespace Tsumiki.Model.Evaluation
{
    /// <summary>
    /// アセンブリ結果 (unitigs/contigs/scaffolds) の品質を大まかに把握するための基本統計量
    /// </summary>
    /// <remarks>
    /// リファレンスなしで計算できる範囲の指標のみを対象とする
    /// </remarks>
    internal readonly record struct アセンブリ統計(int A_配列数, long A_総延長, int A_最大長, int A_最小長, int A_N50, int A_L50, double A_GC率)
    {
        #region 継承メソッド

        public override string ToString()
        {
            return $"count={this.A_配列数}, total_length={this.A_総延長}, " +
                   $"N50={this.A_N50}, L50={this.A_L50}, max={this.A_最大長}, min={this.A_最小長}, " +
                   $"GC%={this.A_GC率:0.00}";
        }

        #endregion
    }
}
