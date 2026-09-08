namespace Tsumiki.Model
{
    /// <summary>
    /// アセンブリ結果(unitigs/contigs/scaffolds)の品質を大まかに把握するための
    /// 基本統計量。リファレンスなしで計算できる範囲の指標のみを対象とする。
    /// </summary>
    internal readonly struct アセンブリ統計(
        int p_配列数,
        long p_総延長,
        int p_最大長,
        int p_最小長,
        int p_N50,
        int p_L50,
        double p_GC率)
    {
        public readonly int A_配列数 = p_配列数;
        public readonly long A_総延長 = p_総延長;
        public readonly int A_最大長 = p_最大長;
        public readonly int A_最小長 = p_最小長;

        /// <summary>長い順に並べて累積長が全長の50%に達した時点の配列長。</summary>
        public readonly int A_N50 = p_N50;

        /// <summary>N50 に達するまでに必要だった配列の本数(1始まり)。</summary>
        public readonly int A_L50 = p_L50;

        public readonly double A_GC率 = p_GC率;

        public override string ToString()
        {
            return $"count={this.A_配列数}, total_length={this.A_総延長}, " +
                   $"N50={this.A_N50}, L50={this.A_L50}, max={this.A_最大長}, min={this.A_最小長}, " +
                   $"GC%={this.A_GC率:0.00}";
        }
    }
}
