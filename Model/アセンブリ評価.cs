namespace Tsumiki.Model
{
    /// <summary>
    /// リファレンス無しで測ったアセンブリの良さ。
    /// 「この配列が正しいなら各 k-mer は何回現れるはずか」をカバレッジから求め、
    /// 実際の出力と突き合わせた結果。
    /// </summary>
    internal record アセンブリ評価(
        long A_期待延べ数,
        long A_欠損延べ数,
        long A_過剰延べ数,
        long A_総延長,
        int A_本数,
        long A_NG50)
    {
        /// <summary>
        /// 出すべき k-mer のうち実際に出せた割合。
        /// 反復配列を飛ばして繋いだ誤アセンブリは、飛ばした領域の k-mer が
        /// 欠損として現れるためここに反映される。
        /// </summary>
        public double A_完全性 => this.A_期待延べ数 == 0
            ? 0
            : 1.0 - ((double)this.A_欠損延べ数 / this.A_期待延べ数);

        /// <summary>
        /// カバレッジが支持する以上に同じ配列を出していない度合い。
        /// </summary>
        public double A_正確性 => this.A_期待延べ数 == 0
            ? 0
            : 1.0 - ((double)this.A_過剰延べ数 / this.A_期待延べ数);

        public override string ToString()
        {
            return $"NG50={this.A_NG50:N0}, completeness={this.A_完全性 * 100:F2}%, " +
                $"accuracy={this.A_正確性 * 100:F2}% " +
                $"({this.A_本数} seq(s), {this.A_総延長:N0} bp)";
        }
    }
}
