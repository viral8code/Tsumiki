namespace Tsumiki.Models.Evaluation
{
    /// <summary>
    /// リファレンス無しで測ったアセンブリの良さ
    /// </summary>
    /// <remarks>
    /// 「この配列が正しいなら各 k-mer は何回現れるはずか」をカバレッジから求め、実際の出力と突き合わせた結果
    /// </remarks>
    /// <param name="A_期待延べ数"></param>
    /// <param name="A_欠損延べ数"></param>
    /// <param name="A_過剰延べ数"></param>
    /// <param name="A_総延長"></param>
    /// <param name="A_本数"></param>
    /// <param name="A_NG50"></param>
    /// <param name="A_環状本数"></param>
    /// <param name="A_環状化率"></param>
    internal record アセンブリ評価(long A_期待延べ数, long A_欠損延べ数, long A_過剰延べ数, long A_総延長, int A_本数, long A_NG50, int A_環状本数 = 0, double A_環状化率 = 0D)
    {
        #region カスタムプロパティ

        /// <summary>
        /// 出すべき k-mer のうち実際に出せた割合
        /// </summary>
        /// <remarks>
        /// 反復配列を飛ばして繋いだ誤アセンブリは、飛ばした領域の k-mer が欠損として現れるためここに反映される
        /// </remarks>
        public double A_完全性 => this.A_期待延べ数 == 0L ? 0D : 1D - ((double)this.A_欠損延べ数 / this.A_期待延べ数);

        /// <summary>
        /// カバレッジが支持する以上に同じ配列を出していない度合い
        /// </summary>
        public double A_正確性 => this.A_期待延べ数 == 0L ? 0D : 1D - ((double)this.A_過剰延べ数 / this.A_期待延べ数);

        #endregion

        #region 継承メソッド

        /// <summary>
        /// (オーバーライド) アセンブリ評価を表す文字列を返す
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return $"circular={this.A_環状本数} ({this.A_環状化率 * 100D:F1}% of genome), " +
                $"NG50={this.A_NG50:N0}, completeness={this.A_完全性 * 100D:F2}%, " +
                $"accuracy={this.A_正確性 * 100D:F2}% " +
                $"({this.A_本数} seq(s), {this.A_総延長:N0} bp)";
        }

        #endregion
    }
}
