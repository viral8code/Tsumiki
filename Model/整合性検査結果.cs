namespace Tsumiki.Model
{
    /// <summary>
    /// アセンブリが観測された k-mer とその出現回数に対して辻褄が合っているかの検査結果。
    /// </summary>
    internal readonly record struct 整合性検査結果(
        long A_信頼kmer数,
        long A_アセンブリ内の延べ数,
        long A_アセンブリ内の種類数,
        long A_取りこぼし数,
        long A_出しすぎkmer種類数,
        long A_余分な延べ数)
    {
        /// <summary>信頼できる k-mer のうち、アセンブリに現れなかった割合。</summary>
        public double A_取りこぼし率 => this.A_信頼kmer数 == 0 ? 0 : 100.0 * this.A_取りこぼし数 / this.A_信頼kmer数;

        /// <summary>
        /// アセンブリ中の k-mer 延べ数のうち、コピー数の推定を超えて
        /// 余分に現れている分の割合。総延長の水増し量にほぼ対応する。
        /// </summary>
        public double A_出しすぎ率 => this.A_アセンブリ内の延べ数 == 0 ? 0 : 100.0 * this.A_余分な延べ数 / this.A_アセンブリ内の延べ数;
    }
}
