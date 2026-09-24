namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// FASTQ の 1 リード
    /// </summary>
    internal class リードデータ
    {
        #region プロパティ

        /// <summary>
        /// リード ID
        /// </summary>
        public required string A_ID { get; set; }

        /// <summary>
        /// 曖昧塩基を許容する経路 (Program.V_読込_リードファイル_曖昧塩基あり) 向け
        /// </summary>
        public List<byte[]>? A_塩基候補列 { get; set; }

        /// <summary>
        /// 曖昧塩基を無視する経路 (KmerCounting.V_読込_リードファイル) 向けの軽量表現
        /// </summary>
        public byte[]? A_塩基列 { get; set; }

        /// <summary>
        /// 生リード
        /// </summary>
        public required string A_生リード { get; set; }

        /// <summary>
        /// クオリティ
        /// </summary>
        public required string A_クオリティ { get; set; }

        #endregion

        #region 継承メソッド

        /// <summary>
        /// (オーバーライド) リードを FASTQ の 4 行として返す
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return $"""

                ID      : {this.A_ID}
                read    : {this.A_塩基候補列}
                quality : {this.A_クオリティ}

                """;
        }

        #endregion
    }
}
