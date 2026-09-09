namespace Tsumiki.Model.Foundation
{
    /// <summary>
    /// FASTQ の 1 リード
    /// </summary>
    internal class リードデータ
    {
        /// <summary>
        /// リード ID
        /// </summary>
        public required string A_ID { get; set; }

        /// <summary>
        /// 曖昧塩基を許容する経路 (Program.V_読込_リードファイル_曖昧塩基あり) 向け
        /// </summary>
        /// <remarks>
        /// 各塩基が取りうる ID の候補リスト
        /// </remarks>
        public List<byte[]>? A_塩基候補列 { get; set; }

        /// <summary>
        /// 曖昧塩基を無視する経路 (KmerCounting.V_読込_リードファイル) 向けの軽量表現
        /// </summary>
        /// <remarks>
        /// 各塩基を 1 バイトの ID に変換したもの<br/>
        /// A/C/G/T 以外は Consts.無効な塩基<br/>
        /// LINQ や per-base の配列アロケーションを避けるための専用フィールド
        /// </remarks>
        public byte[]? A_塩基列 { get; set; }

        /// <summary>
        /// 生リード
        /// </summary>
        public required string A_生リード { get; set; }

        /// <summary>
        /// クオリティ
        /// </summary>
        public required string A_クオリティ { get; set; }

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
    }
}
