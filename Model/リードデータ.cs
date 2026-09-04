namespace Tsumiki.Model
{
    internal class リードデータ
    {
        public required string A_ID { get; set; }

        /// <summary>
        /// 曖昧塩基を許容する経路(Program.V_読込_リードファイル_曖昧塩基あり)向け。
        /// 各塩基が取りうる ID の候補リスト。
        /// </summary>
        public List<byte[]>? A_塩基候補列 { get; set; }

        /// <summary>
        /// 曖昧塩基を無視する経路(KmerCounting.V_読込_リードファイル)向けの軽量表現。
        /// 各塩基を1バイトの ID に変換したもの。A/C/G/T 以外は Consts.無効な塩基。
        /// LINQ や per-base の配列アロケーションを避けるための専用フィールド。
        /// </summary>
        public byte[]? A_塩基列 { get; set; }

        public required string A_生リード { get; set; }

        public required string A_クオリティ { get; set; }

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
