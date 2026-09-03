namespace Tsumiki.Model
{
    internal class ReadData
    {
        public required string ID { get; set; }

        /// <summary>
        /// 曖昧塩基を許容する経路(Program.LoadReadFileWithAmbiguousBases)向け。
        /// 各塩基が取りうる ID の候補リスト。
        /// </summary>
        public List<byte[]>? Read { get; set; }

        /// <summary>
        /// 曖昧塩基を無視する経路(KmerCounting.LoadReadFile)向けの軽量表現。
        /// 各塩基を1バイトの ID に変換したもの。A/C/G/T 以外は Consts.InvalidBase。
        /// LINQ や per-base の配列アロケーションを避けるための専用フィールド。
        /// </summary>
        public byte[]? SimpleRead { get; set; }

        public required string RowRead { get; set; }

        public required string Quality { get; set; }

        public override string ToString()
        {
            return $"""

                ID      : {this.ID}
                read    : {this.Read}
                quality : {this.Quality}

                """;
        }
    }
}