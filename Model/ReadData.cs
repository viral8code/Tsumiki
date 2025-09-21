namespace Tsumiki.Model
{
    internal class ReadData
    {
        public required string ID { get; set; }

        public required List<byte[]> Read { get; set; }

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
