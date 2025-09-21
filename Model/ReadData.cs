namespace Tsumiki.Model
{
    internal class ReadData
    {
        public required string ID { get; set; }

        public required List<byte[]> Read { get; set; }

        public required string Quality { get; set; }

        public override string ToString()
        {
            return $@"
                ID      : {ID}
                read    : {Read}
                quality : {Quality}";
        }
    }
}
