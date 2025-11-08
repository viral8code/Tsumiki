namespace Tsumiki.Model
{
    internal readonly struct Sequence(string id, string seq)
    {
        public readonly string ID = id;
        public readonly string Seq = seq;
    }
}
