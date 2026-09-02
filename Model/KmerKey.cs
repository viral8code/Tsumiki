using Tsumiki.Common;

namespace Tsumiki.Model
{
    internal readonly struct KmerKey : IEquatable<KmerKey>
    {
        public readonly ulong[] Data;

        public KmerKey(ReadOnlySpan<char> kmer)
        {
            this.Data = new ulong[(kmer.Length + 31) >> 5];
            for (var i = 0; i < kmer.Length; i++)
            {
                var index = i >> 5;
                var shift = (31 ^ (i & 31)) << 1;
                // GetNucleotideIDs は曖昧塩基対応のため List<int> を確保するが、
                // ContigMaker 側では badBase/revBadBase によって曖昧塩基を含む
                // 区間はそもそも KmerKey 化されない(呼ばれない)ため、
                // ここでは List 確保のない軽量な単一塩基変換で十分。
                // (曖昧塩基が来た場合は GetSimpleNucleotideID が InvalidBase(5) を
                //  返すが、そのようなケースは呼び出し元で事前に除外されている前提。)
                var val = (ulong)Util.GetSimpleNucleotideID(kmer[i]) - 1;
                // 32塩基ごとに同じ ulong 要素(2bit x 32 = 64bit)を共有するため、
                // 代入(=)ではなく OR(|=)で詰め込まないと、直前までに書き込んだ
                // 塩基の情報が上書きで消えてしまう。
                // (この不具合により、同じ ulong 要素に収まる k-mer 同士が
                //  実質「末尾の数文字だけで同一視される」形になっていた。)
                this.Data[index] |= val << shift;
            }
        }

        /// <summary>
        /// Consts.NucleotideID(1=A,2=C,3=G,4=T)のバイト列から直接構築する版。
        /// UnitigMaker/TrustedKmerIndex はbyte-ID空間で動作しているため、
        /// char経由の変換を挟まずに済む(ホットパス向け)。
        /// </summary>
        public KmerKey(ReadOnlySpan<byte> kmer)
        {
            this.Data = new ulong[(kmer.Length + 31) >> 5];
            for (var i = 0; i < kmer.Length; i++)
            {
                var index = i >> 5;
                var shift = (31 ^ (i & 31)) << 1;
                var val = (ulong)kmer[i] - 1;
                this.Data[index] |= val << shift;
            }
        }

        private KmerKey(ulong[] Data)
        {
            this.Data = Data;
        }

        /// <summary>
        /// この k-mer とその逆相補のうち、Data を辞書式順序で比較して小さい方を返す。
        /// 挿入時・検索時の双方でこれを使えば、順鎖/逆鎖どちらから見ても
        /// 同一のキーに正規化されるため、逆相補を別途リトライする必要がなくなる。
        /// </summary>
        public KmerKey Canonical()
        {
            var rev = this.ReverseComprement();
            return CompareData(this.Data, rev.Data) <= 0 ? this : rev;
        }

        private static int CompareData(ulong[] a, ulong[] b)
        {
            for (var i = 0; i < a.Length; i++)
            {
                if (a[i] != b[i])
                {
                    return a[i] < b[i] ? -1 : 1;
                }
            }
            return 0;
        }

        /// <summary>
        /// 以前はビット単位の反転("ReverseBit")とNOTの組み合わせで逆相補を
        /// 計算していたが、64bit全体を単純にビット反転すると各2bitコドン
        /// (塩基1個分)の内部のビット順まで入れ替わってしまい
        /// (例: コドン順序は正しく反転されるが、C(01)とG(10)のような
        /// 「2bit内の上位/下位」を持つ塩基同士で値が化けていた)、
        /// 実際には正しい逆相補になっていなかった(k=4,31,33,64のいずれでも
        /// Util.ReverseComprement(string)の結果と一致しないことをテストで確認)。
        /// このバグは Core/ContigMaker.cs の kmerDict 構築(逆鎖k-merの登録)で
        /// 使われており、逆鎖側の読み取りマッピングを広範囲で壊していた。
        ///
        /// 塩基ID列へいったんデコードし、既に実績のある
        /// Util.ReverseComprement(Span&lt;byte&gt;) で逆相補を取ってから
        /// 再エンコードすることで、確実に正しい結果にする。
        /// </summary>
        public KmerKey ReverseComprement()
        {
            var kmerLength = ConfigurationManager.Arguments.Kmer;
            var bytes = new byte[kmerLength];
            for (var i = 0; i < kmerLength; i++)
            {
                var index = i >> 5;
                var shift = (31 ^ (i & 31)) << 1;
                var val = (byte)((this.Data[index] >> shift) & 0x3UL);
                bytes[i] = (byte)(val + 1);
            }
            var revBytes = Util.ReverseComprement(bytes.AsSpan());
            return new KmerKey(revBytes);
        }

        public bool Equals(KmerKey other)
        {
            if (this.Data.Length != other.Data.Length)
            {
                return false;
            }

            for (var i = 0; i < this.Data.Length; i++)
            {
                if (this.Data[i] != other.Data[i])
                {
                    return false;
                }
            }
            return true;
        }

        public override bool Equals(object? obj)
        {
            return obj is KmerKey other && this.Equals(other);
        }

        public override int GetHashCode()
        {
            var hash = 1469598103934665603UL;
            foreach (var v in this.Data)
            {
                hash ^= v;
                hash *= 1099511628211UL;
            }
            return (int)(hash ^ (hash >> 32));
        }
    }
}