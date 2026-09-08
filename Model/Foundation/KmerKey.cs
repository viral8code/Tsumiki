using Tsumiki.Common;

namespace Tsumiki.Model
{
    internal readonly struct KmerKey : IEquatable<KmerKey>
    {
        public readonly ulong[] A_パック済みデータ;

        public KmerKey(ReadOnlySpan<char> p_kmer)
        {
            this.A_パック済みデータ = new ulong[(p_kmer.Length + 31) >> 5];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                var l_要素位置 = i >> 5;
                var l_シフト量 = (31 ^ (i & 31)) << 1;
                // Get_塩基ID候補 は曖昧塩基対応のため List を確保するが、
                // ContigMaker 側では曖昧塩基を含む区間はそもそも KmerKey 化されない
                // (呼ばれない)ため、ここでは List 確保のない軽量な単一塩基変換で十分。
                var l_値 = (ulong)Util.Get_塩基ID(p_kmer[i]) - 1;
                // 32塩基ごとに同じ ulong 要素(2bit x 32 = 64bit)を共有するため、
                // 代入ではなく OR で詰め込まないと、直前までに書き込んだ
                // 塩基の情報が上書きで消えてしまう。
                // (この不具合により、同じ ulong 要素に収まる k-mer 同士が
                //  実質「末尾の数文字だけで同一視される」形になっていた。)
                this.A_パック済みデータ[l_要素位置] |= l_値 << l_シフト量;
            }
        }

        /// <summary>
        /// 塩基ID(1=A,2=C,3=G,4=T)のバイト列から直接構築する版。
        /// UnitigMaker/TrustedKmerIndex はバイトID空間で動作しているため、
        /// char経由の変換を挟まずに済む(ホットパス向け)。
        /// </summary>
        public KmerKey(ReadOnlySpan<byte> p_kmer)
        {
            this.A_パック済みデータ = new ulong[(p_kmer.Length + 31) >> 5];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                var l_要素位置 = i >> 5;
                var l_シフト量 = (31 ^ (i & 31)) << 1;
                var l_値 = (ulong)p_kmer[i] - 1;
                this.A_パック済みデータ[l_要素位置] |= l_値 << l_シフト量;
            }
        }

        private KmerKey(ulong[] p_パック済みデータ)
        {
            this.A_パック済みデータ = p_パック済みデータ;
        }

        /// <summary>
        /// この k-mer とその逆相補のうち、パック済みデータを辞書式順序で比較して
        /// 小さい方を返す。挿入時・検索時の双方でこれを使えば、順鎖/逆鎖どちらから
        /// 見ても同一のキーに正規化されるため、逆相補を別途リトライする必要がなくなる。
        /// </summary>
        public KmerKey Get_正規形()
        {
            var l_逆相補 = this.Get_逆相補();
            return Get_比較結果(this.A_パック済みデータ, l_逆相補.A_パック済みデータ) <= 0 ? this : l_逆相補;
        }

        private static int Get_比較結果(ulong[] p_左, ulong[] p_右)
        {
            for (var i = 0; i < p_左.Length; i++)
            {
                if (p_左[i] != p_右[i])
                {
                    return p_左[i] < p_右[i] ? -1 : 1;
                }
            }
            return 0;
        }

        /// <summary>
        /// 塩基ID列へデコードしてから逆相補を取り、再エンコードする。
        /// 64bit 全体のビット反転で済ませてはいけない。2bit コドン内部の
        /// ビット順まで入れ替わり、C(01) と G(10) のような塩基で値が化ける。
        /// </summary>
        public KmerKey Get_逆相補()
        {
            var l_逆相補 = Util.V_逆相補(this.Get_塩基列(ConfigurationManager.A_実行時引数.A_k長).AsSpan());
            return new KmerKey(l_逆相補);
        }

        /// <summary>
        /// パック済みデータを、塩基ID(1=A,2=C,3=G,4=T)のバイト列へデコードする。
        /// p_長さ は元のk-mer長(コンストラクタに渡した長さ)を指定する。
        /// </summary>
        public byte[] Get_塩基列(int p_長さ)
        {
            var l_塩基列 = new byte[p_長さ];
            for (var i = 0; i < p_長さ; i++)
            {
                var l_要素位置 = i >> 5;
                var l_シフト量 = (31 ^ (i & 31)) << 1;
                var l_値 = (byte)((this.A_パック済みデータ[l_要素位置] >> l_シフト量) & 0x3UL);
                l_塩基列[i] = (byte)(l_値 + 1);
            }
            return l_塩基列;
        }

        public bool Equals(KmerKey other)
        {
            if (this.A_パック済みデータ.Length != other.A_パック済みデータ.Length)
            {
                return false;
            }

            for (var i = 0; i < this.A_パック済みデータ.Length; i++)
            {
                if (this.A_パック済みデータ[i] != other.A_パック済みデータ[i])
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
            var l_ハッシュ = 1469598103934665603UL;
            foreach (var l_要素 in this.A_パック済みデータ)
            {
                l_ハッシュ ^= l_要素;
                l_ハッシュ *= 1099511628211UL;
            }
            return (int)(l_ハッシュ ^ (l_ハッシュ >> 32));
        }
    }
}
