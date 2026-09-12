using Tsumiki.Commons;

namespace Tsumiki.Models.Foundation
{
    /// <summary>
    /// k-mer をパックした固定長キー
    /// </summary>
    internal readonly struct KmerKey : IEquatable<KmerKey>
    {
        #region 内部変数

        /// <summary>
        /// 作成時の塩基数
        /// </summary>
        private readonly int _長さ;

        /// <summary>
        /// パック済みデータ
        /// </summary>
        public readonly ulong[] A_パック済みデータ;

        #endregion

        #region コンストラクタ

        /// <summary>
        /// char の k-mer 文字列からパック済みデータを構築する
        /// </summary>
        /// <param name="p_kmer">パックする k-mer 文字列</param>
        public KmerKey(ReadOnlySpan<char> p_kmer)
        {
            this._長さ = p_kmer.Length;
            this.A_パック済みデータ = new ulong[(p_kmer.Length + 31) >> 5];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                var l_要素位置 = i >> 5;
                var l_シフト量 = (31 ^ (i & 31)) << 1;
                // Get_塩基 ID 候補 は曖昧塩基対応のため List を確保するが、
                // ContigMaker 側では曖昧塩基を含む区間はそもそも KmerKey 化されない
                // (呼ばれない) ため、ここでは List 確保のない軽量な単一塩基変換で十分
                var l_値 = Util.Get_塩基ID(p_kmer[i]) - 1UL;
                // 32 塩基ごとに同じ ulong 要素 (2 bit x 32 = 64 bit) を共有するため、
                // 代入ではなく OR で詰め込まないと、直前までに書き込んだ
                // 塩基の情報が上書きで消えてしまう
                // (この不具合により、同じ ulong 要素に収まる k-mer 同士が
                // 実質「末尾の数文字だけで同一視される」形になっていた)
                this.A_パック済みデータ[l_要素位置] |= l_値 << l_シフト量;
            }
        }

        /// <summary>
        /// 塩基 ID (1=A,2=C,3=G,4=T) のバイト列から直接構築する版
        /// </summary>
        /// <param name="p_kmer">パックする塩基 ID 列</param>
        /// <remarks>
        /// UnitigMaker/TrustedKmerIndex はバイト ID 空間で動作しているため、char 経由の変換を挟まずに済む (ホットパス向け)
        /// </remarks>
        public KmerKey(ReadOnlySpan<byte> p_kmer)
        {
            this._長さ = p_kmer.Length;
            this.A_パック済みデータ = new ulong[(p_kmer.Length + 31) >> 5];
            for (var i = 0; i < p_kmer.Length; i++)
            {
                var l_要素位置 = i >> 5;
                var l_シフト量 = (31 ^ (i & 31)) << 1;
                var l_値 = p_kmer[i] - 1UL;
                this.A_パック済みデータ[l_要素位置] |= l_値 << l_シフト量;
            }
        }

        #endregion

        #region 公開メソッド

        /// <summary>
        /// この k-mer とその逆相補のうち、パック済みデータを辞書式順序で比較して小さい方を返す
        /// </summary>
        /// <remarks>
        /// 挿入時・検索時の双方でこれを使えば、順鎖/逆鎖どちらから見ても同一のキーに正規化されるため、逆相補を別途リトライする必要がなくなる
        /// </remarks>
        /// <returns></returns>
        public KmerKey Get_正規形()
        {
            var l_逆相補 = this.Get_逆相補();
            return Get_比較結果(this.A_パック済みデータ, l_逆相補.A_パック済みデータ) <= 0 ? this : l_逆相補;
        }

        /// <summary>
        /// 塩基 ID 列へデコードしてから逆相補を取り、再エンコードする
        /// </summary>
        /// <remarks>
        /// 64 bit 全体のビット反転で済ませてはいけない<br/>
        /// 2 bit コドン内部のビット順まで入れ替わり、C (01) と G (10) のような塩基で値が化ける
        /// </remarks>
        /// <returns></returns>
        public KmerKey Get_逆相補()
        {
            var l_逆相補 = Util.V_逆相補(this.Get_塩基列(this._長さ).AsSpan());
            return new KmerKey(l_逆相補);
        }

        /// <summary>
        /// パック済みデータを、塩基 ID (1=A,2=C,3=G,4=T) のバイト列へデコードする
        /// </summary>
        /// <param name="p_長さ">元の k-mer 長 (コンストラクタに渡した長さ) </param>
        /// <returns></returns>
        public byte[] Get_塩基列(int p_長さ)
        {
            if (p_長さ < 0 || p_長さ > this._長さ)
            {
                throw new ArgumentOutOfRangeException(nameof(p_長さ));
            }
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

        /// <summary>
        /// 同じ k-mer か
        /// </summary>
        /// <param name="p_比較対象">比べる k-mer</param>
        /// <returns>同じなら true</returns>
        public bool Equals(KmerKey p_比較対象)
        {
            return this._長さ == p_比較対象._長さ && this.A_パック済みデータ.AsSpan().SequenceEqual(p_比較対象.A_パック済みデータ);
        }

        /// <summary>
        /// (オーバーライド) 同じ k-mer か
        /// </summary>
        /// <param name="p_対象">比べる対象</param>
        /// <returns></returns>
        public override bool Equals(object? p_対象)
        {
            return p_対象 is KmerKey l_比較対象 && this.Equals(l_比較対象);
        }

        /// <summary>
        /// (オーバーライド) ハッシュ値を返す
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            var l_ハッシュ = 1469598103934665603UL ^ (ulong)this._長さ;
            foreach (var l_要素 in this.A_パック済みデータ.AsSpan())
            {
                l_ハッシュ ^= l_要素;
                l_ハッシュ *= 1099511628211UL;
            }
            return (int)(l_ハッシュ ^ (l_ハッシュ >> 32));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// パック済みデータを辞書式順序で比べる
        /// </summary>
        /// <param name="p_左">比べるパック済みデータ</param>
        /// <param name="p_右">比べるパック済みデータ</param>
        /// <returns>左が小さければ -1、大きければ 1、等しければ 0</returns>
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

        #endregion
    }
}
