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
        /// 先頭の 32 塩基
        /// </summary>
        private readonly ulong _先頭語;

        /// <summary>
        /// 続く 32 塩基
        /// </summary>
        private readonly ulong _第2語;

        /// <summary>
        /// 64 塩基を超える場合のパック済みデータ
        /// </summary>
        private readonly ulong[]? _長いパック済みデータ;

        #endregion

        #region プロパティ

        /// <summary>
        /// パック済みデータ
        /// </summary>
        public ulong[] A_パック済みデータ => this._長いパック済みデータ ?? (this._長さ <= 32 ? [this._先頭語] : [this._先頭語, this._第2語]);

        #endregion

        #region コンストラクタ

        /// <summary>
        /// char の k-mer 文字列からパック済みデータを構築する
        /// </summary>
        /// <param name="p_kmer">パックする k-mer 文字列</param>
        public KmerKey(ReadOnlySpan<char> p_kmer)
        {
            this._長さ = p_kmer.Length;
            this._先頭語 = 0UL;
            this._第2語 = 0UL;
            this._長いパック済みデータ = p_kmer.Length > 64 ? new ulong[(p_kmer.Length + 31) >> 5] : null;
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
                if (this._長いパック済みデータ is not null)
                {
                    this._長いパック済みデータ[l_要素位置] |= l_値 << l_シフト量;
                }
                else if (l_要素位置 == 0)
                {
                    this._先頭語 |= l_値 << l_シフト量;
                }
                else
                {
                    this._第2語 |= l_値 << l_シフト量;
                }
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
            this._先頭語 = 0UL;
            this._第2語 = 0UL;
            this._長いパック済みデータ = p_kmer.Length > 64 ? new ulong[(p_kmer.Length + 31) >> 5] : null;
            for (var i = 0; i < p_kmer.Length; i++)
            {
                var l_要素位置 = i >> 5;
                var l_シフト量 = (31 ^ (i & 31)) << 1;
                var l_値 = p_kmer[i] - 1UL;
                if (this._長いパック済みデータ is not null)
                {
                    this._長いパック済みデータ[l_要素位置] |= l_値 << l_シフト量;
                }
                else if (l_要素位置 == 0)
                {
                    this._先頭語 |= l_値 << l_シフト量;
                }
                else
                {
                    this._第2語 |= l_値 << l_シフト量;
                }
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
            return Get_比較結果(this, l_逆相補) <= 0 ? this : l_逆相補;
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
            var l_逆相補 = this._長さ <= 128 ? stackalloc byte[this._長さ] : new byte[this._長さ];
            for (var i = 0; i < this._長さ; i++)
            {
                l_逆相補[this._長さ - 1 - i] = Util.Get_相補塩基ID(this.Get_塩基ID(i));
            }
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
                l_塩基列[i] = this.Get_塩基ID(i);
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
            if (this._長さ != p_比較対象._長さ)
            {
                return false;
            }
            var l_語数 = (this._長さ + 31) >> 5;
            for (var i = 0; i < l_語数; i++)
            {
                if (this.Get_語(i) != p_比較対象.Get_語(i))
                {
                    return false;
                }
            }
            return true;
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
            var l_ハッシュ = 1_469_598_103_934_665_603UL ^ (ulong)this._長さ;
            var l_語数 = (this._長さ + 31) >> 5;
            for (var i = 0; i < l_語数; i++)
            {
                l_ハッシュ ^= this.Get_語(i);
                l_ハッシュ *= 1_099_511_628_211UL;
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
        private static int Get_比較結果(KmerKey p_左, KmerKey p_右)
        {
            var l_語数 = (p_左._長さ + 31) >> 5;
            for (var i = 0; i < l_語数; i++)
            {
                var l_左語 = p_左.Get_語(i);
                var l_右語 = p_右.Get_語(i);
                if (l_左語 != l_右語)
                {
                    return l_左語 < l_右語 ? -1 : 1;
                }
            }
            return 0;
        }

        /// <summary>
        /// 指定位置のパック済み語を返す
        /// </summary>
        /// <param name="p_位置">語の位置</param>
        /// <returns>パック済み語</returns>
        private ulong Get_語(int p_位置)
        {
            return this._長いパック済みデータ is not null ? this._長いパック済みデータ[p_位置] : p_位置 == 0 ? this._先頭語 : this._第2語;
        }

        /// <summary>
        /// 指定位置の塩基 ID を返す
        /// </summary>
        /// <param name="p_位置">塩基の位置</param>
        /// <returns>塩基 ID</returns>
        private byte Get_塩基ID(int p_位置)
        {
            var l_シフト量 = (31 ^ (p_位置 & 31)) << 1;
            return (byte)(((this.Get_語(p_位置 >> 5) >> l_シフト量) & 0x3UL) + 1UL);
        }

        #endregion
    }
}
