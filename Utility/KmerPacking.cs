using Tsumiki.Common;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer を 2bit へパックし、逆相補と比べて小さいほう(正規形)へ寄せる。
    /// 集計のキーに文字列を使うとアセンブリ規模で 1GB を超えるため、
    /// k-mer を数える処理は常にこのパック済みの値をキーにする。
    /// パックは k &lt;= 64 でしか使えないので、それを超える長さには
    /// Get_正規化キー(正規形の 128bit ハッシュ)を使う。
    /// </summary>
    internal static class KmerPacking
    {
        /// <summary>
        /// 配列の位置 p_開始位置 から p_k長 塩基を 2bit パックし、正規形を返す。
        /// 曖昧塩基(N など)を含む場合は false を返す。
        /// </summary>
        public static bool Get_正規化パック(string p_配列, int p_開始位置, int p_k長, out UInt128 p_正規形)
        {
            UInt128 l_順鎖 = 0;
            for (var i = 0; i < p_k長; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_配列[p_開始位置 + i]);
                if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    p_正規形 = 0;
                    return false;
                }
                l_順鎖 = (l_順鎖 << 2) | (UInt128)(l_塩基ID - 1);
            }
            p_正規形 = Get_小さいほう(l_順鎖, p_k長);
            return true;
        }

        /// <summary>
        /// 塩基ID列(1..4)を 2bit パックし、正規形を返す。
        /// </summary>
        public static UInt128 Get_正規化パック(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_順鎖 = 0;
            foreach (var l_塩基ID in p_kmer)
            {
                l_順鎖 = (l_順鎖 << 2) | (UInt128)(l_塩基ID - 1);
            }
            return Get_小さいほう(l_順鎖, p_kmer.Length);
        }

        /// <summary>
        /// 長さに縛られない集計キー。k &lt;= 64 ならパック済みの正規形そのもの、
        /// それを超えるなら正規形の 128bit ハッシュを返す。
        ///
        /// ハッシュに落とすのは、数える・存在を問うだけで配列を戻さない用途に
        /// 限る。128bit なら数千万種類でも衝突は事実上起きない。
        /// </summary>
        public static UInt128 Get_正規化キー(ReadOnlySpan<byte> p_kmer)
        {
            return p_kmer.Length <= 64 ? Get_正規化パック(p_kmer) : Get_正規化ハッシュ(p_kmer);
        }

        /// <summary>
        /// 配列の位置 p_開始位置 から p_k長 塩基の集計キー。
        /// 曖昧塩基を含む場合は false を返す。
        /// </summary>
        public static bool Get_正規化キー(string p_配列, int p_開始位置, int p_k長, out UInt128 p_キー)
        {
            if (p_k長 <= 64)
            {
                return Get_正規化パック(p_配列, p_開始位置, p_k長, out p_キー);
            }

            var l_塩基列 = new byte[p_k長];
            for (var i = 0; i < p_k長; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_配列[p_開始位置 + i]);
                if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    p_キー = 0;
                    return false;
                }
                l_塩基列[i] = l_塩基ID;
            }
            p_キー = Get_正規化ハッシュ(l_塩基列);
            return true;
        }

        /// <summary>
        /// 正規形(順鎖と逆相補のうち辞書順で小さいほう)を 64bit へ畳む。
        /// 32 塩基以下なら 2bit パックそのもので、衝突は起きない。
        /// </summary>
        public static ulong Get_正規化ハッシュ_64(ReadOnlySpan<byte> p_kmer)
        {
            if (p_kmer.Length <= 32)
            {
                return (ulong)Get_正規化パック(p_kmer);
            }
            return (ulong)Get_正規化ハッシュ(p_kmer);
        }

        /// <summary>
        /// 正規形のバイト列を 128bit へ畳む。逆相補は実際には作らず、
        /// どちら向きが小さいかを決めてからその向きで畳む。
        /// </summary>
        private static UInt128 Get_正規化ハッシュ(ReadOnlySpan<byte> p_kmer)
        {
            var l_順鎖が小さいか = Get_順鎖が小さいか(p_kmer);
            var l_上位 = 14695981039346656037UL;
            var l_下位 = 1099511628211UL;
            for (var i = 0; i < p_kmer.Length; i++)
            {
                var l_塩基 = l_順鎖が小さいか
                    ? p_kmer[i]
                    : (byte)(5 - p_kmer[p_kmer.Length - 1 - i]);
                l_上位 = (l_上位 ^ l_塩基) * 1099511628211UL;
                l_下位 = (l_下位 ^ l_塩基) * 14695981039346656037UL;
            }
            return ((UInt128)l_上位 << 64) | l_下位;
        }

        /// <summary>
        /// 順鎖と逆相補を辞書順で比べる。塩基IDは A=1..T=4 で相補は 5 - ID。
        /// </summary>
        private static bool Get_順鎖が小さいか(ReadOnlySpan<byte> p_kmer)
        {
            for (var i = 0; i < p_kmer.Length; i++)
            {
                var l_逆 = (byte)(5 - p_kmer[p_kmer.Length - 1 - i]);
                if (p_kmer[i] != l_逆)
                {
                    return p_kmer[i] < l_逆;
                }
            }
            return true;
        }

        /// <summary>
        /// 配列の位置 p_開始位置 から p_k長 塩基を、正規化せず順鎖のまま
        /// 2bit パックする。曖昧塩基を含む場合は false を返す。
        ///
        /// 向きを区別したい索引(どちらの鎖に載ったのかで座標の解釈が変わる場合)
        /// では正規形を使えないため、順鎖と逆相補を別々のキーとして扱う。
        /// </summary>
        public static bool Get_パック(string p_配列, int p_開始位置, int p_k長, out UInt128 p_順鎖)
        {
            UInt128 l_順鎖 = 0;
            for (var i = 0; i < p_k長; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_配列[p_開始位置 + i]);
                if (l_塩基ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
                {
                    p_順鎖 = 0;
                    return false;
                }
                l_順鎖 = (l_順鎖 << 2) | (UInt128)(l_塩基ID - 1);
            }
            p_順鎖 = l_順鎖;
            return true;
        }

        /// <summary>
        /// パック済みの値の逆相補。
        /// </summary>
        public static UInt128 Get_逆相補(UInt128 p_パック済み, int p_長さ)
        {
            var l_残り = p_パック済み;
            UInt128 l_逆相補 = 0;
            for (var i = 0; i < p_長さ; i++)
            {
                var l_コドン = l_残り & 3;
                l_逆相補 = (l_逆相補 << 2) | (l_コドン ^ 3);
                l_残り >>= 2;
            }
            return l_逆相補;
        }

        /// <summary>
        /// パック済みの値とその逆相補のうち小さいほうを返す。
        /// </summary>
        public static UInt128 Get_小さいほう(UInt128 p_パック済み, int p_長さ)
        {
            var l_逆相補 = Get_逆相補(p_パック済み, p_長さ);
            return p_パック済み < l_逆相補 ? p_パック済み : l_逆相補;
        }
    }
}
