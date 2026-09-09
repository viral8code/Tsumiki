using Tsumiki.Common;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer を 2 bit へパックし、逆相補と比べて小さいほう (正規形) へ寄せる
    /// </summary>
    /// <remarks>
    /// 集計のキーに文字列を使うとアセンブリ規模で 1 GB を超えるため、
    /// k-mer を数える処理は常にこのパック済みの値をキーにする<br/>
    /// パックは k &lt;= 64 でしか使えないので、
    /// それを超える長さには <see cref="Get_正規化キー(ReadOnlySpan{byte})"/> を使う
    /// </remarks>
    internal static class KmerPacking
    {
        /// <summary>
        /// 配列の位置から k 塩基を 2 bit パックし、正規形を返す
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_開始位置">パックを始める位置</param>
        /// <param name="p_k長">パックする長さ</param>
        /// <param name="p_正規形">パックした正規形</param>
        /// <returns>曖昧塩基 (N など) を含む場合は false</returns>
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
        /// 塩基 ID 列を 2 bit パックし、正規形を返す
        /// </summary>
        /// <param name="p_kmer">塩基 ID 列 (A = 1 .. T = 4)</param>
        /// <returns>パックした正規形</returns>
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
        /// 長さに縛られない集計キー
        /// </summary>
        /// <remarks>
        /// ハッシュに落とすのは、数える、存在を問うだけで配列を戻さない用途に限る<br/>
        /// 128 bit なら数千万種類でも衝突は事実上起きない
        /// </remarks>
        /// <param name="p_kmer">塩基 ID 列 (A = 1 .. T = 4)</param>
        /// <returns>k &lt;= 64 ならパック済みの正規形そのもの、それを超えるなら正規形の 128 bit ハッシュ</returns>
        public static UInt128 Get_正規化キー(ReadOnlySpan<byte> p_kmer)
        {
            return p_kmer.Length <= 64 ? Get_正規化パック(p_kmer) : Get_正規化ハッシュ(p_kmer);
        }

        /// <summary>
        /// 配列の位置から k 塩基ぶんの集計キーを作る
        /// </summary>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_開始位置">キーを作り始める位置</param>
        /// <param name="p_k長">キーにする長さ</param>
        /// <param name="p_キー">作った集計キー</param>
        /// <returns>曖昧塩基を含む場合は false</returns>
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
        /// 正規形 (順鎖と逆相補のうち辞書順で小さいほう) を 64 bit へ畳む
        /// </summary>
        /// <remarks>
        /// 32 塩基以下なら 2 bit パックそのもので、衝突は起きない
        /// </remarks>
        /// <param name="p_kmer">塩基 ID 列 (A = 1 .. T = 4)</param>
        /// <returns>畳んだ値</returns>
        public static ulong Get_正規化ハッシュ_64(ReadOnlySpan<byte> p_kmer)
        {
            return p_kmer.Length <= 32 ? (ulong)Get_正規化パック(p_kmer) : (ulong)Get_正規化ハッシュ(p_kmer);
        }

        /// <summary>
        /// 正規形のバイト列を 128 bit へ畳む
        /// </summary>
        /// <remarks>
        /// 逆相補は実際には作らず、どちら向きが小さいかを決めてからその向きで畳む
        /// </remarks>
        /// <param name="p_kmer">塩基 ID 列 (A = 1 .. T = 4)</param>
        /// <returns>畳んだ値</returns>
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
        /// 順鎖と逆相補を辞書順で比べる
        /// </summary>
        /// <remarks>
        /// 塩基 ID は A = 1 .. T = 4 で、相補は 5 - ID になる
        /// </remarks>
        /// <param name="p_kmer">塩基 ID 列 (A = 1 .. T = 4)</param>
        /// <returns>順鎖のほうが小さいか等しければ true</returns>
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
        /// 配列の位置から k 塩基を、正規化せず順鎖のまま 2 bit パックする
        /// </summary>
        /// <remarks>
        /// 向きを区別したい索引 (どちらの鎖に載ったのかで座標の解釈が変わる場合) では正規形を使えないため、
        /// 順鎖と逆相補を別々のキーとして扱う
        /// </remarks>
        /// <param name="p_配列">元の配列</param>
        /// <param name="p_開始位置">パックを始める位置</param>
        /// <param name="p_k長">パックする長さ</param>
        /// <param name="p_順鎖">パックした順鎖の値</param>
        /// <returns>曖昧塩基を含む場合は false</returns>
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
        /// パック済みの値の逆相補を返す
        /// </summary>
        /// <param name="p_パック済み">パック済みの値</param>
        /// <param name="p_長さ">パックした塩基の数</param>
        /// <returns>逆相補をパックした値</returns>
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
        /// パック済みの値とその逆相補のうち小さいほうを返す
        /// </summary>
        /// <param name="p_パック済み">パック済みの値</param>
        /// <param name="p_長さ">パックした塩基の数</param>
        /// <returns>正規形</returns>
        public static UInt128 Get_小さいほう(UInt128 p_パック済み, int p_長さ)
        {
            var l_逆相補 = Get_逆相補(p_パック済み, p_長さ);
            return p_パック済み < l_逆相補 ? p_パック済み : l_逆相補;
        }
    }
}
