using Tsumiki.Common;

namespace Tsumiki.Utility
{
    /// <summary>
    /// k-mer を 2bit へパックし、逆相補と比べて小さいほう(正規形)へ寄せる。
    /// 集計のキーに文字列を使うとアセンブリ規模で 1GB を超えるため、
    /// k-mer を数える処理は常にこのパック済みの値をキーにする。
    /// k &lt;= 64 でしか使えない。
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

        /// <summary>塩基ID列(1..4)を 2bit パックし、正規形を返す。</summary>
        public static UInt128 Get_正規化パック(ReadOnlySpan<byte> p_kmer)
        {
            UInt128 l_順鎖 = 0;
            foreach (var l_塩基ID in p_kmer)
            {
                l_順鎖 = (l_順鎖 << 2) | (UInt128)(l_塩基ID - 1);
            }
            return Get_小さいほう(l_順鎖, p_kmer.Length);
        }

        /// <summary>パック済みの値とその逆相補のうち小さいほうを返す。</summary>
        public static UInt128 Get_小さいほう(UInt128 p_パック済み, int p_長さ)
        {
            var l_残り = p_パック済み;
            UInt128 l_逆相補 = 0;
            for (var i = 0; i < p_長さ; i++)
            {
                var l_コドン = l_残り & 3;
                l_逆相補 = (l_逆相補 << 2) | (l_コドン ^ 3);
                l_残り >>= 2;
            }
            return p_パック済み < l_逆相補 ? p_パック済み : l_逆相補;
        }
    }
}
