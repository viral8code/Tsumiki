namespace Tsumiki.Utility
{
    /// <summary>
    /// 「見たことがあるか」だけを答える確率的な集合。
    ///
    /// 偽陽性(見ていないものを見たと答える)は起きるが、偽陰性は起きない。
    /// これを使う側は「リードに現れない経路を棄却する拒否権」なので、
    /// 偽陽性は棄却し損ねる方向にしか働かない。長い r-mer では厳密な集合が
    /// 数GBになる一方、この向きの誤りなら安全側に倒れる。
    /// </summary>
    internal sealed class BloomFilter(long p_ビット数, int p_ハッシュ数)
    {
        private readonly ulong[] _ビット = new ulong[(p_ビット数 + 63) / 64];
        private readonly long _ビット数 = p_ビット数;
        private readonly int _ハッシュ数 = p_ハッシュ数;

        public void V_登録(ulong p_値)
        {
            foreach (var l_位置 in this.Get_位置列(p_値))
            {
                this._ビット[l_位置 >> 6] |= 1UL << (int)(l_位置 & 63);
            }
        }

        public bool Get_含まれるか(ulong p_値)
        {
            foreach (var l_位置 in this.Get_位置列(p_値))
            {
                if ((this._ビット[l_位置 >> 6] & (1UL << (int)(l_位置 & 63))) == 0)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// 1つの値から複数のビット位置を作る。2つの独立なハッシュの
        /// 線形結合で済ませる (Kirsch-Mitzenmacher)。
        /// </summary>
        private IEnumerable<long> Get_位置列(ulong p_値)
        {
            var l_h1 = Get_混ぜる(p_値);
            var l_h2 = Get_混ぜる(p_値 ^ 0x9E3779B97F4A7C15UL) | 1UL;
            for (var i = 0; i < this._ハッシュ数; i++)
            {
                yield return (long)((l_h1 + ((ulong)i * l_h2)) % (ulong)this._ビット数);
            }
        }

        private static ulong Get_混ぜる(ulong p_値)
        {
            var l_値 = p_値;
            l_値 ^= l_値 >> 33;
            l_値 *= 0xFF51AFD7ED558CCDUL;
            l_値 ^= l_値 >> 33;
            l_値 *= 0xC4CEB9FE1A85EC53UL;
            l_値 ^= l_値 >> 33;
            return l_値;
        }
    }
}
