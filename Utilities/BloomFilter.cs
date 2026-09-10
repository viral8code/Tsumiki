namespace Tsumiki.Utilities
{
    /// <summary>
    /// 「見たことがあるか」だけを答える確率的な集合
    /// </summary>
    /// <remarks>
    /// 偽陽性 (見ていないものを見たと答える) は起きるが、偽陰性は起きない<br/>
    /// これを使う側はリードに現れない経路を棄却する拒否権なので、
    /// 偽陽性は棄却し損ねる方向にしか働かない<br/>
    /// 長い r-mer では厳密な集合が数 GB になる一方、この向きの誤りなら安全側に倒れる
    /// </remarks>
    /// <param name="p_ビット数">確保するビット数</param>
    /// <param name="p_ハッシュ数">1 つの値に対して立てるビットの数</param>
    internal sealed class BloomFilter(long p_ビット数, int p_ハッシュ数)
    {
        #region 内部変数

        /// <summary>
        /// ビット列の実体
        /// </summary>
        private readonly ulong[] _ビット = new ulong[((p_ビット数 + 63L) / 64L)];

        /// <summary>
        /// 確保したビット数
        /// </summary>
        private readonly long _ビット数 = p_ビット数;

        /// <summary>
        /// 1 つの値に対して立てるビットの数
        /// </summary>
        private readonly int _ハッシュ数 = p_ハッシュ数;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 値を登録する
        /// </summary>
        /// <param name="p_値">登録する値</param>
        public void V_登録(ulong p_値)
        {
            foreach (var l_位置 in this.Get_位置列(p_値))
            {
                this._ビット[l_位置 >> 6] |= 1UL << (int)(l_位置 & 63L);
            }
        }

        /// <summary>
        /// 値が登録されている可能性があるか
        /// </summary>
        /// <param name="p_値">調べる値</param>
        /// <returns>登録されていれば必ず true、未登録でも偽陽性で true になりうる</returns>
        public bool Get_含まれるか(ulong p_値)
        {
            foreach (var l_位置 in this.Get_位置列(p_値))
            {
                if ((this._ビット[l_位置 >> 6] & (1UL << (int)(l_位置 & 63L))) == 0UL)
                {
                    return false;
                }
            }
            return true;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 つの値から複数のビット位置を作る
        /// </summary>
        /// <remarks>
        /// 2 つの独立なハッシュの線形結合で済ませる (Kirsch-Mitzenmacher)
        /// </remarks>
        /// <param name="p_値">元の値</param>
        /// <returns>立てる、あるいは調べるビットの位置</returns>
        private IEnumerable<long> Get_位置列(ulong p_値)
        {
            var l_ハッシュ1 = Get_ハッシュ値(p_値);
            var l_ハッシュ2 = Get_ハッシュ値(p_値 ^ 0x9E3779B97F4A7C15UL) | 1UL;
            for (var i = 0; i < this._ハッシュ数; i++)
            {
                yield return (long)((l_ハッシュ1 + ((ulong)i * l_ハッシュ2)) % (ulong)this._ビット数);
            }
        }

        /// <summary>
        /// ビットが偏らないように値をかき混ぜる
        /// </summary>
        /// <param name="p_値">元の値</param>
        /// <returns>かき混ぜた値</returns>
        private static ulong Get_ハッシュ値(ulong p_値)
        {
            var l_値 = p_値;
            l_値 ^= l_値 >> 33;
            l_値 *= 0xFF51AFD7ED558CCDUL;
            l_値 ^= l_値 >> 33;
            l_値 *= 0xC4CEB9FE1A85EC53UL;
            l_値 ^= l_値 >> 33;
            return l_値;
        }

        #endregion
    }
}
