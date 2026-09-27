namespace Tsumiki.Utilities
{
    /// <summary>
    /// 集合に無いことを速く答えるブルームフィルタ (1 回の問い合わせで 64 バイトの 1 ブロックだけを読む)
    /// </summary>
    /// <param name="p_件数">入れる見込みの件数</param>
    internal sealed class ブルームフィルタ(long p_件数)
    {
        #region 定数

        /// <summary>
        /// 1 件あたりのビット数
        /// </summary>
        private const long C_1件あたりのビット数 = 12L;

        /// <summary>
        /// 1 件ごとに立てるビットの数
        /// </summary>
        private const int C_立てるビット数 = 6;

        /// <summary>
        /// 1 ブロックの語数 (64 バイト)
        /// </summary>
        private const int C_ブロックの語数 = 8;

        /// <summary>
        /// ブロック数の下限
        /// </summary>
        private const long C_ブロック数の下限 = 1_024L;

        #endregion

        #region 内部変数

        /// <summary>
        /// ブロック数
        /// </summary>
        private readonly ulong _ブロック数 = (ulong)Get_ブロック数(p_件数);

        /// <summary>
        /// ビット列
        /// </summary>
        private readonly ulong[] _語 = new ulong[Get_ブロック数(p_件数) * C_ブロックの語数];

        #endregion

        #region 公開メソッド

        /// <summary>
        /// ハッシュ値を入れる (同時に呼ばない)
        /// </summary>
        /// <param name="p_ハッシュ">入れる値のハッシュ</param>
        public void V_追加(ulong p_ハッシュ)
        {
            var l_先頭 = this.Get_ブロック先頭(p_ハッシュ);
            var l_ビット源 = Get_混合(p_ハッシュ);
            for (var i = 0; i < C_立てるビット数; i++)
            {
                var l_ビット = (int)(l_ビット源 >> (9 * i)) & 511;
                this._語[l_先頭 + (l_ビット >> 6)] |= 1UL << (l_ビット & 63);
            }
        }

        /// <summary>
        /// 入っているかもしれないか (false なら確実に入っていない)
        /// </summary>
        /// <param name="p_ハッシュ">調べる値のハッシュ</param>
        /// <returns>入っているかもしれなければ true</returns>
        public bool Has候補(ulong p_ハッシュ)
        {
            var l_先頭 = this.Get_ブロック先頭(p_ハッシュ);
            var l_ビット源 = Get_混合(p_ハッシュ);
            for (var i = 0; i < C_立てるビット数; i++)
            {
                var l_ビット = (int)(l_ビット源 >> (9 * i)) & 511;
                if ((this._語[l_先頭 + (l_ビット >> 6)] & (1UL << (l_ビット & 63))) == 0)
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// 64 bit の値を混ぜてハッシュにする
        /// </summary>
        /// <param name="p_値">混ぜる値</param>
        /// <returns>ハッシュ</returns>
        public static ulong Get_混合(ulong p_値)
        {
            var l_値 = p_値;
            l_値 ^= l_値 >> 33;
            l_値 *= 0xff51afd7ed558ccdUL;
            l_値 ^= l_値 >> 33;
            l_値 *= 0xc4ceb9fe1a85ec53UL;
            l_値 ^= l_値 >> 33;
            return l_値;
        }

        /// <summary>
        /// 128 bit の値のハッシュ
        /// </summary>
        /// <param name="p_値">値</param>
        /// <returns>ハッシュ</returns>
        public static ulong Get_ハッシュ(UInt128 p_値)
        {
            return Get_混合((ulong)p_値 ^ Get_混合((ulong)(p_値 >> 64)));
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 件数に見合うブロック数
        /// </summary>
        /// <param name="p_件数">入れる見込みの件数</param>
        /// <returns>ブロック数</returns>
        private static long Get_ブロック数(long p_件数)
        {
            return Math.Max(C_ブロック数の下限, ((p_件数 * C_1件あたりのビット数) + 511L) / 512L);
        }

        /// <summary>
        /// ハッシュが入るブロックの先頭の語の位置
        /// </summary>
        /// <param name="p_ハッシュ">ハッシュ</param>
        /// <returns>語の位置</returns>
        private int Get_ブロック先頭(ulong p_ハッシュ)
        {
            return (int)Math.BigMul(p_ハッシュ, this._ブロック数, out _) * C_ブロックの語数;
        }

        #endregion
    }
}
