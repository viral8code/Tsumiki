using Tsumiki.Commons;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 128 塩基までの正準キーを 1 塩基ずつ更新する
    /// </summary>
    internal struct RollingKmer
    {
        #region 内部変数

        /// <summary>窓の長さ</summary>
        private readonly int _長さ;
        /// <summary>上位語の有効ビット</summary>
        private readonly UInt128 _マスク;
        /// <summary>順鎖の上位語</summary>
        private UInt128 _順上;
        /// <summary>順鎖の下位語</summary>
        private UInt128 _順下;
        /// <summary>逆鎖の上位語</summary>
        private UInt128 _逆上;
        /// <summary>逆鎖の下位語</summary>
        private UInt128 _逆下;
        /// <summary>連続する有効塩基数</summary>
        private int _有効数;

        #endregion

        #region コンストラクタ

        /// <summary>空の窓を作る</summary>
        /// <param name="p_長さ">窓の長さ</param>
        public RollingKmer(int p_長さ)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(p_長さ, 1);
            ArgumentOutOfRangeException.ThrowIfGreaterThan(p_長さ, 128);
            this._長さ = p_長さ;
            var l_ビット数 = 2 * (p_長さ <= 64 ? p_長さ : p_長さ - 64);
            this._マスク = l_ビット数 == 128 ? UInt128.MaxValue : ((UInt128)1 << l_ビット数) - 1;
        }

        #endregion

        #region 公開メソッド

        /// <summary>末尾へ塩基を足し、完成した窓の正準キーを返す</summary>
        /// <param name="p_塩基">追加する塩基</param>
        /// <param name="p_キー">正準キー</param>
        /// <returns>有効塩基だけで窓が埋まれば true</returns>
        public bool Try追加(char p_塩基, out (UInt128 A_上位, UInt128 A_下位, string? A_長い配列) p_キー)
        {
            p_キー = default;
            var l_ID = Util.Get_塩基ID(p_塩基);
            if (l_ID is < 1 or > 4)
            {
                this._順上 = this._順下 = this._逆上 = this._逆下 = 0;
                this._有効数 = 0;
                return false;
            }

            var l_値 = (UInt128)(l_ID - 1);
            if (this._長さ <= 64)
            {
                this._順下 = ((this._順下 << 2) | l_値) & this._マスク;
                this._逆下 = (this._逆下 >> 2) | ((l_値 ^ 3) << (2 * (this._長さ - 1)));
            }
            else
            {
                this._順上 = ((this._順上 << 2) | (this._順下 >> 126)) & this._マスク;
                this._順下 = (this._順下 << 2) | l_値;
                this._逆下 = (this._逆下 >> 2) | (this._逆上 << 126);
                this._逆上 = (this._逆上 >> 2) | ((l_値 ^ 3) << (2 * (this._長さ - 65)));
            }

            this._有効数 = Math.Min(this._長さ, this._有効数 + 1);
            if (this._有効数 < this._長さ)
            {
                return false;
            }

            p_キー = this._順上 < this._逆上 || (this._順上 == this._逆上 && this._順下 <= this._逆下)
                ? (this._順上, this._順下, null) : (this._逆上, this._逆下, null);
            return true;
        }

        #endregion
    }
}
