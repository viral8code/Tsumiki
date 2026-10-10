using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 配列を 1 塩基ずつ送り、向きをそのままにした k-mer のキーを詰め直さずに作る (A・C・G・T 以外で途切れる)
    /// </summary>
    /// <param name="p_長さ">k-mer の長さ</param>
    internal sealed class 順鎖Kmer転がし(int p_長さ)
    {
        #region 内部変数

        /// <summary>
        /// 末尾の塩基を置く位置のシフト量
        /// </summary>
        private readonly int _末尾のシフト量 = (31 ^ ((p_長さ - 1) & 31)) << 1;

        /// <summary>
        /// 塩基を先頭から 2 bit ずつ上位側へ詰めた語
        /// </summary>
        private readonly ulong[] _語 = new ulong[(p_長さ + 31) >> 5];

        /// <summary>
        /// 途切れてから続けて送った塩基の数
        /// </summary>
        private int _有効数;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 1 塩基送る
        /// </summary>
        /// <param name="p_塩基">送る塩基の文字</param>
        /// <param name="p_キー">k 塩基そろったときのキー (k が 64 を超えると作業領域を参照するので、次に送るまでに使い終える)</param>
        /// <returns>k 塩基そろったら true</returns>
        public bool Is成功_追加(char p_塩基, out KmerKey p_キー)
        {
            p_キー = default;
            var l_ID = Util.Get_塩基ID(p_塩基);
            if (l_ID is < Consts.塩基ID.A or > Consts.塩基ID.T)
            {
                Array.Clear(this._語);
                this._有効数 = 0;
                return false;
            }

            var l_末尾 = this._語.Length - 1;
            for (var j = 0; j < l_末尾; j++)
            {
                this._語[j] = (this._語[j] << 2) | (this._語[j + 1] >> 62);
            }

            this._語[l_末尾] = (this._語[l_末尾] << 2) | ((ulong)(l_ID - 1) << this._末尾のシフト量);
            this._有効数 = Math.Min(p_長さ, this._有効数 + 1);
            if (this._有効数 < p_長さ)
            {
                return false;
            }

            p_キー = new KmerKey(p_長さ, this._語);
            return true;
        }

        #endregion
    }
}
