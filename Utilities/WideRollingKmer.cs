using Tsumiki.Commons;
using Tsumiki.Models.Foundation;

namespace Tsumiki.Utilities
{
    /// <summary>
    /// 128 塩基を超える正準キーを 1 塩基ずつ更新する
    /// </summary>
    /// <remarks>
    /// 語の並びは <see cref="KmerKey"/> と同じ (先頭の塩基ほど上位ビット) なので、作業領域をそのままキーとして辞書を引ける<br/>
    /// 窓ごとにキーを作り直すと、詰め直しと逆相補の生成で窓 1 つあたり O (k) の計算と配列の確保が 3 回ずつ要る
    /// </remarks>
    internal sealed class WideRollingKmer
    {
        #region 内部変数

        /// <summary>窓の長さ</summary>
        private readonly int _長さ;
        /// <summary>最後の語で塩基が入っているビット</summary>
        private readonly ulong _末尾語のマスク;
        /// <summary>窓の末尾の塩基を置くシフト量</summary>
        private readonly int _末尾のシフト量;
        /// <summary>順鎖</summary>
        private readonly ulong[] _順;
        /// <summary>逆鎖</summary>
        private readonly ulong[] _逆;
        /// <summary>連続する有効塩基数</summary>
        private int _有効数;

        #endregion

        #region コンストラクタ

        /// <summary>空の窓を作る</summary>
        /// <param name="p_長さ">窓の長さ</param>
        public WideRollingKmer(int p_長さ)
        {
            ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(p_長さ, 128);
            this._長さ = p_長さ;
            var l_語数 = (p_長さ + 31) >> 5;
            this._順 = new ulong[l_語数];
            this._逆 = new ulong[l_語数];
            this._末尾のシフト量 = (31 ^ ((p_長さ - 1) & 31)) << 1;
            this._末尾語のマスク = ulong.MaxValue << this._末尾のシフト量;
        }

        #endregion

        #region 公開メソッド

        /// <summary>末尾へ塩基を足し、完成した窓の正準キーを返す</summary>
        /// <param name="p_塩基">追加する塩基</param>
        /// <param name="p_キー">正準キー、作業領域を指すので次に塩基を足すまでしか使えない</param>
        /// <returns>有効塩基だけで窓が埋まれば true</returns>
        public bool Try追加(char p_塩基, out KmerKey p_キー)
        {
            p_キー = default;
            var l_ID = Util.Get_塩基ID(p_塩基);
            if (l_ID is < 1 or > 4)
            {
                Array.Clear(this._順);
                Array.Clear(this._逆);
                this._有効数 = 0;
                return false;
            }

            var l_値 = (ulong)(l_ID - 1);
            var l_末尾 = this._順.Length - 1;

            // 順鎖は全体を 1 塩基ぶん先頭側へ送り、空いた末尾へ置く
            for (var j = 0; j < l_末尾; j++)
            {
                this._順[j] = (this._順[j] << 2) | (this._順[j + 1] >> 62);
            }
            this._順[l_末尾] = (this._順[l_末尾] << 2) | (l_値 << this._末尾のシフト量);

            // 逆鎖は相補塩基を先頭に置き、全体を末尾側へ送ってはみ出た塩基を落とす
            for (var j = l_末尾; j > 0; j--)
            {
                this._逆[j] = (this._逆[j] >> 2) | (this._逆[j - 1] << 62);
            }
            this._逆[0] = (this._逆[0] >> 2) | ((l_値 ^ 3UL) << 62);
            this._逆[l_末尾] &= this._末尾語のマスク;

            this._有効数 = Math.Min(this._長さ, this._有効数 + 1);
            if (this._有効数 < this._長さ)
            {
                return false;
            }

            p_キー = new KmerKey(this._長さ, Is順が小さい(this._順, this._逆) ? this._順 : this._逆);
            return true;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 語を先頭から比べて、順鎖が逆鎖以下か
        /// </summary>
        /// <param name="p_順"></param>
        /// <param name="p_逆"></param>
        /// <remarks>
        /// 等しいときに順鎖を採るのは <see cref="KmerKey.Get_正規形"/> と揃えるため
        /// </remarks>
        /// <returns></returns>
        private static bool Is順が小さい(ulong[] p_順, ulong[] p_逆)
        {
            for (var i = 0; i < p_順.Length; i++)
            {
                if (p_順[i] != p_逆[i])
                {
                    return p_順[i] < p_逆[i];
                }
            }
            return true;
        }

        #endregion
    }
}
