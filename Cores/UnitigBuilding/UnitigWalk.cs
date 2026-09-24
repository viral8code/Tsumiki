using Tsumiki.Commons;
using Tsumiki.Models.UnitigBuilding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// unitig の walk を、パック値を転がしながら進める実装 (k &lt;= 128 用)
    /// </summary>
    /// <param name="p_kmerインデックス"></param>
    /// <param name="p_k長"></param>
    internal sealed class UnitigWalk(TrustedKmerIndex p_kmerインデックス, int p_k長)
    {
        #region 内部変数

        /// <summary>
        /// 小経路か
        /// </summary>
        private readonly bool _Is小経路 = p_k長 <= 32;

        /// <summary>
        /// マスク
        /// </summary>
        private readonly UInt128 _マスク = p_k長 >= 64 ? UInt128.MaxValue : ((UInt128)1 << (2 * p_k長)) - 1;

        /// <summary>
        /// 128 bit を超える側のマスク (k &gt; 64)
        /// </summary>
        private readonly UInt128 _上位語マスク = p_k長 <= 64 ? 0 : p_k長 >= 128 ? UInt128.MaxValue : ((UInt128)1 << ((2 * p_k長) - 128)) - 1;

        /// <summary>
        /// 上位シフト
        /// </summary>
        private readonly int _上位シフト = (p_k長 - 1) << 1;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// この実装で扱える k かどうか
        /// </summary>
        /// <param name="p_k長">k 長</param>
        /// <returns>扱えれば true</returns>
        public static bool Is対応k長(int p_k長)
        {
            return p_k長 <= TrustedKmerIndex.パック値のk上限;
        }

        /// <summary>
        /// 開始 k-mer から前進 walk して unitig の塩基列を返す (k &lt;= 64)
        /// </summary>
        /// <param name="p_開始kmer">walk の起点となる k-mer</param>
        /// <param name="p_訪問済み">循環検出用の作業集合</param>
        /// <returns>組み立てた塩基列</returns>
        public List<byte> Get_塩基列(ReadOnlySpan<byte> p_開始kmer, HashSet<UInt128> p_訪問済み)
        {
            p_訪問済み.Clear();
            List<byte> l_配列 = new(p_開始kmer.Length * 2);
            UInt128 l_順鎖 = 0;
            UInt128 l_逆鎖 = 0;
            foreach (var l_塩基ID in p_開始kmer)
            {
                l_配列.Add(l_塩基ID);
                (l_順鎖, l_逆鎖) = this.Get_後続(l_順鎖, l_逆鎖, l_塩基ID);
            }

            while (true)
            {
                if (!p_訪問済み.Add(l_順鎖))
                {
                    l_配列.RemoveAt(l_配列.Count - 1);
                    return l_配列;
                }

                byte l_次の塩基 = 0;
                var l_候補数 = 0;
                UInt128 l_次順 = 0;
                UInt128 l_次逆 = 0;
                for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
                {
                    var (l_順, l_逆) = this.Get_後続(l_順鎖, l_逆鎖, i);
                    if (!this.Haskmer(l_順, l_逆))
                    {
                        continue;
                    }

                    if (++l_候補数 > 1)
                    {
                        break;
                    }
                    l_次の塩基 = i;
                    (l_次順, l_次逆) = (l_順, l_逆);
                }

                if (l_候補数 != 1 || !this.Is入次数1(l_次順, l_次逆))
                {
                    return l_配列;
                }

                l_配列.Add(l_次の塩基);
                (l_順鎖, l_逆鎖) = (l_次順, l_次逆);
            }
        }

        /// <summary>
        /// 開始 k-mer から前進 walk して unitig の塩基列を返す (65 &lt;= k &lt;= 128)
        /// </summary>
        /// <param name="p_開始kmer">walk の起点となる k-mer</param>
        /// <param name="p_訪問済み">循環検出用の作業集合</param>
        /// <returns>組み立てた塩基列</returns>
        public List<byte> Get_塩基列_長(ReadOnlySpan<byte> p_開始kmer, HashSet<(UInt128 A_上位, UInt128 A_下位)> p_訪問済み)
        {
            p_訪問済み.Clear();
            List<byte> l_配列 = new(p_開始kmer.Length * 2);
            var l_状態 = default(長状態);
            foreach (var l_塩基ID in p_開始kmer)
            {
                l_配列.Add(l_塩基ID);
                l_状態 = this.Get_後続_長(l_状態, l_塩基ID);
            }

            while (true)
            {
                if (!p_訪問済み.Add((l_状態.A_順上, l_状態.A_順下)))
                {
                    l_配列.RemoveAt(l_配列.Count - 1);
                    return l_配列;
                }

                byte l_次の塩基 = 0;
                var l_候補数 = 0;
                var l_次 = default(長状態);
                for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
                {
                    var l_後続 = this.Get_後続_長(l_状態, i);
                    if (!this.Haskmer_長(l_後続))
                    {
                        continue;
                    }

                    if (++l_候補数 > 1)
                    {
                        break;
                    }
                    l_次の塩基 = i;
                    l_次 = l_後続;
                }

                if (l_候補数 != 1 || !this.Is入次数1_長(l_次))
                {
                    return l_配列;
                }

                l_配列.Add(l_次の塩基);
                l_状態 = l_次;
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 既に通った k-mer か
        /// </summary>
        /// <param name="p_順鎖">順鎖のパック値</param>
        /// <param name="p_逆鎖">逆鎖のパック値</param>
        /// <returns>通っていれば true</returns>
        private bool Haskmer(UInt128 p_順鎖, UInt128 p_逆鎖)
        {
            var l_正規形 = p_順鎖 < p_逆鎖 ? p_順鎖 : p_逆鎖;
            return this._Is小経路
                ? p_kmerインデックス.Haskmer_小((ulong)l_正規形)
                : p_kmerインデックス.Haskmer_中(l_正規形);
        }

        /// <summary>
        /// 末尾に塩基を足した k-mer の順鎖・逆鎖パック値
        /// </summary>
        /// <param name="p_順鎖"></param>
        /// <param name="p_逆鎖"></param>
        /// <param name="p_塩基ID"></param>
        /// <returns></returns>
        private (UInt128 A_順鎖, UInt128 A_逆鎖) Get_後続(UInt128 p_順鎖, UInt128 p_逆鎖, byte p_塩基ID)
        {
            var l_コドン = (UInt128)(p_塩基ID - 1);
            var l_順鎖 = ((p_順鎖 << 2) | l_コドン) & this._マスク;
            var l_逆鎖 = (p_逆鎖 >> 2) | ((3 - l_コドン) << this._上位シフト);
            return (l_順鎖, l_逆鎖);
        }

        /// <summary>
        /// 先頭に塩基を足した (末尾を落とした) k-mer の順鎖・逆鎖パック値
        /// </summary>
        /// <param name="p_順鎖"></param>
        /// <param name="p_逆鎖"></param>
        /// <param name="p_塩基ID"></param>
        /// <returns></returns>
        private (UInt128 A_順鎖, UInt128 A_逆鎖) Get_予測元(UInt128 p_順鎖, UInt128 p_逆鎖, byte p_塩基ID)
        {
            var l_コドン = (UInt128)(p_塩基ID - 1);
            var l_順鎖 = (p_順鎖 >> 2) | (l_コドン << this._上位シフト);
            var l_逆鎖 = ((p_逆鎖 << 2) | (3 - l_コドン)) & this._マスク;
            return (l_順鎖, l_逆鎖);
        }

        /// <summary>
        /// 入次数がちょうど 1 かどうか
        /// </summary>
        /// <param name="p_順鎖"></param>
        /// <param name="p_逆鎖"></param>
        /// <returns></returns>
        private bool Is入次数1(UInt128 p_順鎖, UInt128 p_逆鎖)
        {
            var l_件数 = 0;
            for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                var (l_元順, l_元逆) = this.Get_予測元(p_順鎖, p_逆鎖, i);
                if (this.Haskmer(l_元順, l_元逆) && ++l_件数 > 1)
                {
                    return false;
                }
            }
            return l_件数 == 1;
        }

        /// <summary>
        /// 集合に含まれる k-mer か (k &gt; 64)
        /// </summary>
        /// <param name="p_状態"></param>
        /// <returns></returns>
        private bool Haskmer_長(長状態 p_状態)
        {
            var l_Is順鎖 = p_状態.A_順上 < p_状態.A_逆上 || (p_状態.A_順上 == p_状態.A_逆上 && p_状態.A_順下 <= p_状態.A_逆下);
            return l_Is順鎖
                ? p_kmerインデックス.Haskmer_正規形(p_状態.A_順上, p_状態.A_順下)
                : p_kmerインデックス.Haskmer_正規形(p_状態.A_逆上, p_状態.A_逆下);
        }

        /// <summary>
        /// 末尾に塩基を足した k-mer の順鎖・逆鎖パック値 (k &gt; 64)
        /// </summary>
        /// <param name="p_状態"></param>
        /// <param name="p_塩基ID"></param>
        /// <returns></returns>
        private 長状態 Get_後続_長(長状態 p_状態, byte p_塩基ID)
        {
            var l_コドン = (UInt128)(p_塩基ID - 1);
            var l_順上 = ((p_状態.A_順上 << 2) | (p_状態.A_順下 >> 126)) & this._上位語マスク;
            var l_順下 = (p_状態.A_順下 << 2) | l_コドン;
            var l_逆下 = (p_状態.A_逆下 >> 2) | (p_状態.A_逆上 << 126);
            var l_逆上 = (p_状態.A_逆上 >> 2) | ((3 - l_コドン) << (this._上位シフト - 128));
            return new 長状態(l_順上, l_順下, l_逆上, l_逆下);
        }

        /// <summary>
        /// 先頭に塩基を足した (末尾を落とした) k-mer の順鎖・逆鎖パック値 (k &gt; 64)
        /// </summary>
        /// <param name="p_状態"></param>
        /// <param name="p_塩基ID"></param>
        /// <returns></returns>
        private 長状態 Get_予測元_長(長状態 p_状態, byte p_塩基ID)
        {
            var l_コドン = (UInt128)(p_塩基ID - 1);
            var l_順下 = (p_状態.A_順下 >> 2) | (p_状態.A_順上 << 126);
            var l_順上 = (p_状態.A_順上 >> 2) | (l_コドン << (this._上位シフト - 128));
            var l_逆上 = ((p_状態.A_逆上 << 2) | (p_状態.A_逆下 >> 126)) & this._上位語マスク;
            var l_逆下 = (p_状態.A_逆下 << 2) | (3 - l_コドン);
            return new 長状態(l_順上, l_順下, l_逆上, l_逆下);
        }

        /// <summary>
        /// 入次数がちょうど 1 かどうか (k &gt; 64)
        /// </summary>
        /// <param name="p_状態"></param>
        /// <returns></returns>
        private bool Is入次数1_長(長状態 p_状態)
        {
            var l_件数 = 0;
            for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                if (this.Haskmer_長(this.Get_予測元_長(p_状態, i)) && ++l_件数 > 1)
                {
                    return false;
                }
            }
            return l_件数 == 1;
        }

        #endregion
    }
}
