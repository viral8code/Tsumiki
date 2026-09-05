using Tsumiki.Common;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// unitig の walk を、パック値を転がしながら進める実装(k &lt;= 64 用)。
    ///
    /// walk は1塩基ずつ進むので、k-mer のパック値は前の値からシフトで作れる。
    /// Span から毎回詰め直すと、1歩あたり O(k) の詰め直しが所属判定と
    /// 入次数判定の回数だけ走る。
    /// </summary>
    internal sealed class UnitigWalk
    {
        private readonly TrustedKmerIndex _kmerインデックス;

        private readonly int _k長;

        private readonly bool _小経路か;

        private readonly UInt128 _マスク;

        private readonly int _上位シフト;

        public UnitigWalk(TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            this._kmerインデックス = p_kmerインデックス;
            this._k長 = p_k長;
            this._小経路か = p_k長 <= 32;
            this._マスク = p_k長 >= 64 ? UInt128.MaxValue : ((UInt128)1 << (2 * p_k長)) - 1;
            this._上位シフト = 2 * (p_k長 - 1);
        }

        /// <summary>この実装で扱える k かどうか。</summary>
        public static bool Get_扱えるか(int p_k長) => p_k長 <= 64;

        private bool Get_含まれるか(UInt128 p_順鎖, UInt128 p_逆鎖)
        {
            var l_正規形 = p_順鎖 < p_逆鎖 ? p_順鎖 : p_逆鎖;
            return this._小経路か
                ? this._kmerインデックス.Get_含まれるか_小((ulong)l_正規形)
                : this._kmerインデックス.Get_含まれるか_中(l_正規形);
        }

        /// <summary>末尾に塩基を足した k-mer の順鎖・逆鎖パック値。</summary>
        private (UInt128 A_順鎖, UInt128 A_逆鎖) Get_後続(UInt128 p_順鎖, UInt128 p_逆鎖, byte p_塩基ID)
        {
            var l_コドン = (UInt128)(p_塩基ID - 1);
            var l_順鎖 = ((p_順鎖 << 2) | l_コドン) & this._マスク;
            var l_逆鎖 = (p_逆鎖 >> 2) | ((3 - l_コドン) << this._上位シフト);
            return (l_順鎖, l_逆鎖);
        }

        /// <summary>先頭に塩基を足した(末尾を落とした)k-mer の順鎖・逆鎖パック値。</summary>
        private (UInt128 A_順鎖, UInt128 A_逆鎖) Get_予測元(UInt128 p_順鎖, UInt128 p_逆鎖, byte p_塩基ID)
        {
            var l_コドン = (UInt128)(p_塩基ID - 1);
            var l_順鎖 = (p_順鎖 >> 2) | (l_コドン << this._上位シフト);
            var l_逆鎖 = ((p_逆鎖 << 2) | (3 - l_コドン)) & this._マスク;
            return (l_順鎖, l_逆鎖);
        }

        /// <summary>
        /// 入次数がちょうど1かどうか。
        /// 前進規則が 後続 = kmer[1..] + c である以上、その逆を解くと
        /// 予測元は c + kmer[..^1] になる。
        /// </summary>
        private bool Get_入次数が1か(UInt128 p_順鎖, UInt128 p_逆鎖)
        {
            var l_件数 = 0;
            for (byte i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
            {
                var (l_元順, l_元逆) = this.Get_予測元(p_順鎖, p_逆鎖, i);
                if (this.Get_含まれるか(l_元順, l_元逆) && ++l_件数 > 1)
                {
                    return false;
                }
            }
            return l_件数 == 1;
        }

        /// <summary>
        /// 開始 k-mer から前進 walk して unitig の塩基列を返す。
        /// 循環を検出したら打ち切る。
        /// </summary>
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
                    // 循環。1塩基前で打ち切った場合と同じ配列になるよう末尾を外す。
                    l_配列.RemoveAt(l_配列.Count - 1);
                    return l_配列;
                }

                byte l_次の塩基 = 0;
                var l_候補数 = 0;
                UInt128 l_次順 = 0;
                UInt128 l_次逆 = 0;
                for (byte i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
                {
                    var (l_順, l_逆) = this.Get_後続(l_順鎖, l_逆鎖, i);
                    if (!this.Get_含まれるか(l_順, l_逆))
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

                // 出次数が1でなければ、ここが unitig の終端。
                // 次の k-mer の入次数が2以上なら別の経路が合流しており、
                // そこからは別の unitig が始まるのでやはり終端になる。
                if (l_候補数 != 1 || !this.Get_入次数が1か(l_次順, l_次逆))
                {
                    return l_配列;
                }

                l_配列.Add(l_次の塩基);
                (l_順鎖, l_逆鎖) = (l_次順, l_次逆);
            }
        }
    }
}
