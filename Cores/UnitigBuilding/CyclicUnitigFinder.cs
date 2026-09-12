using Tsumiki.Commons;
using Tsumiki.Models.UnitigBuilding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.UnitigBuilding
{
    /// <summary>
    /// 分岐を 1 つも持たない閉路を拾い、そこからの走査の開始点を返す
    /// </summary>
    /// <remarks>
    /// unitig の開始点は「入次数が 1 でない、または唯一の予測元が分岐している」k-mer として選ぶ<br/>
    /// 閉路の全頂点が入次数 1 ・出次数 1 で、予測元も分岐していない場合、この条件を満たす k-mer が 1 つも存在せず、閉路が丸ごと走査対象から外れる<br/>
    /// エラーの少ない小さなプラスミドや、きれいな環状染色体がそのまま出力から消える<br/>
    /// 閉路には始点が無いので、どの頂点から始めても同じ環を 1 周する<br/>
    /// 覆われずに残った k-mer を 1 つ選んで開始点にすればよい
    /// </remarks>
    internal static class CyclicUnitigFinder
    {
        #region 公開メソッド

        /// <summary>
        /// p_walk結果 が覆えなかった閉路それぞれについて、開始点を 1 つずつ返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_walk結果"></param>
        /// <param name="p_k長"></param>
        /// <returns>覆い残しが無ければ空</returns>
        public static List<byte[]> Get_閉路の開始kmer(TrustedKmerIndex p_kmerインデックス, IReadOnlyList<string> p_walk結果, int p_k長)
        {
            var l_覆済み = new 正規形集合(p_k長);
            foreach (var l_配列 in p_walk結果)
            {
                V_記録_覆った範囲(l_覆済み, l_配列, p_k長);
            }

            List<byte[]> l_開始kmer = [];
            foreach (var l_kmer in p_kmerインデックス.Get_信頼kmer一覧())
            {
                if (l_覆済み.Get_含まれるか(l_kmer))
                {
                    continue;
                }
                if (V_辿る_閉路(p_kmerインデックス, l_kmer, p_k長, l_覆済み))
                {
                    l_開始kmer.Add(l_kmer);
                }
            }
            return l_開始kmer;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 走査で得た配列が覆った k-mer を記録する
        /// </summary>
        /// <param name="p_覆済み"></param>
        /// <param name="p_配列"></param>
        /// <param name="p_k長"></param>
        /// <remarks>
        /// k &lt;= 64 ではパック値を転がして作る (位置ごとに詰め直すと総延長 x k の手間になる)
        /// </remarks>
        private static void V_記録_覆った範囲(正規形集合 p_覆済み, string p_配列, int p_k長)
        {
            if (p_配列.Length < p_k長)
            {
                return;
            }

            if (p_k長 > 64)
            {
                var l_塩基列 = Util.V_変換_塩基列(p_配列);
                for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
                {
                    p_覆済み.V_追加(l_塩基列.AsSpan(i, p_k長));
                }
                return;
            }

            var l_マスク = p_k長 >= 64 ? UInt128.MaxValue : (UInt128.One << (2 * p_k長)) - 1;
            var l_上位シフト = 2 * (p_k長 - 1);
            UInt128 l_順鎖 = 0;
            UInt128 l_逆鎖 = 0;
            var l_直近の曖昧位置 = -1;

            for (var i = 0; i < p_配列.Length; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_配列[i]);
                var l_有効か = l_塩基ID is >= Consts.塩基ID.A and <= Consts.塩基ID.T;
                var l_コドン = (UInt128)(l_有効か ? l_塩基ID - 1 : 0);
                l_順鎖 = ((l_順鎖 << 2) | l_コドン) & l_マスク;
                l_逆鎖 = (l_逆鎖 >> 2) | ((3 - l_コドン) << l_上位シフト);
                if (!l_有効か)
                {
                    l_直近の曖昧位置 = i;
                }

                var l_開始 = i - p_k長 + 1;
                if (l_開始 >= 0 && l_直近の曖昧位置 < l_開始)
                {
                    p_覆済み.V_追加(l_順鎖 < l_逆鎖 ? l_順鎖 : l_逆鎖);
                }
            }
        }

        /// <summary>
        /// 開始 k-mer から前進し、元へ戻ってくれば閉路として true を返す
        /// </summary>
        /// <remarks>
        /// 通った k-mer は覆済みに入れ、同じ閉路の別の k-mer から二度目の走査が始まらないようにする<br/>
        /// 覆われずに残る k-mer は、どの開始点からも到達されない= 予測元を遡ると必ず輪になる、という性質を満たすものだけなので、途中で出次数が 1 でなくなることは無い<br/>
        /// それでも念のため見るのは、上流の判定が変わったときに無限に回り続けないようにするため
        /// </remarks>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_開始kmer"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_覆済み"></param>
        /// <returns>閉路であれば true</returns>
        private static bool V_辿る_閉路(TrustedKmerIndex p_kmerインデックス, byte[] p_開始kmer, int p_k長, 正規形集合 p_覆済み)
        {
            var l_現在 = (byte[])p_開始kmer.Clone();
            var l_次 = new byte[p_k長];

            while (true)
            {
                p_覆済み.V_追加(l_現在);

                l_現在.AsSpan(1).CopyTo(l_次);
                var l_候補数 = 0;
                byte l_次の塩基 = 0;
                for (var i = Consts.塩基ID.A; i <= Consts.塩基ID.T; i++)
                {
                    l_次[^1] = i;
                    if (p_kmerインデックス.Get_含まれるか(l_次))
                    {
                        l_候補数++;
                        l_次の塩基 = i;
                    }
                }
                if (l_候補数 != 1)
                {
                    return false;
                }

                l_次[^1] = l_次の塩基;
                if (p_覆済み.Get_含まれるか(l_次))
                {
                    // 既に通った所へ戻った
                    // それが出発点なら 1 周できている
                    return 正規形集合.Get_同じ座位か(l_次, p_開始kmer, p_k長);
                }
                l_次.CopyTo(l_現在, 0);
            }
        }

        #endregion
    }
}
