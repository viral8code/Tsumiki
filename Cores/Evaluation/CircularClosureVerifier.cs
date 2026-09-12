using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 環状に閉じたと判定された配列について、その閉じ目を元リードが実際に読んでいるかを確かめる
    /// </summary>
    internal static class CircularClosureVerifier
    {
        #region 定数

        /// <summary>
        /// 閉じ目の左右それぞれに要求する踏み込みの長さ
        /// </summary>
        /// <remarks>
        /// 窓はこの 2 倍になる (2 bit パックが UInt128 に収まる範囲に収める)
        /// </remarks>
        private const int 接合フランク長 = 30;

        /// <summary>
        /// 閉じ目を跨いだとみなす窓の長さ
        /// </summary>
        private const int 接合窓長 = 接合フランク長 << 1;

        /// <summary>
        /// 閉じ目を支持されたとみなすのに必要なリード本数
        /// </summary>
        /// <remarks>
        /// 通常の接合点より高く取る<br/>
        /// ここが誤っていると、アセンブリ全体の見え方が「完全長」から「1 本の線状断片」へ変わってしまうため
        /// </remarks>
        private const int 閉じ目に必要なリード数 = 5;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_FASTAパス の環状配列それぞれについて、閉じ目を跨ぐリードを数える
        /// </summary>
        /// <param name="p_FASTAパス">検証する配列を含む FASTA のパス</param>
        /// <param name="p_リード1のパス">支持を数えるリードのパス</param>
        /// <param name="p_リード2のパス">ペアの相方のパス、無ければ null</param>
        /// <returns>配列ごとの環状閉鎖検証結果</returns>
        /// <remarks>
        /// 環状の配列が 1 本も無ければ空を返す
        /// </remarks>
        public static IReadOnlyList<環状閉鎖検証結果> Get_検証結果(string p_FASTAパス, string p_リード1のパス, string? p_リード2のパス)
        {
            var l_エントリ群 = FastaReader.Get_全エントリ(p_FASTAパス);

            List<(string A_ID, int A_長さ)> l_対象 = [];
            Dictionary<UInt128, int> l_接合窓 = [];
            foreach (var (l_ID, l_配列) in l_エントリ群)
            {
                if (!l_ID.Contains(Consts.環状の目印, StringComparison.OrdinalIgnoreCase)
                    || l_配列.Length < 接合窓長)
                {
                    continue;
                }

                // 環状なので末尾の続きは先頭になる
                // その繋ぎ目を跨ぐ窓を作る
                var l_窓 = string.Concat(l_配列.AsSpan(l_配列.Length - 接合フランク長), l_配列.AsSpan(0, 接合フランク長));
                if (!KmerPacking.TryGet_パック(l_窓, 0, 接合窓長, out var l_順鎖))
                {
                    continue;
                }
                l_接合窓[KmerPacking.Get_小さいほう(l_順鎖, 接合窓長)] = l_対象.Count;
                l_対象.Add((l_ID, l_配列.Length));
            }

            if (l_対象.Count == 0)
            {
                return [];
            }

            Logger.V_出力(メッセージID.閉じ目の検証開始, l_対象.Count, 接合窓長);

            var l_支持数 = new int[l_対象.Count];
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            ReadPipeline.V_実行(l_スレッド数, l_スレッド数 * 256, FastqReader.Get_生リード列(p_リード1のパス, p_リード2のパス), (l_リード, _) => V_集計_接合支持(l_リード, l_接合窓, l_支持数));

            List<環状閉鎖検証結果> l_結果 = [];
            for (var i = 0; i < l_対象.Count; i++)
            {
                l_結果.Add(new 環状閉鎖検証結果(l_対象[i].A_ID, l_対象[i].A_長さ, l_支持数[i], 閉じ目に必要なリード数));
            }
            return l_結果;
        }

        /// <summary>
        /// 環状閉鎖の検証結果をログへ出力する
        /// </summary>
        /// <param name="p_結果群">配列ごとの検証結果</param>
        public static void V_出力_検証結果(IReadOnlyList<環状閉鎖検証結果> p_結果群)
        {
            if (p_結果群.Count == 0)
            {
                Logger.V_出力(メッセージID.閉じ目_環状の配列が無い);
                return;
            }

            foreach (var l_結果 in p_結果群)
            {
                Logger.V_出力(l_結果.A_Has支持 ? メッセージID.閉じ目を裏付けた : メッセージID.閉じ目を裏付けられず, l_結果.A_配列ID, l_結果.A_長さ, l_結果.A_跨いだリード数, l_結果.A_必要本数);
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 本のリードが閉じ目の窓を含むかを調べ、含めばその配列の支持を 1 つ増やす
        /// </summary>
        /// <param name="p_リード">検査するリードの配列</param>
        /// <param name="p_接合窓">正規形の窓 -> 配列番号</param>
        /// <param name="p_支持数">配列番号ごとの支持数、見つかれば加算する</param>
        /// <remarks>
        /// 同じリードが同じ配列を何度支持しても 1 本と数える
        /// </remarks>
        private static void V_集計_接合支持(string p_リード, Dictionary<UInt128, int> p_接合窓, int[] p_支持数)
        {
            if (p_リード.Length < 接合窓長)
            {
                return;
            }

            var l_マスク = (UInt128.One << (接合窓長 << 1)) - 1;
            var l_最上位へ = (接合窓長 - 1) << 1;
            UInt128 l_順鎖 = 0;
            UInt128 l_逆鎖 = 0;
            var l_直近の曖昧位置 = -1;
            HashSet<int>? l_数えた配列 = null;

            for (var i = 0; i < p_リード.Length; i++)
            {
                var l_塩基ID = Util.Get_塩基ID(p_リード[i]);
                var l_Is有効 = l_塩基ID is >= Consts.塩基ID.A and <= Consts.塩基ID.T;
                var l_コドン = (UInt128)(l_Is有効 ? l_塩基ID - 1 : 0);
                l_順鎖 = ((l_順鎖 << 2) | l_コドン) & l_マスク;
                l_逆鎖 = (l_逆鎖 >> 2) | ((l_コドン ^ 3) << l_最上位へ);
                if (!l_Is有効)
                {
                    l_直近の曖昧位置 = i;
                }

                var l_開始 = i - 接合窓長 + 1;
                if (l_開始 < 0 || l_直近の曖昧位置 >= l_開始)
                {
                    continue;
                }

                var l_正規形 = l_順鎖 < l_逆鎖 ? l_順鎖 : l_逆鎖;
                if (!p_接合窓.TryGetValue(l_正規形, out var l_配列番号))
                {
                    continue;
                }

                // 1 本のリードは、同じ閉じ目を何度跨いで見えても 1 本の証拠でしかない
                l_数えた配列 ??= [];
                if (!l_数えた配列.Add(l_配列番号))
                {
                    continue;
                }
                _ = Interlocked.Increment(ref p_支持数[l_配列番号]);
            }
        }

        #endregion
    }
}
