using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model.Evaluation;
using Tsumiki.Model.Foundation;
using Tsumiki.Utility;

namespace Tsumiki.Core.Evaluation
{
    /// <summary>
    /// 出した配列の各位置について、そこを跨ぐ r-mer がリードに1度でも
    /// 現れるかを調べる。
    ///
    /// 誤って繋いだ接合は、両側それぞれは正しい配列なので、カバレッジや
    /// コピー数といった局所の量では正当に見える。区別できるのは接合を
    /// 跨ぐ証拠だけで、繋ぎ目を含む r-mer がどのリードにも無いことが
    /// そのまま「この繋ぎ方を見た読みは一つも無い」を意味する。
    /// マッピングによる深度がこれを取り逃がすのは、マッパーが不一致を
    /// 許すため誤った接合を跨いでリードが載ってしまうから。完全一致で
    /// 問う必要がある。
    ///
    /// 表はアセンブリ側に持ち、リードを流して印を付ける。リード側の
    /// k-mer 集合を作るとエラー由来の種類数がゲノムの数十倍に膨らむが、
    /// この向きならアセンブリの長さぶんで済む。
    /// </summary>
    internal static class ReadSupportChecker
    {
        /// <summary>
        /// 1バッチあたりのリード数。まとめて読んで並列に照合する。
        /// </summary>
        private const int 照合バッチサイズ = 20000;

        /// <summary>
        /// 検査結果を返す。r 長は 2bit パックの上限を超えられない。
        /// </summary>
        public static 支持検査結果? Get_検査結果(
            string p_FASTAパス, string? p_リード1のパス, string? p_リード2のパス, int p_r長)
        {
            if (p_r長 is < 1 or > 64 || !File.Exists(p_FASTAパス))
            {
                return null;
            }

            var l_全件 = FastaReader.Get_全エントリ(p_FASTAパス);
            var (l_表, l_位置ごとの番号) = Get_照合表(l_全件, p_r長);
            if (l_表.Count == 0)
            {
                return null;
            }

            var l_見たか = new byte[l_表.Count];
            V_印を付ける(l_表, l_見たか, p_リード1のパス, p_リード2のパス, p_r長);

            return Get_集計(l_全件, l_位置ごとの番号, l_見たか, p_r長);
        }

        /// <summary>
        /// アセンブリの r-mer に通し番号を振った表と、位置ごとの番号を作る。
        /// 曖昧塩基を含む位置は -1 にして検査の対象から外す。
        /// </summary>
        private static (Dictionary<UInt128, int> A_表, Dictionary<string, int[]> A_位置ごとの番号) Get_照合表(
            IReadOnlyList<(string A_ID, string A_配列)> p_全件, int p_r長)
        {
            Dictionary<UInt128, int> l_表 = [];
            Dictionary<string, int[]> l_位置ごと = [];

            foreach (var (l_ID, l_配列) in p_全件)
            {
                var l_数 = l_配列.Length - p_r長 + 1;
                if (l_数 <= 0)
                {
                    l_位置ごと[l_ID] = [];
                    continue;
                }

                var l_番号列 = new int[l_数];
                for (var i = 0; i < l_数; i++)
                {
                    if (!KmerPacking.Get_正規化パック(l_配列, i, p_r長, out var l_正規形))
                    {
                        l_番号列[i] = -1;
                        continue;
                    }
                    if (!l_表.TryGetValue(l_正規形, out var l_番号))
                    {
                        l_番号 = l_表.Count;
                        l_表[l_正規形] = l_番号;
                    }
                    l_番号列[i] = l_番号;
                }
                l_位置ごと[l_ID] = l_番号列;
            }

            return (l_表, l_位置ごと);
        }

        /// <summary>
        /// リードを流して、表にある r-mer に印を付ける。
        /// 書き込みは「0 でなくする」だけなので、並列に走っても取りこぼしは出ない。
        /// </summary>
        private static void V_印を付ける(
            Dictionary<UInt128, int> p_表, byte[] p_見たか,
            string? p_リード1のパス, string? p_リード2のパス, int p_r長)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_バッチ = new string[照合バッチサイズ];
            var l_件数 = 0;

            void V_照合()
            {
                var l_今回 = l_件数;
                _ = Parallel.For(0, l_今回, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, i =>
                {
                    var l_リード = l_バッチ[i];
                    for (var p = 0; p + p_r長 <= l_リード.Length; p++)
                    {
                        if (KmerPacking.Get_正規化パック(l_リード, p, p_r長, out var l_正規形)
                            && p_表.TryGetValue(l_正規形, out var l_番号))
                        {
                            p_見たか[l_番号] = 1;
                        }
                    }
                });
                l_件数 = 0;
            }

            foreach (var l_リード in FastqReader.Get_生リード列(p_リード1のパス, p_リード2のパス))
            {
                l_バッチ[l_件数++] = l_リード;
                if (l_件数 == 照合バッチサイズ)
                {
                    V_照合();
                }
            }
            if (l_件数 > 0)
            {
                V_照合();
            }
        }

        /// <summary>
        /// 印の付かなかった位置を数え、連続しているものを一つの区間にまとめる。
        /// </summary>
        private static 支持検査結果 Get_集計(
            IReadOnlyList<(string A_ID, string A_配列)> p_全件,
            Dictionary<string, int[]> p_位置ごとの番号,
            byte[] p_見たか,
            int p_r長)
        {
            long l_調べた = 0;
            long l_支持なし = 0;
            List<支持のない区間> l_区間 = [];

            foreach (var (l_ID, _) in p_全件)
            {
                var l_番号列 = p_位置ごとの番号[l_ID];
                var l_開始 = -1;
                for (var i = 0; i < l_番号列.Length; i++)
                {
                    var l_番号 = l_番号列[i];
                    if (l_番号 < 0)
                    {
                        // 曖昧塩基を含む位置。ギャップの N がここに来るので、
                        // 支持が無いのではなく問えないものとして飛ばす。
                        V_閉じる(l_区間, l_ID, ref l_開始, i - 1, p_r長);
                        continue;
                    }

                    l_調べた++;
                    if (p_見たか[l_番号] != 0)
                    {
                        V_閉じる(l_区間, l_ID, ref l_開始, i - 1, p_r長);
                        continue;
                    }

                    l_支持なし++;
                    if (l_開始 < 0)
                    {
                        l_開始 = i;
                    }
                }
                V_閉じる(l_区間, l_ID, ref l_開始, l_番号列.Length - 1, p_r長);
            }

            return new 支持検査結果(p_r長, l_調べた, l_支持なし, l_区間);
        }

        /// <summary>
        /// 続いていた支持なしの並びを区間として確定する。
        /// 区間は r-mer の開始位置ではなく、その r-mer が覆う塩基の範囲で表す。
        /// </summary>
        private static void V_閉じる(
            List<支持のない区間> p_区間, string p_ID, ref int p_開始, int p_終わり, int p_r長)
        {
            if (p_開始 < 0)
            {
                return;
            }
            p_区間.Add(new 支持のない区間(p_ID, p_開始 + 1, p_終わり + p_r長));
            p_開始 = -1;
        }

        public static void V_出力_検査結果(支持検査結果? p_結果)
        {
            if (p_結果 is not { } l_結果)
            {
                return;
            }
            Logger.V_出力(
                メッセージID.支持検査の結果,
                l_結果.A_r長,
                l_結果.A_支持のない位置数,
                l_結果.A_調べた位置数,
                $"{l_結果.A_支持のない率:F4}",
                l_結果.A_区間.Count);
        }
    }
}
