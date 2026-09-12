using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Evaluation
{
    /// <summary>
    /// 出した配列の各位置について、そこを跨ぐ r-mer がリードに 1 度でも現れるかを調べる
    /// </summary>
    internal static class ReadSupportChecker
    {
        #region 定数

        /// <summary>
        /// 1 バッチあたりのリード数
        /// </summary>
        /// <remarks>
        /// まとめて読んで並列に照合する
        /// </remarks>
        private const int 照合バッチサイズ = 20_000;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 最終成果物がリードに裏付けられているかを調べる
        /// </summary>
        /// <param name="p_FASTAパス">調べる FASTA のパス</param>
        /// <param name="p_リード1のパス">リード 1 のパス</param>
        /// <param name="p_リード2のパス">リード 2 のパス、単一リードなら null</param>
        /// <param name="p_r長">支持を問う r-mer の長さ、2 bit パックの上限を超えられない</param>
        /// <returns>検査結果、調べられなかった場合は null</returns>
        public static 支持検査結果? Get_検査結果(string p_FASTAパス, string? p_リード1のパス, string? p_リード2のパス, int p_r長)
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

            var l_観測状態 = new byte[l_表.Count];
            V_記録_支持(l_表, l_観測状態, p_リード1のパス, p_リード2のパス, p_r長);

            return Get_集計(l_全件, l_位置ごとの番号, l_観測状態, p_r長);
        }

        /// <summary>
        /// 検査結果をログへ出力する
        /// </summary>
        /// <param name="p_結果">検査結果、調べられなかった場合は null</param>
        public static void V_出力_検査結果(支持検査結果? p_結果)
        {
            if (p_結果 is not { } l_結果)
            {
                return;
            }

            Logger.V_出力(メッセージID.支持検査の結果, l_結果.A_r長, l_結果.A_支持のない位置数, l_結果.A_調べた位置数, $"{l_結果.A_支持のない率:F4}", l_結果.A_区間.Count);
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// アセンブリの r-mer に通し番号を振った表と、位置ごとの番号を作る
        /// </summary>
        /// <remarks>
        /// 曖昧塩基を含む位置は -1 にして検査の対象から外す
        /// </remarks>
        /// <param name="p_全件">アセンブリの全配列</param>
        /// <param name="p_r長">支持を問う r-mer の長さ</param>
        /// <returns>r-mer から通し番号への表と、配列 ID ごとの位置別番号</returns>
        private static (Dictionary<UInt128, int> A_表, Dictionary<string, int[]> A_位置ごとの番号) Get_照合表(IReadOnlyList<(string A_ID, string A_配列)> p_全件, int p_r長)
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
                    if (!KmerPacking.TryGet_正規化キー(l_配列, i, p_r長, out var l_正規形))
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
        /// リードを流して、表にある r-mer に印を付ける
        /// </summary>
        /// <remarks>
        /// 書き込みは 0 でなくするだけなので、並列に走っても取りこぼしは出ない
        /// </remarks>
        /// <param name="p_表">r-mer から通し番号への表</param>
        /// <param name="p_観測状態">通し番号ごとの観測状態</param>
        /// <param name="p_リード1のパス">リード 1 のパス</param>
        /// <param name="p_リード2のパス">リード 2 のパス、単一リードなら null</param>
        /// <param name="p_r長">支持を問う r-mer の長さ</param>
        private static void V_記録_支持(Dictionary<UInt128, int> p_表, byte[] p_観測状態, string? p_リード1のパス, string? p_リード2のパス, int p_r長)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_バッチ = new string[照合バッチサイズ];
            var l_件数 = 0;

            foreach (var l_リード in FastqReader.Get_生リード列(p_リード1のパス, p_リード2のパス))
            {
                l_バッチ[l_件数++] = l_リード;
                if (l_件数 == 照合バッチサイズ)
                {
                    V_照合(p_表, p_観測状態, l_バッチ, l_件数, p_r長, l_スレッド数);
                    l_件数 = 0;
                }
            }
            if (l_件数 > 0)
            {
                V_照合(p_表, p_観測状態, l_バッチ, l_件数, p_r長, l_スレッド数);
            }
        }

        /// <summary>
        /// バッチに溜めたリードを並列に照合し、表にある r-mer に印を付ける
        /// </summary>
        /// <param name="p_表">r-mer から通し番号への表</param>
        /// <param name="p_観測状態">通し番号ごとの観測状態</param>
        /// <param name="p_バッチ">照合対象のリードを溜めた配列</param>
        /// <param name="p_件数">バッチに溜まっている件数</param>
        /// <param name="p_r長">支持を問う r-mer の長さ</param>
        /// <param name="p_スレッド数">並列度</param>
        private static void V_照合(Dictionary<UInt128, int> p_表, byte[] p_観測状態, string[] p_バッチ, int p_件数, int p_r長, int p_スレッド数)
        {
            _ = Parallel.For(0, p_件数, new ParallelOptions { MaxDegreeOfParallelism = p_スレッド数 }, i =>
            {
                var l_リード = p_バッチ[i];
                for (var p = 0; p + p_r長 <= l_リード.Length; p++)
                {
                    if (KmerPacking.TryGet_正規化キー(l_リード, p, p_r長, out var l_正規形)
                        && p_表.TryGetValue(l_正規形, out var l_番号))
                    {
                        p_観測状態[l_番号] = 1;
                    }
                }
            });
        }

        /// <summary>
        /// 印の付かなかった位置を数え、連続しているものを一つの区間にまとめる
        /// </summary>
        /// <param name="p_全件">アセンブリの全配列</param>
        /// <param name="p_位置ごとの番号">配列 ID ごとの位置別番号</param>
        /// <param name="p_観測状態">通し番号ごとの観測状態</param>
        /// <param name="p_r長">支持を問う r-mer の長さ</param>
        /// <returns>検査結果</returns>
        private static 支持検査結果 Get_集計(IReadOnlyList<(string A_ID, string A_配列)> p_全件, Dictionary<string, int[]> p_位置ごとの番号, byte[] p_観測状態, int p_r長)
        {
            var l_調べた = 0L;
            var l_支持なし = 0L;
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
                        // 曖昧塩基を含む位置にはギャップの N が来るので、支持が無いのではなく問えないものとして飛ばす
                        V_閉じる(l_区間, l_ID, ref l_開始, i - 1, p_r長);
                        continue;
                    }

                    l_調べた++;
                    if (p_観測状態[l_番号] != 0)
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
        /// 続いていた支持なしの並びを区間として確定する
        /// </summary>
        /// <remarks>
        /// 区間は r-mer の開始位置ではなく、その r-mer が覆う塩基の範囲で表す
        /// </remarks>
        /// <param name="p_区間">確定した区間の集まり</param>
        /// <param name="p_ID">配列 ID</param>
        /// <param name="p_開始">続いている並びの開始位置、続いていなければ -1</param>
        /// <param name="p_終わり">続いている並びの終了位置</param>
        /// <param name="p_r長">支持を問う r-mer の長さ</param>
        private static void V_閉じる(List<支持のない区間> p_区間, string p_ID, ref int p_開始, int p_終わり, int p_r長)
        {
            if (p_開始 < 0)
            {
                return;
            }
            p_区間.Add(new 支持のない区間(p_ID, p_開始 + 1, p_終わり + p_r長));
            p_開始 = -1;
        }

        #endregion
    }
}
