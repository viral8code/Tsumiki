using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Preprocessing
{
    /// <summary>
    /// 前段の k で組み上がった配列を、次の k の k-mer 集合へ引き継ぐ
    /// </summary>
    internal static class KmerCarryOver
    {
        #region 定数

        /// <summary>
        /// 引き継ぐ配列の最小長
        /// </summary>
        /// <remarks>
        /// 前段で短く切れた断片は連結の役に立たないうえ、エラー由来の残骸である可能性が相対的に高い
        /// </remarks>
        private const int 引き継ぐ配列の最小長 = 500;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 引き継ぎ元の配列とカバレッジを、その k の成果物から作る
        /// </summary>
        /// <param name="p_FASTAパス"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <remarks>
        /// k-mer インデックスが生きているうちにしか作れない
        /// </remarks>
        /// <returns></returns>
        public static List<引き継ぎ配列> Get_引き継ぎ配列(string p_FASTAパス, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            List<引き継ぎ配列> l_結果 = [];
            using var l_読み込み = new FastaReader(p_FASTAパス);

            while (l_読み込み.Has続き())
            {
                var l_配列 = l_読み込み.Get_次の配列().A_配列;
                if (l_配列.Length < Math.Max(引き継ぐ配列の最小長, p_k長))
                {
                    continue;
                }

                var l_塩基列 = l_配列.Select(Util.Get_塩基ID).ToArray();
                var l_カバレッジ = new int[l_配列.Length - p_k長 + 1];
                for (var i = 0; i < l_カバレッジ.Length; i++)
                {
                    l_カバレッジ[i] = (int)Math.Min(int.MaxValue, p_kmerインデックス.Get_カバレッジ(l_塩基列.AsSpan(i, p_k長)));
                }
                l_結果.Add(new 引き継ぎ配列(l_配列, l_カバレッジ, p_k長, A_Is確定経路: true));
            }
            return l_結果;
        }

        /// <summary>
        /// 引き継ぎ配列のうち、この k の集合に無い k-mer を足す
        /// </summary>
        /// <param name="p_引き継ぎ配列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 既にある k-mer は触らない (実際のリード由来の観測を優先する) <br/>
        /// 戻り値は足した k-mer の数
        /// </remarks>
        /// <returns></returns>
        public static int V_引き継ぎ(IReadOnlyList<引き継ぎ配列> p_引き継ぎ配列, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長)
        {
            using var l_計測 = new StageTimer($"carry-over k={p_k長}");
            return p_k長 <= TrustedKmerIndex.パック値のk上限
                ? Get_追加数_パック値(p_引き継ぎ配列, p_kmerインデックス, p_k長, p_リード長)
                : Get_追加数_塩基列(p_引き継ぎ配列, p_kmerインデックス, p_k長, p_リード長);
        }

        /// <summary>
        /// この k-mer に与えるカバレッジ
        /// </summary>
        /// <param name="p_引き継ぎ"></param>
        /// <param name="p_位置"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <returns></returns>
        public static ulong Get_引き継ぐカバレッジ(引き継ぎ配列 p_引き継ぎ, int p_位置, int p_k長, int? p_リード長)
        {
            var l_終端 = Math.Min(p_引き継ぎ.A_カバレッジ.Length - 1, p_位置 + p_k長 - p_引き継ぎ.A_k長);
            if (p_位置 > l_終端)
            {
                return 0UL;
            }

            var l_最小 = int.MaxValue;
            for (var i = p_位置; i <= l_終端; i++)
            {
                l_最小 = Math.Min(l_最小, p_引き継ぎ.A_カバレッジ[i]);
            }
            return Get_換算カバレッジ(l_最小, p_引き継ぎ.A_k長, p_k長, p_リード長);
        }

        /// <summary>各窓の最小値を単調な待ち行列で求める</summary>
        /// <param name="p_値">元のカバレッジ</param>
        /// <param name="p_幅">窓の幅</param>
        /// <param name="p_結果">位置ごとの最小値</param>
        /// <param name="p_待ち行列">元のカバレッジ以上の長さの作業領域</param>
        internal static void V_計算_最小値列(int[] p_値, int p_幅, Span<int> p_結果, Span<int> p_待ち行列)
        {
            var l_先頭 = 0;
            var l_末尾 = 0;
            var l_読込位置 = 0;
            for (var i = 0; i < p_結果.Length; i++)
            {
                while (l_先頭 < l_末尾 && p_待ち行列[l_先頭] < i)
                {
                    l_先頭++;
                }
                var l_終端 = Math.Min(p_値.Length, i + p_幅);
                while (l_読込位置 < l_終端)
                {
                    while (l_先頭 < l_末尾 && p_値[p_待ち行列[l_末尾 - 1]] >= p_値[l_読込位置])
                    {
                        l_末尾--;
                    }
                    p_待ち行列[l_末尾++] = l_読込位置++;
                }
                p_結果[i] = p_幅 <= 0 || i >= p_値.Length || l_先頭 == l_末尾 ? 0 : p_値[p_待ち行列[l_先頭]];
            }
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// パック値で窓を転がし、集合に無い k-mer を足す (k &lt;= 128)
        /// </summary>
        /// <param name="p_引き継ぎ配列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <remarks>
        /// 集合に有るかの判定は読むだけなので配列ごとに並列に行い、足すのは配列の順に 1 本のスレッドで行う<br/>
        /// 同じ k-mer を複数の配列が持つとき先に来た配列のカバレッジが残る点は、逐次に足した場合と変わらない
        /// </remarks>
        /// <returns>足した k-mer の数</returns>
        private static int Get_追加数_パック値(IReadOnlyList<引き継ぎ配列> p_引き継ぎ配列, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長)
        {
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_バッチ長 = (int)Consts.進捗ログ間隔;
            var l_未登録群 = new List<(UInt128 A_上位, UInt128 A_下位, ulong A_カバレッジ)>?[Math.Min(l_バッチ長, p_引き継ぎ配列.Count)];
            var l_追加数 = 0;

            for (var l_開始 = 0; l_開始 < p_引き継ぎ配列.Count; l_開始 += l_バッチ長)
            {
                var l_件数 = Math.Min(l_バッチ長, p_引き継ぎ配列.Count - l_開始);
                _ = Parallel.For(0, l_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, () => (A_最小値列: Array.Empty<int>(), A_待ち行列: Array.Empty<int>()), (i, _, l_作業域) =>
                {
                    l_未登録群[i] = Get_未登録kmer(p_引き継ぎ配列[l_開始 + i], p_kmerインデックス, p_k長, p_リード長, ref l_作業域.A_最小値列, ref l_作業域.A_待ち行列);
                    return l_作業域;
                }, _ => { });

                for (var i = 0; i < l_件数; i++)
                {
                    if (l_未登録群[i] is not { } l_未登録)
                    {
                        continue;
                    }
                    foreach (var (l_上位, l_下位, l_カバレッジ) in l_未登録)
                    {
                        if (p_kmerインデックス.Try追加_信頼kmer_正規形(l_上位, l_下位, l_カバレッジ))
                        {
                            l_追加数++;
                        }
                    }
                    l_未登録群[i] = null;
                }

                // -sr を使うと引き継ぎ配列はリードペアの数まで増える
                // 全 k-mer を足し終えるまで無言だと止まったように見える
                if (l_件数 == l_バッチ長)
                {
                    Logger.V_出力(メッセージID.引き継ぎの統合進捗, l_開始 + l_件数, p_引き継ぎ配列.Count);
                }
            }
            return l_追加数;
        }

        /// <summary>
        /// 引き継ぎ配列 1 本のうち、集合に無く、与えるカバレッジが正の k-mer を配列の順に返す (k &lt;= 128)
        /// </summary>
        /// <param name="p_引き継ぎ"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <param name="p_最小値列">窓ごとの最小カバレッジの作業領域、足りなければ確保し直す</param>
        /// <param name="p_待ち行列">最小値を求める待ち行列の作業領域、足りなければ確保し直す</param>
        /// <returns>1 つも無ければ null</returns>
        private static List<(UInt128 A_上位, UInt128 A_下位, ulong A_カバレッジ)>? Get_未登録kmer(引き継ぎ配列 p_引き継ぎ, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長, ref int[] p_最小値列, ref int[] p_待ち行列)
        {
            if (p_引き継ぎ.A_配列.Length < p_k長)
            {
                return null;
            }

            var l_窓数 = p_引き継ぎ.A_配列.Length - p_k長 + 1;
            if (p_最小値列.Length < l_窓数)
            {
                p_最小値列 = new int[l_窓数];
            }
            if (p_待ち行列.Length < p_引き継ぎ.A_カバレッジ.Length)
            {
                p_待ち行列 = new int[p_引き継ぎ.A_カバレッジ.Length];
            }
            V_計算_最小値列(p_引き継ぎ.A_カバレッジ, p_k長 - p_引き継ぎ.A_k長 + 1, p_最小値列.AsSpan(0, l_窓数), p_待ち行列);

            List<(UInt128 A_上位, UInt128 A_下位, ulong A_カバレッジ)>? l_未登録 = null;
            var l_窓 = new RollingKmer(p_k長);
            for (var i = 0; i < p_引き継ぎ.A_配列.Length; i++)
            {
                if (!l_窓.Try追加(p_引き継ぎ.A_配列[i], out var l_キー) || p_kmerインデックス.Haskmer_正規形(l_キー.A_上位, l_キー.A_下位))
                {
                    continue;
                }
                var l_カバレッジ = Get_換算カバレッジ(p_最小値列[i - p_k長 + 1], p_引き継ぎ.A_k長, p_k長, p_リード長);
                if (l_カバレッジ > 0UL)
                {
                    (l_未登録 ??= []).Add((l_キー.A_上位, l_キー.A_下位, l_カバレッジ));
                }
            }
            return l_未登録;
        }

        /// <summary>
        /// 塩基列のまま窓を見て、集合に無い k-mer を足す (k &gt; 128)
        /// </summary>
        /// <param name="p_引き継ぎ配列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_リード長"></param>
        /// <returns>足した k-mer の数</returns>
        private static int Get_追加数_塩基列(IReadOnlyList<引き継ぎ配列> p_引き継ぎ配列, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_リード長)
        {
            var l_追加数 = 0;
            var l_処理数 = 0;
            var l_出力済みの区切り = 0UL;
            int[] l_最小値列 = [];
            int[] l_待ち行列 = [];
            foreach (var l_引き継ぎ in p_引き継ぎ配列)
            {
                l_処理数++;
                var l_区切り = (ulong)l_処理数 / Consts.進捗ログ間隔;
                if (l_区切り > l_出力済みの区切り)
                {
                    l_出力済みの区切り = l_区切り;
                    Logger.V_出力(メッセージID.引き継ぎの統合進捗, l_処理数, p_引き継ぎ配列.Count);
                }

                if (l_引き継ぎ.A_配列.Length < p_k長)
                {
                    continue;
                }

                // Array.IndexOf を窓ごとに呼ぶと窓 1 つあたり O (k) かかり、
                // 配列全体では O (n*k) になる
                // 窓をスライドさせる際に新しく入る 1 塩基だけを見て
                // 「直近に見た無効塩基の位置」を更新すれば、
                // その位置が現在の窓の左端以降にある間は判定を使い回せる
                // (無効塩基は稀なので償却 O (n) で済む)
                var l_塩基列 = Util.V_変換_塩基列(l_引き継ぎ.A_配列);
                var l_窓数 = l_塩基列.Length - p_k長 + 1;
                if (l_最小値列.Length < l_窓数)
                {
                    l_最小値列 = new int[l_窓数];
                }
                if (l_待ち行列.Length < l_引き継ぎ.A_カバレッジ.Length)
                {
                    l_待ち行列 = new int[l_引き継ぎ.A_カバレッジ.Length];
                }
                V_計算_最小値列(l_引き継ぎ.A_カバレッジ, p_k長 - l_引き継ぎ.A_k長 + 1, l_最小値列.AsSpan(0, l_窓数), l_待ち行列);

                var l_直近の無効塩基位置 = -1;
                for (var i = 0; i + p_k長 <= l_塩基列.Length; i++)
                {
                    var l_新規末尾 = i + p_k長 - 1;
                    if (i == 0)
                    {
                        for (var j = 0; j < p_k長; j++)
                        {
                            if (l_塩基列[j] == Consts.無効な塩基)
                            {
                                l_直近の無効塩基位置 = j;
                            }
                        }
                    }
                    else if (l_塩基列[l_新規末尾] == Consts.無効な塩基)
                    {
                        l_直近の無効塩基位置 = l_新規末尾;
                    }

                    if (l_直近の無効塩基位置 >= i)
                    {
                        continue;
                    }

                    var l_カバレッジ = Get_換算カバレッジ(l_最小値列[i], l_引き継ぎ.A_k長, p_k長, p_リード長);
                    if (l_カバレッジ > 0UL
                        && p_kmerインデックス.Try追加_信頼kmer(l_塩基列.AsSpan(i, p_k長), l_カバレッジ))
                    {
                        l_追加数++;
                    }
                }
            }
            return l_追加数;
        }

        /// <summary>最小カバレッジを次の k の観測本数へ換算する</summary>
        /// <param name="p_最小">最小カバレッジ</param>
        /// <param name="p_前k長">前段の k</param>
        /// <param name="p_k長">次段の k</param>
        /// <param name="p_リード長">リード長</param>
        /// <returns>換算したカバレッジ</returns>
        private static ulong Get_換算カバレッジ(int p_最小, int p_前k長, int p_k長, int? p_リード長)
        {
            if (p_最小 <= 0)
            {
                return 0UL;
            }

            if (p_リード長 is not { } l_リード長 || l_リード長 <= p_k長)
            {
                return (ulong)p_最小;
            }

            var l_前段の本数 = l_リード長 - p_前k長 + 1;
            var l_今の本数 = l_リード長 - p_k長 + 1;
            return l_前段の本数 <= 0 || l_今の本数 <= 0 ? (ulong)p_最小 : (ulong)Math.Max(1L, (long)Math.Round((double)p_最小 * l_今の本数 / l_前段の本数));
        }

        #endregion
    }
}
