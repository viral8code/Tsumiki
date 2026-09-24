using Tsumiki.Commons;
using Tsumiki.IO;
using Tsumiki.Models.Foundation;
using Tsumiki.Models.Preprocessing;
using Tsumiki.Models.Scaffolding;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Preprocessing
{
    /// <summary>
    /// ペアエンドの 2 本を、間の未読区間ごと 1 本の合成リード (SuperRead) に統合する
    /// </summary>
    internal static class SuperReadJoiner
    {
        #region 定数

        /// <summary>
        /// ペア結合に必要な最小重なり長
        /// </summary>
        private const int ペア結合の最小重なり長 = 15;

        /// <summary>
        /// ペア結合で許す不一致率
        /// </summary>
        private const double ペア結合の許容不一致率 = 0.05D;

        /// <summary>
        /// 橋渡しする長さの上限
        /// </summary>
        private const int 橋渡し長の上限 = 500;

        /// <summary>
        /// 書き出す合成リードに付けるクオリティ文字
        /// </summary>
        private const char 合成リードのクオリティ = 'I';

        /// <summary>
        /// 1 バッチあたりのペア数
        /// </summary>
        private const int バッチサイズ = 5_000;

        /// <summary>
        /// 1 ペアあたりに展開してよい探索状態の上限
        /// </summary>
        private const int 橋渡しの状態数上限 = 20_000;

        /// <summary>
        /// -i でインサートサイズが分かっているときに、そこから見積もった橋渡し長に掛ける許容比
        /// </summary>
        private const double インサートサイズの許容比 = 1.5D;

        #endregion

        #region 公開メソッド

        /// <summary>
        /// ペアの FASTQ を読み込み、統合できたペアを合成配列 (引き継ぎ配列) として返す
        /// </summary>
        /// <param name="p_リード1のパス">リード 1 のパス</param>
        /// <param name="p_リード2のパス">リード 2 のパス</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_統計">統合の内訳</param>
        /// <param name="p_推定断片長上限">-i が無いときに橋渡し長の上限を見積もる断片長、分からなければ null</param>
        /// <param name="p_推定断片長下限">重なりを探すオフセットの下限を決める断片長、分からなければ null</param>
        /// <param name="p_重なり結合の出力パス">渡すと、重なりで繋いだ断片だけを FASTQ として書き出す</param>
        /// <returns>橋渡しで繋いだ合成配列 (重なりで繋いだものは p_重なり結合の出力パス へ書き出し、引き継ぎには含めない)</returns>
        public static List<引き継ぎ配列> Get_合成リード(string p_リード1のパス, string p_リード2のパス, TrustedKmerIndex p_kmerインデックス, int p_k長, out SuperRead統計 p_統計, int? p_推定断片長上限 = null, int? p_推定断片長下限 = null, string? p_重なり結合の出力パス = null)
        {
            using var l_計測 = new StageTimer($"superread k={p_k長}");
            var l_スレッド数 = Math.Max(1, ConfigurationManager.A_実行時引数.A_スレッド数);
            var l_インサートサイズ = ConfigurationManager.A_実行時引数.A_インサートサイズ ?? p_推定断片長上限;

            var l_総ペア数 = 0;
            var l_統合数 = 0;
            var l_重なり結合数 = 0;
            var l_曖昧で捨てた数 = 0;
            var l_出力済みの区切り = 0UL;
            List<引き継ぎ配列> l_結果 = [];

            Logger.V_出力(メッセージID.SuperRead橋渡し開始);

            using var l_読み込み1 = new FastqReader(p_リード1のパス);
            using var l_読み込み2 = new FastqReader(p_リード2のパス);
            using var l_書き出し = p_重なり結合の出力パス is null ? null : new FastqWriter(p_重なり結合の出力パス);

            var l_配列1群 = new string[バッチサイズ];
            var l_配列2群 = new string[バッチサイズ];
            var l_統合結果群 = new (string? A_配列, bool A_Is重なり結合, int A_曖昧で捨てた数)[バッチサイズ];

            while (l_読み込み1.Has続き() && l_読み込み2.Has続き())
            {
                var l_件数 = 0;
                while (l_件数 < バッチサイズ && l_読み込み1.Has続き() && l_読み込み2.Has続き())
                {
                    l_配列1群[l_件数] = l_読み込み1.Get_次のリード_軽量().A_生リード;
                    l_配列2群[l_件数] = l_読み込み2.Get_次のリード_軽量().A_生リード;
                    l_件数++;
                }
                l_総ペア数 += l_件数;

                _ = Parallel.For(0, l_件数, new ParallelOptions { MaxDegreeOfParallelism = l_スレッド数 }, i =>
                {
                    l_統合結果群[i] = Get_合成配列_内訳つき(l_配列1群[i], l_配列2群[i], p_kmerインデックス, p_k長, l_インサートサイズ, p_推定断片長下限, p_推定断片長上限);
                });

                for (var i = 0; i < l_件数; i++)
                {
                    l_曖昧で捨てた数 += l_統合結果群[i].A_曖昧で捨てた数;
                    if (l_統合結果群[i].A_配列 is not { } l_配列)
                    {
                        continue;
                    }
                    l_統合数++;
                    if (l_統合結果群[i].A_Is重なり結合)
                    {
                        l_重なり結合数++;
                        l_書き出し?.V_書き込み(FormattableString.Invariant($"@F{l_重なり結合数}"), l_配列, new string(合成リードのクオリティ, l_配列.Length));
                        continue;
                    }

                    l_結果.Add(Get_引き継ぎ配列(l_配列, p_kmerインデックス, p_k長));
                }

                var l_区切り = (ulong)l_総ペア数 / Consts.進捗ログ間隔;
                if (l_区切り > l_出力済みの区切り)
                {
                    l_出力済みの区切り = l_区切り;
                    Logger.V_出力(メッセージID.SuperRead橋渡し進捗, l_総ペア数, l_統合数);
                }
            }

            p_統計 = new SuperRead統計(l_総ペア数, l_統合数, l_重なり結合数, l_曖昧で捨てた数);
            return l_結果;
        }

        /// <summary>
        /// 統合の内訳をログへ出力する
        /// </summary>
        /// <param name="p_統計">統合の内訳</param>
        public static void V_出力_統計(SuperRead統計 p_統計)
        {
            Logger.V_出力(メッセージID.SuperRead統計, p_統計.A_統合数, p_統計.A_総ペア数);
            Logger.V_出力(メッセージID.SuperRead統計_重なり, p_統計.A_重なり結合数, p_統計.A_橋渡し数);
            if (p_統計.A_曖昧で捨てた数 > 0)
            {
                Logger.V_出力(メッセージID.SuperRead統計_曖昧, p_統計.A_曖昧で捨てた数);
            }
        }

        #endregion

        #region テストメソッド

        /// <summary>
        /// 1 ペア分の統合を試みる
        /// </summary>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_インサートサイズ">-i で指定されたインサートサイズ、未指定なら null</param>
        /// <returns>合成配列、統合できなかった場合は null</returns>
        internal static string? Get_合成配列(string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_インサートサイズ = null)
        {
            return Get_合成配列_内訳つき(p_配列1, p_配列2, p_kmerインデックス, p_k長, p_インサートサイズ).A_配列;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 1 ペア分の統合を、どちらの手段で繋いだかと一緒に返す
        /// </summary>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_インサートサイズ">-i で指定されたインサートサイズ、未指定なら null</param>
        /// <param name="p_断片長下限">重なりを探す断片長の下限、分からなければ null</param>
        /// <param name="p_断片長上限">重なりを探す断片長の上限、分からなければ null</param>
        /// <returns>合成配列と、重なりで繋いだかどうか、重なりが曖昧で捨てた数</returns>
        internal static (string? A_配列, bool A_Is重なり結合, int A_曖昧で捨てた数) Get_合成配列_内訳つき(string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_インサートサイズ, int? p_断片長下限 = null, int? p_断片長上限 = null)
        {
            var l_曖昧 = 0;
            if (Get_重なりで結合(p_配列1, p_配列2, p_kmerインデックス, p_k長, ref l_曖昧, p_断片長下限, p_断片長上限) is { } l_重なり結合)
            {
                return (l_重なり結合, true, 0);
            }
            return (Get_橋渡しで結合(p_配列1, p_配列2, p_kmerインデックス, p_k長, p_インサートサイズ), false, l_曖昧);
        }

        /// <summary>
        /// read1 と RC (read2) の重なりから断片を復元する
        /// </summary>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_曖昧で捨てた数">重なりが一つに定まらず捨てた数、該当すれば 1 加算する</param>
        /// <param name="p_断片長下限">重なりを探す断片長の下限、分からなければ null</param>
        /// <param name="p_断片長上限">重なりを探す断片長の上限、分からなければ null</param>
        /// <returns>
        /// 復元した断片<br/>
        /// 重なりが見つからない (断片がリード長の 2 倍を超える) 場合と、繋いでも伸びない場合は null
        /// </returns>
        private static string? Get_重なりで結合(string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長, ref int p_曖昧で捨てた数, int? p_断片長下限 = null, int? p_断片長上限 = null)
        {
            var l_RC配列2 = Util.V_逆相補_曖昧塩基あり(p_配列2);

            var l_最小オフセット = p_断片長下限 is { } l_下限 ? l_下限 - p_配列2.Length : (int?)null;
            var l_最大オフセット = p_断片長上限 is { } l_上限 ? l_上限 - p_配列2.Length : (int?)null;
            var l_重なり = Preprocessor.Get_最適オーバーラップ(Util.V_変換_塩基列(p_配列1), Util.V_変換_塩基列(l_RC配列2), ペア結合の最小重なり長, ペア結合の許容不一致率, out var l_対抗馬があるか, l_最小オフセット, l_最大オフセット);
            if (l_重なり is not { } l_位置合わせ || l_位置合わせ.A_offset < 0)
            {
                return null;
            }

            if (l_対抗馬があるか)
            {
                p_曖昧で捨てた数++;
                return null;
            }

            var l_フラグメント長 = l_位置合わせ.A_offset + p_配列2.Length;
            if (l_フラグメント長 <= p_配列1.Length)
            {
                return null;
            }

            var l_合成 = p_配列1 + l_RC配列2[(p_配列1.Length - l_位置合わせ.A_offset)..];
            return l_合成.Length >= p_k長 && Has継ぎ目支持(l_合成, p_配列1.Length, p_kmerインデックス, p_k長)
                ? l_合成
                : null;
        }

        /// <summary>
        /// 合成配列の継ぎ目 (read1 から RC (read2) へ切り替わる位置) を跨ぐ k-mer が、信頼できる k-mer 集合に入っているか
        /// </summary>
        /// <param name="p_合成">継ぎ目を含む合成配列</param>
        /// <param name="p_継ぎ目">read1 から RC (read2) へ切り替わる位置</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <returns>継ぎ目を跨ぐ k-mer がすべて集合にあれば true</returns>
        private static bool Has継ぎ目支持(string p_合成, int p_継ぎ目, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_開始 = Math.Max(0, p_継ぎ目 - p_k長 + 1);
            var l_終了 = Math.Min(p_継ぎ目, p_合成.Length - p_k長);
            for (var l_位置 = l_開始; l_位置 <= l_終了; l_位置++)
            {
                var l_kmer = Util.V_変換_塩基列(p_合成.Substring(l_位置, p_k長));
                if (Array.IndexOf(l_kmer, Consts.無効な塩基) >= 0
                    || !p_kmerインデックス.Haskmer(l_kmer))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// read1 の末尾 k-mer から RC (read2) の先頭 k-mer まで、信頼できる k-mer 集合の中で経路を探して繋ぐ
        /// </summary>
        /// <param name="p_配列1">read1 の配列</param>
        /// <param name="p_配列2">read2 の配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <param name="p_インサートサイズ">-i で指定されたインサートサイズ、未指定なら null</param>
        /// <returns>
        /// 繋いだ配列<br/>
        /// どちらかの端の k-mer が集合に無い、あるいは経路が一意に定まらない場合は null
        /// </returns>
        private static string? Get_橋渡しで結合(string p_配列1, string p_配列2, TrustedKmerIndex p_kmerインデックス, int p_k長, int? p_インサートサイズ)
        {
            if (p_配列1.Length < p_k長 || p_配列2.Length < p_k長)
            {
                return null;
            }

            var l_最大長 = Get_橋渡し長上限(p_配列1.Length, p_配列2.Length, p_インサートサイズ);
            if (l_最大長 < 0)
            {
                return null;
            }

            var l_左のkmer = Util.V_変換_塩基列(p_配列1[^p_k長..]);
            if (Array.IndexOf(l_左のkmer, Consts.無効な塩基) >= 0
                || !p_kmerインデックス.Haskmer(l_左のkmer))
            {
                return null;
            }

            var l_目標kmer = Util.V_変換_塩基列(Util.V_逆相補_曖昧塩基あり(p_配列2[^p_k長..]));
            if (Array.IndexOf(l_目標kmer, Consts.無効な塩基) >= 0
                || !p_kmerインデックス.Haskmer(l_目標kmer))
            {
                return null;
            }

            if (Util.V_変換_塩基列(p_配列1).AsSpan().IndexOf(l_目標kmer) >= 0)
            {
                return null;
            }

            var (l_橋渡し配列, l_判定) = ConstrainedPathFinder.Get_経路(l_左のkmer, l_目標kmer, p_最小長: 0, p_最大長: l_最大長, p_kmerインデックス, p_k長, 橋渡しの状態数上限);
            return l_判定 == ギャップ充填判定.充填済み
                ? p_配列1 + l_橋渡し配列 + Util.V_逆相補_曖昧塩基あり(p_配列2)
                : null;
        }

        /// <summary>
        /// このペアで探索してよい橋渡し長の上限
        /// </summary>
        /// <param name="p_長さ1">read1 の長さ</param>
        /// <param name="p_長さ2">read2 の長さ</param>
        /// <param name="p_インサートサイズ">-i で指定されたインサートサイズ、未指定なら null</param>
        /// <returns>探索してよい橋渡し長の上限</returns>
        private static int Get_橋渡し長上限(int p_長さ1, int p_長さ2, int? p_インサートサイズ)
        {
            if (p_インサートサイズ is not { } l_インサートサイズ)
            {
                return 橋渡し長の上限;
            }
            var l_見積もり = (int)(l_インサートサイズ * インサートサイズの許容比) - p_長さ1 - p_長さ2;
            return Math.Min(橋渡し長の上限, l_見積もり);
        }

        /// <summary>
        /// 合成配列を、位置ごとのカバレッジ付きの引き継ぎ配列にする
        /// </summary>
        /// <param name="p_配列">合成配列</param>
        /// <param name="p_kmerインデックス">この k の信頼できる k-mer 集合</param>
        /// <param name="p_k長">この k の長さ</param>
        /// <returns>次の k へ渡す引き継ぎ配列</returns>
        private static 引き継ぎ配列 Get_引き継ぎ配列(string p_配列, TrustedKmerIndex p_kmerインデックス, int p_k長)
        {
            var l_塩基列 = p_配列.Select(Util.Get_塩基ID).ToArray();
            var l_カバレッジ = new int[p_配列.Length - p_k長 + 1];
            for (var i = 0; i < l_カバレッジ.Length; i++)
            {
                l_カバレッジ[i] = (int)Math.Min(int.MaxValue, p_kmerインデックス.Get_カバレッジ(l_塩基列.AsSpan(i, p_k長)));
            }
            return new 引き継ぎ配列(p_配列, l_カバレッジ, p_k長);
        }

        #endregion
    }
}
