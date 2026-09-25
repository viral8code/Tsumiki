using Tsumiki.Commons;
using Tsumiki.Core;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Cores.Preprocessing;
using Tsumiki.Cores.Scaffolding;
using Tsumiki.Cores.UnitigBuilding;
using Tsumiki.IO;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Foundation;
using Tsumiki.Utilities;

namespace Tsumiki.Cores.Pipeline
{
    /// <summary>
    /// 指定した k 長で、k-mer カウントから scaffold までを一通り実行する
    /// </summary>
    internal static class AssemblyPipeline
    {
        #region 定数

        /// <summary>
        /// r-mer 長を k より長くする量
        /// </summary>
        private const int C_rMer長のk超過分の既定値 = 10;

        /// <summary>
        /// r-mer 検証に必要な窓数
        /// </summary>
        private const int C_rMer検証に必要な窓数 = 20;

        /// <summary>
        /// unitig ファイル名
        /// </summary>
        private const string C_Unitigファイル名 = "unitigs.fasta";

        /// <summary>
        /// contig ファイル名
        /// </summary>
        internal const string C_Contigファイル名 = "contigs.fasta";

        /// <summary>
        /// 最終アセンブリファイル名
        /// </summary>
        private const string C_最終アセンブリファイル名 = "assembly.fasta";

        /// <summary>
        /// unitig 数の上限
        /// </summary>
        private const int C_Unitig数の上限 = 100_000;

        /// <summary>
        /// 合成リードの橋渡し長の上限を見積もる断片長の分位
        /// </summary>
        private const double C_橋渡しに使う断片長の分位 = 0.99D;

        /// <summary>
        /// 重なりを探すオフセットの下限を決める断片長の分位
        /// </summary>
        private const double C_重なりに使う断片長の分位 = 0.01D;

        /// <summary>
        /// 重なりで繋いだ断片を書き出すファイル名
        /// </summary>
        private const string C_合成リードファイル名 = "fragments.fq";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// p_k長 でアセンブリを実行し、生成物のパスを返す
        /// </summary>
        /// <param name="p_引数"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_一時ディレクトリ"></param>
        /// <param name="p_リード長"></param>
        /// <param name="p_引き継ぎ"></param>
        /// <param name="p_次への引き継ぎ"></param>
        /// <param name="p_合成リードの控え"></param>
        /// <param name="p_原入力">反復検査に使う加工前の入力</param>
        /// <returns></returns>
        public static アセンブリ実行結果? Get_実行結果(Parameters p_引数, int p_k長, string p_一時ディレクトリ, int? p_リード長, IReadOnlyList<引き継ぎ配列>? p_引き継ぎ = null, List<引き継ぎ配列>? p_次への引き継ぎ = null, List<引き継ぎ配列>? p_合成リードの控え = null, Parameters? p_原入力 = null)
        {
            using var l_計測 = new StageTimer($"assembly k={p_k長}");
            if (p_引数.A_k長 != p_k長)
            {
                p_引数.Set_推定k長(p_k長);
            }

            var l_作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"k{p_k長}");
            _ = Directory.CreateDirectory(l_作業ディレクトリ);
            var l_断片パス = Path.Combine(p_一時ディレクトリ, C_合成リードファイル名);

            AmbiguityRecorder.V_開始(p_k長);

            var l_unitigパス = Path.Combine(l_作業ディレクトリ, C_Unitigファイル名);
            var l_contigパス = Path.Combine(l_作業ディレクトリ, C_Contigファイル名);
            var l_scaffoldパス = Path.Combine(l_作業ディレクトリ, Consts.Scaffoldファイル名);
            var l_GFAパス = Path.Combine(l_作業ディレクトリ, Consts.GFAファイル名);

            if (File.Exists(l_scaffoldパス))
            {
                File.Delete(l_scaffoldパス);
            }

            Logger.V_出力(メッセージID.kmerインデックス構築開始);
            using var l_kmerインデックス = new TrustedKmerIndex(l_作業ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_kmerインデックス;

            KmerCounting.V_読込_リードペア(p_引数, l_kmerインデックス, p_Is進行状況出力: true);

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.kmerカットオフ適用);
            using (new StageTimer($"cutoff k={p_k長}"))
            {
                KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_kmerインデックス);

                l_kmerインデックス.V_適用_カットオフ(p_引数.A_kmerカットオフ, p_引数.A_Is救済kmer使用 ? LowCoverageBridger.C_控えの最小出現回数 : 0UL);
            }

            KmerHistogram.V_出力_スペクトル(l_kmerインデックス.A_出現回数ヒストグラム, p_k長, p_リード長);

            if (p_引数.A_Is救済kmer使用)
            {
                Logger.V_出力(メッセージID.救済kmerの開始);
                _ = MercyKmerRescuer.Get_救済数(p_引数, l_kmerインデックス, p_k長);
            }

            if (p_引き継ぎ is { Count: > 0 })
            {
                Logger.V_出力(メッセージID.引き継ぎの統合開始, p_引き継ぎ.Count);
                var l_追加数 = KmerCarryOver.V_引き継ぎ(p_引き継ぎ, l_kmerインデックス, p_k長, p_リード長);
                Logger.V_出力(メッセージID.引き継ぎで追加したkmer数, l_追加数);
            }

            if (p_引数.A_Is救済kmer使用)
            {
                Logger.V_出力(メッセージID.低カバレッジ架橋開始);
                _ = LowCoverageBridger.Get_架橋kmer数(l_kmerインデックス, p_k長, p_リード長);
            }

            l_kmerインデックス.V_解放_控え();

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.tip除去開始);

            List<byte[]> l_開始kmer;
            using (new StageTimer($"graph-simplify k={p_k長}"))
            {
                l_開始kmer = GraphSimplifier.V_除去_tip(l_kmerインデックス, p_k長, p_リード長, p_Is低カバレッジ端トリミング: p_引数.A_Is低カバレッジ端トリミング);
            }

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.Unitig構築開始);
            var l_unitig配列 = Get_Unitig(l_kmerインデックス, l_開始kmer, p_k長, l_unitigパス, out var l_上限に達したか);

            AssemblyStatsReporter.V_出力_統計("unitigs", l_unitigパス);

            if (l_上限に達したか)
            {
                Logger.V_出力(メッセージID.グラフが複雑すぎる, p_k長, C_Unitig数の上限);
                return null;
            }

            Logger.V_出力(メッセージID.リードのマッピング開始);
            var l_contig構築 = new ContigMaker(l_unitigパス);

            var l_unitig長 = l_unitig配列.ToDictionary(x => x.Key, x => x.Value.Length);
            var l_カバレッジ = CopyNumberEstimator.Get_カバレッジ(l_kmerインデックス, l_unitig配列, p_k長);
            var l_グラフ = l_contig構築.Get_グラフ();
            var l_コピー数推定 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_unitig長, l_グラフ, p_引数.A_コピー数基準の出所);
            CopyNumberEstimator.V_出力_推定結果(l_コピー数推定, l_unitig長);

            Logger.V_出力_タイムスタンプ();

            using (new StageTimer($"read-mapping k={p_k長}"))
            {
                for (var i = 0; i < p_引数.A_ライブラリ数; i++)
                {
                    var (A_リード1, A_リード2) = p_引数.A_ライブラリ群[i];
                    Logger.V_出力(メッセージID.リードファイルのパス, A_リード1);
                    if (string.IsNullOrWhiteSpace(A_リード2))
                    {
                        l_contig構築.V_マッピング_リード(A_リード1);
                    }
                    else
                    {
                        Logger.V_出力(メッセージID.リードファイルのパス, A_リード2);
                        l_contig構築.V_マッピング_ペアリード(A_リード1, A_リード2, i);
                    }
                }

                foreach (var l_パス in Get_断片パス群(l_断片パス, p_引数.A_ライブラリ数).Where(中間データ置き場.Is存在))
                {
                    Logger.V_出力(メッセージID.リードファイルのパス, l_パス);
                    l_contig構築.V_マッピング_リード(l_パス);
                }
            }

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.Unitig結合開始);

            List<string> l_バブル敗者 = [];

            RepeatRMerVerifier? l_r_mer検証器 = null;
            if (p_引数.A_Is反復rMer検証)
            {
                var l_r長 = p_k長 + C_rMer長のk超過分の既定値;

                var l_窓数 = (p_リード長 ?? 0) - l_r長 + 1;
                if (l_窓数 >= C_rMer検証に必要な窓数)
                {
                    l_r_mer検証器 = RepeatRMerVerifier.V_構築(Get_全リードパス(p_原入力 ?? p_引数), l_r長, l_kmerインデックス, p_k長);
                }
                else
                {
                    Logger.V_出力(メッセージID.rMer検証の見送り, p_k長, l_r長);
                }
            }

            var l_引き継ぎ経路群 = p_引き継ぎ?.Where(x => x.A_Is確定経路).Select(x => x.A_配列).ToList();
            using (new StageTimer($"contig-join k={p_k長}"))
            {
                l_contig構築.V_結合_Contig(l_contigパス, p_引数.A_ペア結合閾値, p_引数.A_ペア支持数閾値, l_コピー数推定.A_コピー数, l_バブル敗者, p_リード長, l_r_mer検証器, p_引数.A_IsGFA出力 ? l_GFAパス : null, l_引き継ぎ経路群, l_コピー数推定.A_コピー数区間);
            }

            Logger.V_出力(メッセージID.Contig構築完了);
            AssemblyStatsReporter.V_出力_統計("contigs", l_contigパス);

            Logger.V_出力_タイムスタンプ();

            var l_IsScaffold作成済み = false;

            if (p_引数.Hasペア)
            {
                Logger.V_出力(メッセージID.Scaffolding開始);
                using (new StageTimer($"scaffolding k={p_k長}"))
                {
                    var l_scaffold構築 = new Scaffolder(l_contig構築, l_contigパス, p_リード長);
                    l_scaffold構築.V_実行(l_scaffoldパス);
                }

                l_IsScaffold作成済み = File.Exists(l_scaffoldパス);
            }

            if (!l_IsScaffold作成済み)
            {
                var l_contig検査 = AssemblyValidator.Get_検査結果(l_contigパス, l_kmerインデックス, p_k長, l_コピー数推定.A_単一コピー基準値);
                AssemblyValidator.V_出力_検査結果("contigs", l_contig検査);
                Logger.V_出力_タイムスタンプ();

                V_用意_次段引き継ぎ(p_次への引き継ぎ, l_contigパス, l_kmerインデックス, p_k長, p_引数, l_バブル敗者, p_合成リードの控え, l_contig構築.A_同一unitig標本, l_contig構築.A_分岐の継ぎ目, l_断片パス);
                var l_contigのみの結果 = new アセンブリ実行結果(p_k長, l_unitigパス, l_contigパス, null, p_引数.A_kmerカットオフ, l_コピー数推定.A_単一コピー基準値, p_引数.A_IsGFA出力 ? l_GFAパス : null, l_contig検査, l_コピー数推定.A_基準の出所);
                V_保存_チェックポイント(l_作業ディレクトリ, p_k長);
                return l_contigのみの結果;
            }

            AssemblyStatsReporter.V_出力_統計("scaffolds", l_scaffoldパス);

            Logger.V_出力(メッセージID.ギャップ充填開始);
            var l_ギャップ統計 = GapFiller.V_充填_ギャップ(l_scaffoldパス, l_kmerインデックス, p_k長);
            GapFiller.V_出力_充填統計(l_ギャップ統計);
            if (l_ギャップ統計.A_埋めたギャップ数 > 0)
            {
                AssemblyStatsReporter.V_出力_統計("scaffolds (gaps filled)", l_scaffoldパス);
            }

            if (p_引数.A_Is局所アセンブリ)
            {
                var l_局所統計 = LocalAssembler.V_充填_ギャップ(l_scaffoldパス, p_引数.A_ライブラリ群, p_k長);
                LocalAssembler.V_出力_統計(l_局所統計);
                if (l_局所統計.A_埋めたギャップ数 > 0)
                {
                    AssemblyStatsReporter.V_出力_統計("scaffolds (local assembly)", l_scaffoldパス);
                }
            }

            var l_scaffoldの検査 = AssemblyValidator.Get_検査結果(l_scaffoldパス, l_kmerインデックス, p_k長, l_コピー数推定.A_単一コピー基準値);
            AssemblyValidator.V_出力_検査結果("scaffolds", l_scaffoldの検査);

            Logger.V_出力_タイムスタンプ();

            V_用意_次段引き継ぎ(p_次への引き継ぎ, l_scaffoldパス, l_kmerインデックス, p_k長, p_引数, l_バブル敗者, p_合成リードの控え, l_contig構築.A_同一unitig標本, l_contig構築.A_分岐の継ぎ目, l_断片パス, p_原入力, p_リード長);
            var l_結果 = new アセンブリ実行結果(p_k長, l_unitigパス, l_contigパス, l_scaffoldパス, p_引数.A_kmerカットオフ, l_コピー数推定.A_単一コピー基準値, p_引数.A_IsGFA出力 ? l_GFAパス : null, l_scaffoldの検査, l_コピー数推定.A_基準の出所);
            V_保存_チェックポイント(l_作業ディレクトリ, p_k長);
            return l_結果;
        }

        /// <summary>
        /// 採用した結果を作業ディレクトリの直下へ複製する
        /// </summary>
        /// <param name="p_結果"></param>
        /// <param name="p_出力ディレクトリ"></param>
        /// <returns></returns>
        public static string V_複製_最終成果物(アセンブリ実行結果 p_結果, string p_出力ディレクトリ)
        {
            V_複製(p_結果.A_unitigパス, Path.Combine(p_出力ディレクトリ, C_Unitigファイル名));
            V_複製(p_結果.A_contigパス, Path.Combine(p_出力ディレクトリ, C_Contigファイル名));
            if (p_結果.A_scaffoldパス is { } l_scaffoldパス)
            {
                V_複製(l_scaffoldパス, Path.Combine(p_出力ディレクトリ, Consts.Scaffoldファイル名));
            }

            if (p_結果.A_GFAパス is { } l_GFAパス)
            {
                V_複製(l_GFAパス, Path.Combine(p_出力ディレクトリ, Consts.GFAファイル名));
            }

            var l_最終パス = Path.Combine(p_出力ディレクトリ, C_最終アセンブリファイル名);
            V_複製(p_結果.A_最終パス, l_最終パス);
            return l_最終パス;
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// この k で決めきれなかった箇所を残す
        /// </summary>
        /// <param name="p_作業ディレクトリ"></param>
        /// <param name="p_k長"></param>
        private static void V_保存_チェックポイント(string p_作業ディレクトリ, int p_k長)
        {
            AmbiguityRecorder.V_保存(p_作業ディレクトリ, p_k長);
        }

        /// <summary>
        /// ファイルを複製する
        /// </summary>
        /// <param name="p_元">複製元</param>
        /// <param name="p_先">複製先</param>
        private static void V_複製(string p_元, string p_先)
        {
            if (p_元 != p_先 && File.Exists(p_元))
            {
                File.Copy(p_元, p_先, overwrite: true);
            }
        }

        /// <summary>
        /// 次の k へ渡す配列とカバレッジを用意する
        /// </summary>
        /// <param name="p_次への引き継ぎ"></param>
        /// <param name="p_FASTAパス"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_引数"></param>
        /// <param name="p_バブル敗者"></param>
        /// <param name="p_合成リードの控え"></param>
        /// <param name="p_断片長標本">合成リードの橋渡し長の上限を見積もる断片長</param>
        /// <param name="p_分岐の継ぎ目">この k の contig が分岐のある継ぎ目で通った辺 ((k+1)-mer、両向き)</param>
        /// <param name="p_原入力">r-mer の裏付けを数える元のリードを持つ設定</param>
        /// <param name="p_リード長">代表リード長</param>
        private static void V_用意_次段引き継ぎ(List<引き継ぎ配列>? p_次への引き継ぎ, string p_FASTAパス, TrustedKmerIndex p_kmerインデックス, int p_k長, Parameters p_引数, IReadOnlyList<string> p_バブル敗者, List<引き継ぎ配列>? p_合成リードの控え, List<int> p_断片長標本, HashSet<string> p_分岐の継ぎ目, string p_断片パス, Parameters? p_原入力 = null, int? p_リード長 = null)
        {
            if (p_次への引き継ぎ is null)
            {
                return;
            }

            Logger.V_出力(メッセージID.引き継ぎの準備開始);

            RepeatRMerVerifier? l_持ち越し検証器 = null;
            if ((p_リード長 ?? 0) - KmerCarryOver.C_持ち越し検証のr長 + 1 >= C_rMer検証に必要な窓数)
            {
                l_持ち越し検証器 = RepeatRMerVerifier.V_構築(
                    Get_全リードパス(p_原入力 ?? p_引数),
                    KmerCarryOver.C_持ち越し検証のr長,
                    p_問い合わせ配列: Get_配列列(p_FASTAパス));
            }

            p_次への引き継ぎ.Clear();
            p_次への引き継ぎ.AddRange(KmerCarryOver.Get_引き継ぎ配列(p_FASTAパス, p_kmerインデックス, p_k長, p_分岐の継ぎ目, l_持ち越し検証器));
            Logger.V_出力_そのまま(FormattableString.Invariant($"[Carry-over] {p_分岐の継ぎ目.Count / 2:N0} branch junction edge(s) are not carried into longer k-mers"));

            foreach (var l_配列 in p_バブル敗者)
            {
                if (l_配列.Length < p_k長)
                {
                    continue;
                }

                p_次への引き継ぎ.Add(Get_引き継ぎ配列(l_配列, p_kmerインデックス, p_k長));
            }

            Logger.V_出力(メッセージID.引き継ぎの準備完了, p_次への引き継ぎ.Count);
            Logger.V_出力_タイムスタンプ();

            if (!p_引数.A_IsSuperRead作成 || !p_引数.Hasペア)
            {
                return;
            }

            if (p_合成リードの控え is { Count: > 0 })
            {
                Logger.V_出力(メッセージID.合成リードを再利用, p_合成リードの控え.Count);
                p_次への引き継ぎ.AddRange(p_合成リードの控え);
                return;
            }

            var l_整列した断片長 = p_断片長標本.Count > 0 ? p_断片長標本.Order().ToArray() : null;
            int? l_断片長上限 = l_整列した断片長 is null ? null : StatsUtil.Get_分位点(l_整列した断片長, C_橋渡しに使う断片長の分位);
            int? l_断片長下限 = l_整列した断片長 is null ? null : StatsUtil.Get_分位点(l_整列した断片長, C_重なりに使う断片長の分位);
            for (var i = 0; i < p_引数.A_ライブラリ数; i++)
            {
                var (A_リード1, A_リード2) = p_引数.A_ライブラリ群[i];
                if (string.IsNullOrWhiteSpace(A_リード2))
                {
                    continue;
                }

                var l_合成リード = SuperReadJoiner.Get_合成リード(A_リード1, A_リード2, p_kmerインデックス, p_k長, out var l_統計, l_断片長上限, l_断片長下限, Get_断片パス(p_断片パス, i));
                SuperReadJoiner.V_出力_統計(l_統計);
                p_次への引き継ぎ.AddRange(l_合成リード);
                p_合成リードの控え?.AddRange(l_合成リード);
            }

            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// 配列を、位置ごとのカバレッジ付きの引き継ぎ配列にする
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
        /// <summary>
        /// FASTA の配列を順に返す
        /// </summary>
        /// <param name="p_FASTAパス">読む FASTA</param>
        /// <returns></returns>
        private static IEnumerable<string> Get_配列列(string p_FASTAパス)
        {
            using FastaReader l_読み込み = new(p_FASTAパス);
            while (l_読み込み.Has続き())
            {
                yield return l_読み込み.Get_次の配列().A_配列;
            }
        }

        /// <summary>
        /// 全ライブラリのリードのパス
        /// </summary>
        /// <param name="p_引数">実行時引数</param>
        /// <returns></returns>
        private static string[] Get_全リードパス(Parameters p_引数)
        {
            return [.. p_引数.A_ライブラリ群.SelectMany(x => new[] { x.A_リード1, x.A_リード2 })];
        }

        /// <summary>
        /// ライブラリごとの断片ファイルのパス
        /// </summary>
        /// <param name="p_基準パス">先頭ライブラリのパス</param>
        /// <param name="p_ライブラリ番号">0 起点のライブラリ番号</param>
        /// <returns></returns>
        private static string Get_断片パス(string p_基準パス, int p_ライブラリ番号)
        {
            if (p_ライブラリ番号 == 0)
            {
                return p_基準パス;
            }

            var l_ディレクトリ = Path.GetDirectoryName(p_基準パス) ?? string.Empty;
            var l_幹 = Path.GetFileNameWithoutExtension(p_基準パス);
            return Path.Combine(l_ディレクトリ, FormattableString.Invariant($"{l_幹}.lib{p_ライブラリ番号 + 1}{Path.GetExtension(p_基準パス)}"));
        }

        /// <summary>
        /// ライブラリの数だけ断片ファイルのパスを並べる
        /// </summary>
        /// <param name="p_基準パス">先頭ライブラリのパス</param>
        /// <param name="p_ライブラリ数">ライブラリの数</param>
        /// <returns></returns>
        private static IEnumerable<string> Get_断片パス群(string p_基準パス, int p_ライブラリ数)
        {
            for (var i = 0; i < p_ライブラリ数; i++)
            {
                yield return Get_断片パス(p_基準パス, i);
            }
        }

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

        /// <summary>
        /// unitig を構築して FASTA へ書き出し、ID -> 配列 の対応を返す
        /// </summary>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_開始kmer"></param>
        /// <param name="p_k長"></param>
        /// <param name="p_出力パス"></param>
        /// <param name="p_Is上限到達"></param>
        /// <returns></returns>
        private static Dictionary<int, string> Get_Unitig(TrustedKmerIndex p_kmerインデックス, List<byte[]> p_開始kmer, int p_k長, string p_出力パス, out bool p_Is上限到達)
        {
            var l_walk結果 = UnitigMaker.Get_walk結果(p_kmerインデックス, p_開始kmer);

            var l_閉路の開始kmer = CyclicUnitigFinder.Get_閉路開始kmer(p_kmerインデックス, l_walk結果, p_k長);
            if (l_閉路の開始kmer.Count > 0)
            {
                Logger.V_出力(メッセージID.分岐のない閉路, l_閉路の開始kmer.Count);
                l_walk結果 =
                    [.. l_walk結果, .. UnitigMaker.Get_walk結果(p_kmerインデックス, l_閉路の開始kmer)];
            }

            HashSet<string> l_既出 = [];
            Dictionary<int, string> l_unitig配列 = [];
            var l_ID = 1;

            using (var l_書き込み = new FastaWriter(p_出力パス))
            {
                foreach (var l_配列 in l_walk結果)
                {
                    if (l_既出.Contains(l_配列) || l_既出.Contains(Util.V_逆相補(l_配列)))
                    {
                        continue;
                    }

                    _ = l_既出.Add(l_配列);
                    _ = l_既出.Add(Util.V_逆相補(l_配列));
                    l_unitig配列[l_ID] = l_配列;
                    l_書き込み.V_書き込み(l_ID++, l_配列);

                    if (l_ID > C_Unitig数の上限)
                    {
                        break;
                    }
                }
            }

            p_Is上限到達 = l_ID > C_Unitig数の上限;
            return l_unitig配列;
        }

        #endregion
    }
}
