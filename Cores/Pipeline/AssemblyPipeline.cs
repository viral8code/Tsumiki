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
        private const int rMer長のk超過分の既定値 = 10;

        /// <summary>
        /// r-mer 検証に必要な窓数
        /// </summary>
        private const int rMer検証に必要な窓数 = 20;

        /// <summary>
        /// unitig ファイル名
        /// </summary>
        private const string Unitigファイル名 = "unitigs.fasta";

        /// <summary>
        /// contig ファイル名
        /// </summary>
        private const string Contigファイル名 = "contigs.fasta";

        /// <summary>
        /// 最終アセンブリファイル名
        /// </summary>
        private const string 最終アセンブリファイル名 = "assembly.fasta";

        /// <summary>
        /// unitig 数の上限
        /// </summary>
        private const int Unitig数の上限 = 100_000;

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
        /// <remarks>
        /// unitig 数が上限を超えた場合は null<br/>
        /// 生成物は k ごとの作業ディレクトリに置く<br/>
        /// 最終的に採用したものだけを V_複製_最終成果物 が作業ディレクトリの直下へ複製する
        /// </remarks>
        public static アセンブリ実行結果? Get_実行結果(Parameters p_引数, int p_k長, string p_一時ディレクトリ, int? p_リード長, IReadOnlyList<引き継ぎ配列>? p_引き継ぎ = null, List<引き継ぎ配列>? p_次への引き継ぎ = null, List<引き継ぎ配列>? p_合成リードの控え = null, Parameters? p_原入力 = null)
        {
            // 以降の全処理は ConfigurationManager 経由で k 長を参照する
            // 明示指定の印は立てない (自動選択された値のままとして扱う)
            if (p_引数.A_k長 != p_k長)
            {
                p_引数.Set_推定k長(p_k長);
            }

            var l_作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"k{p_k長}");
            _ = Directory.CreateDirectory(l_作業ディレクトリ);

            AmbiguityRecorder.V_開始(p_k長);

            var l_unitigパス = Path.Combine(l_作業ディレクトリ, Unitigファイル名);
            var l_contigパス = Path.Combine(l_作業ディレクトリ, Contigファイル名);
            var l_scaffoldパス = Path.Combine(l_作業ディレクトリ, Consts.Scaffoldファイル名);
            var l_GFAパス = Path.Combine(l_作業ディレクトリ, Consts.GFAファイル名);

            // 再実行で今回生成されなかった前回の scaffold を採用しない
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
            KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_kmerインデックス);
            _ = l_kmerインデックス.V_カットオフ(p_引数.A_kmerカットオフ);

            KmerHistogram.V_出力_スペクトル(l_kmerインデックス.A_出現回数ヒストグラム, p_k長, p_リード長);

            // 救済はカットオフの直後に行う
            // スペクトルを出し終えてからにするのは、
            // 救済した k-mer を混ぜたヒストグラムがカバレッジ推定の材料に
            // 使われないようにするため
            if (p_引数.A_Is救済kmer使用)
            {
                Logger.V_出力(メッセージID.救済kmerの開始);
                _ = MercyKmerRescuer.Get_救済数(p_引数, l_kmerインデックス, p_k長);
            }

            // 引き継ぎはカットオフの後に行う
            // カウントとカットオフを通常どおり
            // 済ませてから足すことで、スペクトルが実データのまま保たれ、
            // ゲノムサイズとカバレッジの推定が歪まない
            if (p_引き継ぎ is { Count: > 0 })
            {
                Logger.V_出力(メッセージID.引き継ぎの統合開始, p_引き継ぎ.Count);
                var l_追加数 = KmerCarryOver.V_引き継ぎ(p_引き継ぎ, l_kmerインデックス, p_k長, p_リード長);
                Logger.V_出力(メッセージID.引き継ぎで追加したkmer数, l_追加数);
            }

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.tip除去開始);

            // tip 除去は k-mer 集合を縮小するため、開始点はその後の状態で
            // 数え直す必要がある
            // 除去側が最終状態のものを返す
            var l_開始kmer = GraphSimplifier.V_除去_tip(l_kmerインデックス, p_k長, p_リード長);

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.Unitig構築開始);
            var l_unitig配列 = Get_Unitig(l_kmerインデックス, l_開始kmer, p_k長, l_unitigパス, out var l_上限に達したか);

            AssemblyStatsReporter.V_出力_統計("unitigs", l_unitigパス);

            if (l_上限に達したか)
            {
                Logger.V_出力(メッセージID.グラフが複雑すぎる, p_k長, Unitig数の上限);
                return null;
            }

            Logger.V_出力(メッセージID.リードのマッピング開始);
            var l_contig構築 = new ContigMaker(l_unitigパス);

            // 反復配列かどうかをグラフの形ではなく量的な根拠で判定するための
            // コピー数推定
            // k-mer インデックスが生きている今しか計算できない
            // 接続構造 (排他的な鎖) による補正のため、ContigMaker が厳密な
            // de Bruijn グラフから構築した隣接情報も使う (コピー数推定専用に
            // 作る使い捨てのグラフで、V_結合_contig が後で作るものとは別)
            var l_unitig長 = l_unitig配列.ToDictionary(x => x.Key, x => x.Value.Length);
            var l_カバレッジ = CopyNumberEstimator.Get_カバレッジ(l_kmerインデックス, l_unitig配列, p_k長);
            var l_グラフ = l_contig構築.Get_グラフ();
            var l_コピー数推定 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_unitig長, l_グラフ);
            CopyNumberEstimator.V_出力_推定結果(l_コピー数推定, l_unitig長);

            Logger.V_出力_タイムスタンプ();

            if (string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                Logger.V_出力(メッセージID.リードファイルのパス, p_引数.A_リード1のパス);
                l_contig構築.V_マッピング_リード(p_引数.A_リード1のパス);
            }
            else
            {
                // ペアエンドの場合、read1/read2 を同時に読み進めて
                // インサートサイズによる隣接検出も行う
                Logger.V_出力(メッセージID.リードファイルのパス, p_引数.A_リード1のパス);
                Logger.V_出力(メッセージID.リードファイルのパス, p_引数.A_リード2のパス);
                l_contig構築.V_マッピング_ペアリード(p_引数.A_リード1のパス, p_引数.A_リード2のパス);
            }

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.Unitig結合開始);

            // careful_bubble: バブル除去で外れた側の配列も、この k では
            // 敗者と判断しただけであって存在しないわけではない
            // 捨てずに
            // 次の k への引き継ぎ候補として集めておく
            List<string> l_バブル敗者 = [];

            // 短い反復解決の拒否権 (-rv)
            // 生リードを 1 回走査し、集計された
            // ペア支持だけでは決めきれない「本当にその接合点を跨いだリードが
            // あるか」を確かめる
            // r はこの k の k-1 重なりより確実に長く
            // 取らないと共有区間の内側に収まってしまい判定にならないため、
            // k + 余剰分で決め、リード内に十分な窓を取れる場合に検証する
            RepeatRMerVerifier? l_r_mer検証器 = null;
            if (p_引数.A_Is反復rMer検証)
            {
                var l_r長 = p_k長 + rMer長のk超過分の既定値;

                // r は k より長くないと拒否権として働かない (窓が k-1 の共有区間に
                // 収まってしまい常に真になる)
                // 一方でリードから取れる窓の数は
                // リード長 - r + 1 なので、r がリード長に近づくと r-mer の
                // カバレッジが痩せて正しい経路まで棄却しはじめる
                var l_窓数 = (p_リード長 ?? 0) - l_r長 + 1;
                if (l_窓数 >= rMer検証に必要な窓数)
                {
                    l_r_mer検証器 = RepeatRMerVerifier.V_構築([(p_原入力 ?? p_引数).A_リード1のパス, (p_原入力 ?? p_引数).A_リード2のパス], l_r長);
                }
                else
                {
                    Logger.V_出力(メッセージID.rMer検証の見送り, p_k長, l_r長);
                }
            }

            l_contig構築.V_結合_Contig(l_contigパス, p_引数.A_ペア結合閾値, p_引数.A_ペア支持数閾値, l_コピー数推定.A_コピー数, l_バブル敗者, p_リード長, l_r_mer検証器, p_引数.A_IsGFA出力 ? l_GFAパス : null);
            Logger.V_出力(メッセージID.Contig構築完了);
            AssemblyStatsReporter.V_出力_統計("contigs", l_contigパス);

            Logger.V_出力_タイムスタンプ();

            // scaffolding はペアエンド情報を前提とする
            // インサートサイズが推定できず作られないこともある
            var l_IsScaffold作成済み = false;

            if (!string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                Logger.V_出力(メッセージID.Scaffolding開始);
                var l_scaffold構築 = new Scaffolder(l_contig構築, l_contigパス, p_リード長);
                l_scaffold構築.V_実行(l_scaffoldパス);
                l_IsScaffold作成済み = File.Exists(l_scaffoldパス);
            }

            if (!l_IsScaffold作成済み)
            {
                var l_contig検査 = AssemblyValidator.Get_検査結果(l_contigパス, l_kmerインデックス, p_k長, l_コピー数推定.A_単一コピー基準値);
                AssemblyValidator.V_出力_検査結果("contigs", l_contig検査);
                Logger.V_出力_タイムスタンプ();

                V_用意_次段引き継ぎ(p_次への引き継ぎ, l_contigパス, l_kmerインデックス, p_k長, p_引数, l_バブル敗者, p_合成リードの控え);
                var l_contigのみの結果 = new アセンブリ実行結果(p_k長, l_unitigパス, l_contigパス, null, p_引数.A_kmerカットオフ, l_コピー数推定.A_単一コピー基準値, p_引数.A_IsGFA出力 ? l_GFAパス : null, l_contig検査);
                V_保存_チェックポイント(l_作業ディレクトリ, p_k長);
                return l_contigのみの結果;
            }

            AssemblyStatsReporter.V_出力_統計("scaffolds", l_scaffoldパス);

            // contig が途切れる原因は配列の不在より分岐の未解決が多く、
            // その場合ギャップを埋める配列はグラフ上に実在する
            Logger.V_出力(メッセージID.ギャップ充填開始);
            var l_ギャップ統計 = GapFiller.V_充填_ギャップ(l_scaffoldパス, l_kmerインデックス, p_k長);
            GapFiller.V_出力_充填統計(l_ギャップ統計);
            if (l_ギャップ統計.A_埋めたギャップ数 > 0)
            {
                AssemblyStatsReporter.V_出力_統計("scaffolds (gaps filled)", l_scaffoldパス);
            }

            // GapFiller が埋めきれなかった残りを、その両端に実際にマップされた
            // 局所リードだけの使い捨てミニアセンブリで埋める (-la、-mg の安全な代替)
            if (p_引数.A_Is局所アセンブリ)
            {
                var l_局所統計 = LocalAssembler.V_充填_ギャップ(l_scaffoldパス, p_引数.A_リード1のパス, p_引数.A_リード2のパス, p_k長);
                LocalAssembler.V_出力_統計(l_局所統計);
                if (l_局所統計.A_埋めたギャップ数 > 0)
                {
                    AssemblyStatsReporter.V_出力_統計("scaffolds (local assembly)", l_scaffoldパス);
                }
            }

            var l_scaffoldの検査 = AssemblyValidator.Get_検査結果(l_scaffoldパス, l_kmerインデックス, p_k長, l_コピー数推定.A_単一コピー基準値);
            AssemblyValidator.V_出力_検査結果("scaffolds", l_scaffoldの検査);

            Logger.V_出力_タイムスタンプ();

            V_用意_次段引き継ぎ(p_次への引き継ぎ, l_scaffoldパス, l_kmerインデックス, p_k長, p_引数, l_バブル敗者, p_合成リードの控え);
            var l_結果 = new アセンブリ実行結果(p_k長, l_unitigパス, l_contigパス, l_scaffoldパス, p_引数.A_kmerカットオフ, l_コピー数推定.A_単一コピー基準値, p_引数.A_IsGFA出力 ? l_GFAパス : null, l_scaffoldの検査);
            V_保存_チェックポイント(l_作業ディレクトリ, p_k長);
            return l_結果;
        }

        /// <summary>
        /// 採用した結果を作業ディレクトリの直下へ複製する
        /// </summary>
        /// <param name="p_結果"></param>
        /// <param name="p_出力ディレクトリ"></param>
        /// <returns></returns>
        /// <remarks>
        /// k ごとの成果物は k のサブディレクトリに残したまま、利用者が受け取る 1 組だけを上へ出す<br/>
        /// unitigs/contigs/scaffolds はここまでの各段階の出力で、最後に手が入る前の姿<br/>
        /// 利用者が使うべき 1 本は assembly.fasta のほうになる
        /// </remarks>
        public static string V_複製_最終成果物(アセンブリ実行結果 p_結果, string p_出力ディレクトリ)
        {
            V_複製(p_結果.A_unitigパス, Path.Combine(p_出力ディレクトリ, Unitigファイル名));
            V_複製(p_結果.A_contigパス, Path.Combine(p_出力ディレクトリ, Contigファイル名));
            if (p_結果.A_scaffoldパス is { } l_scaffoldパス)
            {
                V_複製(l_scaffoldパス, Path.Combine(p_出力ディレクトリ, Consts.Scaffoldファイル名));
            }
            if (p_結果.A_GFAパス is { } l_GFAパス)
            {
                V_複製(l_GFAパス, Path.Combine(p_出力ディレクトリ, Consts.GFAファイル名));
            }

            var l_最終パス = Path.Combine(p_出力ディレクトリ, 最終アセンブリファイル名);
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
        private static void V_用意_次段引き継ぎ(List<引き継ぎ配列>? p_次への引き継ぎ, string p_FASTAパス, TrustedKmerIndex p_kmerインデックス, int p_k長, Parameters p_引数, IReadOnlyList<string> p_バブル敗者, List<引き継ぎ配列>? p_合成リードの控え)
        {
            if (p_次への引き継ぎ is null)
            {
                return;
            }

            // ここから次の k が始まるまでは出力の無い区間が続く
            // 何をしている
            // ところなのかが分からないと止まったように見えるため、区切りを出す
            Logger.V_出力(メッセージID.引き継ぎの準備開始);

            p_次への引き継ぎ.Clear();
            p_次への引き継ぎ.AddRange(KmerCarryOver.Get_引き継ぎ配列(p_FASTAパス, p_kmerインデックス, p_k長));

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

            if (!p_引数.A_IsSuperRead作成 || string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                return;
            }

            // 合成リードは最初の k で作ったものを以降の k でも使い回す
            // 橋渡しはその k の信頼できる k-mer 集合を通るので k ごとに
            // 作り直していたが、実データでは本数がほとんど動かなかった
            // (7.4 Mbp ・ 170 x で 556,352 -> 557,761 -> 557,867 -> 558,026 -> 557,931)
            // 一方で費用は k とともに増え、6 つの k の合計で実行時間の
            // 3 分の 1 を占めていた
            // 狙いは次の k のためにリードを実効的に
            // 伸ばすことなので、最も繋がりやすい最小の k で 1 度作れば足りる
            if (p_合成リードの控え is { Count: > 0 })
            {
                Logger.V_出力(メッセージID.合成リードを再利用, p_合成リードの控え.Count);
                p_次への引き継ぎ.AddRange(p_合成リードの控え);
                return;
            }

            var l_合成リード = SuperReadJoiner.Get_合成リード(p_引数.A_リード1のパス, p_引数.A_リード2のパス, p_kmerインデックス, p_k長, out var l_統計);
            SuperReadJoiner.V_出力_統計(l_統計);
            p_次への引き継ぎ.AddRange(l_合成リード);
            p_合成リードの控え?.AddRange(l_合成リード);
            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// 配列を、位置ごとのカバレッジ付きの引き継ぎ配列にする
        /// </summary>
        /// <param name="p_配列"></param>
        /// <param name="p_kmerインデックス"></param>
        /// <param name="p_k長"></param>
        /// <returns></returns>
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
        /// <remarks>
        /// 同じ配列を順鎖・逆鎖の両方で出さないよう既出集合で弾く
        /// </remarks>
        private static Dictionary<int, string> Get_Unitig(TrustedKmerIndex p_kmerインデックス, List<byte[]> p_開始kmer, int p_k長, string p_出力パス, out bool p_Is上限到達)
        {
            var l_walk結果 = UnitigMaker.Get_walk結果(p_kmerインデックス, p_開始kmer);

            // 分岐を 1 つも持たない閉路は開始点の条件を満たす k-mer を持たず、
            // ここまでの走査から丸ごと漏れる
            // 覆い残しを拾って足す
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

                    if (l_ID > Unitig数の上限)
                    {
                        break;
                    }
                }
            }

            p_Is上限到達 = l_ID > Unitig数の上限;
            return l_unitig配列;
        }

        #endregion
    }
}
