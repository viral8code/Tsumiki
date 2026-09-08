using Tsumiki.Common;
using Tsumiki.IO;
using Tsumiki.Model;
using Tsumiki.Utility;

namespace Tsumiki.Core
{
    /// <summary>
    /// 指定した k 長で、k-mer カウントからスキャフォールドまでを一通り実行する。
    /// </summary>
    internal static class AssemblyPipeline
    {
        /// <summary>
        /// p_k長 でアセンブリを実行し、生成物のパスを返す。
        /// unitig 数が上限を超えた場合は null。
        /// 生成物は k ごとの作業ディレクトリに置く。最終的に採用したものだけを
        /// V_複製_最終成果物 が実行ディレクトリへ複製する。
        /// </summary>
        public static アセンブリ実行結果? Get_実行結果(
            Parameters p_引数, int p_k長, string p_一時ディレクトリ, int? p_リード長,
            IReadOnlyList<引き継ぎ配列>? p_引き継ぎ = null,
            List<引き継ぎ配列>? p_次への引き継ぎ = null)
        {
            // 以降の全処理は ConfigurationManager 経由で k 長を参照する。
            // 明示指定の印は立てない(自動選択された値のままとして扱う)。
            if (p_引数.A_k長 != p_k長)
            {
                p_引数.Set_推定k長(p_k長);
            }

            var l_作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"k{p_k長}");
            _ = Directory.CreateDirectory(l_作業ディレクトリ);

            // 同じ条件で作り終えているなら、この k は丸ごと飛ばす。
            var l_署名 = CheckpointStore.Get_署名(p_引数, p_k長);
            if (p_引数.A_再開するか
                && CheckpointStore.Get_再開結果(l_作業ディレクトリ, l_署名, p_次への引き継ぎ) is { } l_再開結果)
            {
                Logger.V_出力(メッセージID.再開_kを飛ばした, p_k長);
                AmbiguityRecorder.V_読み込み(l_作業ディレクトリ, p_k長);
                return l_再開結果;
            }

            AmbiguityRecorder.V_開始(p_k長);

            var l_ユニティグパス = Path.Combine(l_作業ディレクトリ, Consts.ユニティグファイル名);
            var l_コンティグパス = Path.Combine(l_作業ディレクトリ, Consts.コンティグファイル名);
            var l_スキャフォールドパス = Path.Combine(l_作業ディレクトリ, Consts.スキャフォールドファイル名);
            var l_GFAパス = Path.Combine(l_作業ディレクトリ, Consts.GFAファイル名);

            Logger.V_出力(メッセージID.kmerインデックス構築開始);
            using var l_kmerインデックス = new TrustedKmerIndex(l_作業ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_kmerインデックス;

            KmerCounting.V_読込_リードペア(p_引数, l_kmerインデックス, p_進行状況を出力するか: true);

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.kmerカットオフ適用);
            KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_kmerインデックス);
            _ = l_kmerインデックス.V_カットオフ(p_引数.A_kmerカットオフ);

            KmerHistogram.V_出力_スペクトル(l_kmerインデックス.A_出現回数ヒストグラム, p_k長, p_リード長);

            // 救済はカットオフの直後に行う。スペクトルを出し終えてからにするのは、
            // 救済した k-mer を混ぜたヒストグラムがカバレッジ推定の材料に
            // 使われないようにするため。
            if (p_引数.A_救済kmerを使うか)
            {
                Logger.V_出力(メッセージID.救済kmerの開始);
                _ = MercyKmerRescuer.Get_救済数(p_引数, l_kmerインデックス, p_k長);
                Logger.V_出力_タイムスタンプ();
            }

            // 引き継ぎはカットオフの後に行う。カウントとカットオフを通常どおり
            // 済ませてから足すことで、スペクトルが実データのまま保たれ、
            // ゲノムサイズとカバレッジの推定が歪まない。
            if (p_引き継ぎ is { Count: > 0 })
            {
                Logger.V_出力(メッセージID.引き継ぎの統合開始, p_引き継ぎ.Count);
                var l_追加数 = KmerCarryOver.V_引き継ぎ(p_引き継ぎ, l_kmerインデックス, p_k長, p_リード長);
                Logger.V_出力(メッセージID.引き継ぎで追加したkmer数, l_追加数);
            }

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.tip除去開始);
            // tip 除去は k-mer 集合を縮小するため、開始点はその後の状態で
            // 数え直す必要がある。除去側が最終状態のものを返す。
            var l_開始kmer = GraphSimplifier.V_除去_tip(l_kmerインデックス, p_k長, p_リード長);

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.ユニティグ構築開始);
            var l_ユニティグ配列 = Get_ユニティグ(
                l_kmerインデックス, l_開始kmer, p_k長, l_ユニティグパス, out var l_上限に達したか);

            AssemblyStatsReporter.V_出力_統計("unitigs", l_ユニティグパス);

            if (l_上限に達したか)
            {
                Logger.V_出力(メッセージID.グラフが複雑すぎる, p_k長, Consts.ユニティグ数の上限);
                return null;
            }

            Logger.V_出力(メッセージID.リードのマッピング開始);
            var l_コンティグ構築 = new ContigMaker(l_ユニティグパス);

            // 反復配列かどうかをグラフの形ではなく量的な根拠で判定するための
            // コピー数推定。k-mer インデックスが生きている今しか計算できない。
            // 接続構造(排他的な鎖)による補正のため、ContigMaker が厳密な
            // de Bruijn グラフから構築した隣接情報も使う(コピー数推定専用に
            // 作る使い捨てのグラフで、V_結合_コンティグ が後で作るものとは別)。
            var l_ユニティグ長 = l_ユニティグ配列.ToDictionary(x => x.Key, x => x.Value.Length);
            var l_カバレッジ = CopyNumberEstimator.Get_カバレッジ(l_kmerインデックス, l_ユニティグ配列, p_k長);
            var l_グラフ = l_コンティグ構築.Get_グラフ();
            var l_コピー数推定 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_ユニティグ長, l_グラフ);
            CopyNumberEstimator.V_出力_推定結果(l_コピー数推定, l_ユニティグ長);

            Logger.V_出力_タイムスタンプ();

            if (string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                Logger.V_出力(メッセージID.リードファイルのパス, p_引数.A_リード1のパス);
                l_コンティグ構築.V_マッピング_リード(p_引数.A_リード1のパス);
            }
            else
            {
                // ペアエンドの場合、read1/read2 を同時に読み進めて
                // インサートサイズによる隣接検出も行う。
                Logger.V_出力(メッセージID.リードファイルのパス, p_引数.A_リード1のパス);
                Logger.V_出力(メッセージID.リードファイルのパス, p_引数.A_リード2のパス);
                l_コンティグ構築.V_マッピング_ペアリード(p_引数.A_リード1のパス, p_引数.A_リード2のパス);
            }

            Logger.V_出力_タイムスタンプ();

            Logger.V_出力(メッセージID.ユニティグ結合開始);
            // careful_bubble: バブル除去で外れた側の配列も、この k では
            // 敗者と判断しただけであって存在しないわけではない。捨てずに
            // 次の k への引き継ぎ候補として集めておく。
            List<string> l_バブル敗者 = [];

            // 短い反復解決の拒否権(-rv)。生リードを1回走査し、集計された
            // ペア支持だけでは決めきれない「本当にその接合点を跨いだリードが
            // あるか」を確かめる。r はこの k の k-1 重なりより確実に長く
            // 取らないと共有区間の内側に収まってしまい判定にならないため、
            // k + 余剰分で決める(32塩基を超えるkでは意味を持てないため見送る)。
            RepeatRMerVerifier? l_r_mer検証器 = null;
            if (p_引数.A_反復をrMerで検証するか)
            {
                var l_r長 = p_k長 + Consts.rMer長のk超過分の既定値;
                if (l_r長 <= 32)
                {
                    l_r_mer検証器 = RepeatRMerVerifier.V_構築(
                        [p_引数.A_リード1のパス, p_引数.A_リード2のパス], l_r長);
                }
                else
                {
                    Logger.V_出力(メッセージID.rMer検証の見送り, p_k長, l_r長);
                }
            }

            l_コンティグ構築.V_結合_コンティグ(
                l_コンティグパス, p_引数.A_ペア結合閾値, p_引数.A_ペア支持数閾値, l_コピー数推定.A_コピー数,
                l_バブル敗者, p_リード長, l_r_mer検証器, p_引数.A_GFAを出力するか ? l_GFAパス : null);
            Logger.V_出力(メッセージID.コンティグ構築完了);
            AssemblyStatsReporter.V_出力_統計("contigs", l_コンティグパス);

            Logger.V_出力_タイムスタンプ();

            // スキャフォールディングはペアエンド情報を前提とする。
            // インサートサイズが推定できず作られないこともある。
            var l_スキャフォールドを作ったか = false;
            if (!string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                Logger.V_出力(メッセージID.スキャフォールディング開始);
                var l_スキャフォールド構築 = new Scaffolder(l_コンティグ構築, l_コンティグパス, p_リード長);
                l_スキャフォールド構築.V_実行(l_スキャフォールドパス);
                l_スキャフォールドを作ったか = File.Exists(l_スキャフォールドパス);
            }

            if (!l_スキャフォールドを作ったか)
            {
                var l_コンティグの検査 = AssemblyValidator.Get_検査結果(
                    l_コンティグパス, l_kmerインデックス, p_k長, l_コピー数推定.A_単一コピー基準値);
                AssemblyValidator.V_出力_検査結果("contigs", l_コンティグの検査);
                Logger.V_出力_タイムスタンプ();

                V_用意_次への引き継ぎ(p_次への引き継ぎ, l_コンティグパス, l_kmerインデックス, p_k長, p_引数, l_バブル敗者);
                var l_コンティグのみの結果 = new アセンブリ実行結果(
                    p_k長, l_ユニティグパス, l_コンティグパス, null,
                    p_引数.A_kmerカットオフ, l_コピー数推定.A_単一コピー基準値,
                    p_引数.A_GFAを出力するか ? l_GFAパス : null,
                    l_コンティグの検査);
                V_保存_チェックポイント(
                    l_作業ディレクトリ, l_署名, l_コンティグのみの結果, p_次への引き継ぎ, p_k長);
                return l_コンティグのみの結果;
            }

            AssemblyStatsReporter.V_出力_統計("scaffolds", l_スキャフォールドパス);

            // contig が途切れる原因は配列の不在より分岐の未解決が多く、
            // その場合ギャップを埋める配列はグラフ上に実在する。
            Logger.V_出力(メッセージID.ギャップ充填開始);
            var l_ギャップ統計 = GapFiller.V_充填_ギャップ(l_スキャフォールドパス, l_kmerインデックス, p_k長);
            GapFiller.V_出力_充填統計(l_ギャップ統計);
            if (l_ギャップ統計.A_埋めたギャップ数 > 0)
            {
                AssemblyStatsReporter.V_出力_統計("scaffolds (gaps filled)", l_スキャフォールドパス);
            }

            // GapFiller が埋めきれなかった残りを、その両端に実際にマップされた
            // 局所リードだけの使い捨てミニアセンブリで埋める(-la、-mgの安全な代替)。
            if (p_引数.A_局所アセンブリするか)
            {
                var l_局所統計 = LocalAssembler.V_充填_ギャップ(
                    l_スキャフォールドパス, p_引数.A_リード1のパス, p_引数.A_リード2のパス, p_k長, l_作業ディレクトリ);
                LocalAssembler.V_出力_統計(l_局所統計);
                if (l_局所統計.A_埋めたギャップ数 > 0)
                {
                    AssemblyStatsReporter.V_出力_統計("scaffolds (local assembly)", l_スキャフォールドパス);
                }
            }

            var l_スキャフォールドの検査 = AssemblyValidator.Get_検査結果(
                l_スキャフォールドパス, l_kmerインデックス, p_k長, l_コピー数推定.A_単一コピー基準値);
            AssemblyValidator.V_出力_検査結果("scaffolds", l_スキャフォールドの検査);

            Logger.V_出力_タイムスタンプ();

            V_用意_次への引き継ぎ(p_次への引き継ぎ, l_スキャフォールドパス, l_kmerインデックス, p_k長, p_引数, l_バブル敗者);
            var l_結果 = new アセンブリ実行結果(
                p_k長, l_ユニティグパス, l_コンティグパス, l_スキャフォールドパス,
                p_引数.A_kmerカットオフ, l_コピー数推定.A_単一コピー基準値,
                p_引数.A_GFAを出力するか ? l_GFAパス : null,
                l_スキャフォールドの検査);
            V_保存_チェックポイント(l_作業ディレクトリ, l_署名, l_結果, p_次への引き継ぎ, p_k長);
            return l_結果;
        }

        /// <summary>
        /// この k を作り終えたことを記録する。決めきれなかった箇所の控えも
        /// 一緒に残す(再開でこの k を飛ばしたときに、レポートから消えないように)。
        /// </summary>
        private static void V_保存_チェックポイント(
            string p_作業ディレクトリ, string p_署名, アセンブリ実行結果 p_結果,
            IReadOnlyList<引き継ぎ配列>? p_次への引き継ぎ, int p_k長)
        {
            AmbiguityRecorder.V_保存(p_作業ディレクトリ, p_k長);
            CheckpointStore.V_保存(p_作業ディレクトリ, p_署名, p_結果, p_次への引き継ぎ);
        }

        /// <summary>
        /// 採用した結果を実行ディレクトリへ複製する。k ごとの成果物は
        /// 一時ディレクトリに残したまま、利用者が受け取る1組だけを外へ出す。
        /// </summary>
        public static void V_複製_最終成果物(アセンブリ実行結果 p_結果)
        {
            V_複製(p_結果.A_ユニティグパス, Consts.ユニティグファイル名);
            V_複製(p_結果.A_コンティグパス, Consts.コンティグファイル名);
            if (p_結果.A_スキャフォールドパス is { } l_スキャフォールドパス)
            {
                V_複製(l_スキャフォールドパス, Consts.スキャフォールドファイル名);
            }
            if (p_結果.A_GFAパス is { } l_GFAパス)
            {
                V_複製(l_GFAパス, Consts.GFAファイル名);
            }
        }

        private static void V_複製(string p_元, string p_先)
        {
            if (p_元 != p_先 && File.Exists(p_元))
            {
                File.Copy(p_元, p_先, overwrite: true);
            }
        }

        /// <summary>
        /// 次の k へ渡す配列とカバレッジを用意する。
        /// k-mer インデックスが破棄される前でなければ作れない。
        ///
        /// -sr が有効なら、この k の信頼できる k-mer 集合の中でペアを橋渡しして
        /// 作った合成リード(SuperRead)も足す。元のリードの2〜4倍の長さを持つため、
        /// マルチ k の上限(リード長で頭打ちになる)を実効的に外せる。
        ///
        /// バブル除去で外れた側の配列(careful_bubble)も足す。この k での
        /// 敗者判定は次の k を拘束しない。
        /// </summary>
        private static void V_用意_次への引き継ぎ(
            List<引き継ぎ配列>? p_次への引き継ぎ, string p_FASTAパス,
            TrustedKmerIndex p_kmerインデックス, int p_k長, Parameters p_引数,
            IReadOnlyList<string> p_バブル敗者)
        {
            if (p_次への引き継ぎ is null)
            {
                return;
            }
            // ここから次の k が始まるまでは出力の無い区間が続く。何をしている
            // ところなのかが分からないと止まったように見えるため、区切りを出す。
            Logger.V_出力(メッセージID.引き継ぎの準備開始);

            p_次への引き継ぎ.Clear();
            p_次への引き継ぎ.AddRange(
                KmerCarryOver.Get_引き継ぎ配列(p_FASTAパス, p_kmerインデックス, p_k長));

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

            if (p_引数.A_SuperReadを作るか && !string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                var l_合成リード = SuperReadJoiner.Get_合成リード(
                    p_引数.A_リード1のパス, p_引数.A_リード2のパス, p_kmerインデックス, p_k長, out var l_統計);
                SuperReadJoiner.V_出力_統計(l_統計);
                p_次への引き継ぎ.AddRange(l_合成リード);
                Logger.V_出力_タイムスタンプ();
            }
        }

        /// <summary>配列を、位置ごとのカバレッジ付きの引き継ぎ配列にする。</summary>
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
        /// unitig を構築して FASTA へ書き出し、ID -> 配列 の対応を返す。
        /// 同じ配列を順鎖・逆鎖の両方で出さないよう既出集合で弾く。
        /// </summary>
        private static Dictionary<int, string> Get_ユニティグ(
            TrustedKmerIndex p_kmerインデックス, List<byte[]> p_開始kmer, int p_k長,
            string p_出力パス, out bool p_上限に達したか)
        {
            var l_walk結果 = UnitigMaker.Get_walk結果(p_kmerインデックス, p_開始kmer);

            // 分岐を1つも持たない閉路は開始点の条件を満たす k-mer を持たず、
            // ここまでの走査から丸ごと漏れる。覆い残しを拾って足す。
            var l_閉路の開始kmer = CyclicUnitigFinder.Get_閉路の開始kmer(
                p_kmerインデックス, l_walk結果, p_k長);
            if (l_閉路の開始kmer.Count > 0)
            {
                Logger.V_出力(メッセージID.分岐のない閉路, l_閉路の開始kmer.Count);
                l_walk結果 =
                    [.. l_walk結果, .. UnitigMaker.Get_walk結果(p_kmerインデックス, l_閉路の開始kmer)];
            }

            HashSet<string> l_既出 = [];
            Dictionary<int, string> l_ユニティグ配列 = [];
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
                    l_ユニティグ配列[l_ID] = l_配列;
                    l_書き込み.V_書き込み(l_ID++, l_配列);
                    if (l_ID > Consts.ユニティグ数の上限)
                    {
                        break;
                    }
                }
            }

            p_上限に達したか = l_ID > Consts.ユニティグ数の上限;
            return l_ユニティグ配列;
        }
    }
}
