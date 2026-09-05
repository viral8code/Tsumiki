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
        /// p_出力接頭辞 は出力ファイル名の先頭に付く(単一 k では空文字)。
        /// </summary>
        public static アセンブリ実行結果? Get_実行結果(
            Parameters p_引数, int p_k長, string p_一時ディレクトリ, string p_出力接頭辞, int? p_リード長)
        {
            // 以降の全処理は ConfigurationManager 経由で k 長を参照する。
            // 明示指定の印は立てない(自動選択された値のままとして扱う)。
            if (p_引数.A_k長 != p_k長)
            {
                p_引数.Set_推定k長(p_k長);
            }

            var l_作業ディレクトリ = Path.Combine(p_一時ディレクトリ, $"k{p_k長}");
            _ = Directory.CreateDirectory(l_作業ディレクトリ);

            var l_ユニティグパス = p_出力接頭辞 + Consts.ユニティグファイル名;
            var l_コンティグパス = p_出力接頭辞 + Consts.コンティグファイル名;
            var l_スキャフォールドパス = p_出力接頭辞 + Consts.スキャフォールドファイル名;

            Console.WriteLine("Start construction k-mer index");
            using var l_kmerインデックス = new TrustedKmerIndex(l_作業ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_kmerインデックス;

            V_読込_リード(p_引数, l_kmerインデックス);

            Logger.V_出力_タイムスタンプ();

            Console.WriteLine("Applying k-mer cutoff");
            KmerCutoffSelector.V_解決_kmerカットオフ(p_引数, l_kmerインデックス);
            _ = l_kmerインデックス.V_カットオフ(p_引数.A_kmerカットオフ);

            KmerHistogram.V_出力_スペクトル(l_kmerインデックス.A_出現回数ヒストグラム, p_k長, p_リード長);

            Logger.V_出力_タイムスタンプ();

            Console.WriteLine("Clipping short tips");
            // tip 除去は k-mer 集合を縮小するため、開始点はその後の状態で
            // 数え直す必要がある。除去側が最終状態のものを返す。
            var l_開始kmer = GraphSimplifier.V_除去_tip(l_kmerインデックス, p_k長);

            Logger.V_出力_タイムスタンプ();

            Console.WriteLine("Make unitigs");
            var l_ユニティグ配列 = Get_ユニティグ(
                l_kmerインデックス, l_開始kmer, l_ユニティグパス, out var l_上限に達したか);

            AssemblyStatsReporter.V_出力_統計("unitigs", l_ユニティグパス);

            // 反復配列かどうかをグラフの形ではなく量的な根拠で判定するための
            // コピー数推定。k-mer インデックスが生きている今しか計算できない。
            var l_ユニティグ長 = l_ユニティグ配列.ToDictionary(x => x.Key, x => x.Value.Length);
            var l_カバレッジ = CopyNumberEstimator.Get_カバレッジ(l_kmerインデックス, l_ユニティグ配列, p_k長);
            var l_コピー数推定 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_ユニティグ長);
            CopyNumberEstimator.V_出力_推定結果(l_コピー数推定, l_ユニティグ長);

            Logger.V_出力_タイムスタンプ();

            if (l_上限に達したか)
            {
                Console.WriteLine($"[Warning] The graph is too complex to assemble at k={p_k長} " +
                    $"(unitig count exceeded {Consts.ユニティグ数の上限}). Skipping this k.");
                return null;
            }

            Console.WriteLine("Map reads to unitigs");
            var l_コンティグ構築 = new ContigMaker(l_ユニティグパス);
            if (string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                Console.WriteLine(p_引数.A_リード1のパス);
                l_コンティグ構築.V_マッピング_リード(p_引数.A_リード1のパス);
            }
            else
            {
                // ペアエンドの場合、read1/read2 を同時に読み進めて
                // インサートサイズによる隣接検出も行う。
                Console.WriteLine(p_引数.A_リード1のパス);
                Console.WriteLine(p_引数.A_リード2のパス);
                l_コンティグ構築.V_マッピング_ペアリード(p_引数.A_リード1のパス, p_引数.A_リード2のパス);
            }

            Logger.V_出力_タイムスタンプ();

            Console.WriteLine("unite unitigs");
            l_コンティグ構築.V_結合_コンティグ(
                l_コンティグパス, p_引数.A_ペア結合閾値, p_引数.A_ペア支持数閾値, l_コピー数推定.A_コピー数);
            Console.WriteLine("Maked contigs");
            AssemblyStatsReporter.V_出力_統計("contigs", l_コンティグパス);

            Logger.V_出力_タイムスタンプ();

            // スキャフォールディングはペアエンド情報を前提とする。
            // インサートサイズが推定できず作られないこともある。
            var l_スキャフォールドを作ったか = false;
            if (!string.IsNullOrWhiteSpace(p_引数.A_リード2のパス))
            {
                Console.WriteLine("Scaffolding contigs");
                var l_スキャフォールド構築 = new Scaffolder(l_コンティグ構築, l_コンティグパス);
                l_スキャフォールド構築.V_実行(l_スキャフォールドパス);
                l_スキャフォールドを作ったか = File.Exists(l_スキャフォールドパス);
            }

            if (!l_スキャフォールドを作ったか)
            {
                AssemblyValidator.V_出力_検査結果(
                    "contigs",
                    AssemblyValidator.Get_検査結果(
                        l_コンティグパス, l_kmerインデックス, p_k長, l_コピー数推定.A_単一コピー基準値));
                Logger.V_出力_タイムスタンプ();

                return new アセンブリ実行結果(
                    p_k長, l_ユニティグパス, l_コンティグパス, null,
                    p_引数.A_kmerカットオフ, l_コピー数推定.A_単一コピー基準値);
            }

            AssemblyStatsReporter.V_出力_統計("scaffolds", l_スキャフォールドパス);

            // contig が途切れる原因は配列の不在より分岐の未解決が多く、
            // その場合ギャップを埋める配列はグラフ上に実在する。
            Console.WriteLine("Filling scaffold gaps");
            var l_ギャップ統計 = GapFiller.V_充填_ギャップ(l_スキャフォールドパス, l_kmerインデックス, p_k長);
            GapFiller.V_出力_充填統計(l_ギャップ統計);
            if (l_ギャップ統計.A_埋めたギャップ数 > 0)
            {
                AssemblyStatsReporter.V_出力_統計("scaffolds (gaps filled)", l_スキャフォールドパス);
            }

            AssemblyValidator.V_出力_検査結果(
                "scaffolds",
                AssemblyValidator.Get_検査結果(
                    l_スキャフォールドパス, l_kmerインデックス, p_k長, l_コピー数推定.A_単一コピー基準値));

            Logger.V_出力_タイムスタンプ();

            return new アセンブリ実行結果(
                p_k長, l_ユニティグパス, l_コンティグパス, l_スキャフォールドパス,
                p_引数.A_kmerカットオフ, l_コピー数推定.A_単一コピー基準値);
        }

        private static void V_読込_リード(Parameters p_引数, TrustedKmerIndex p_kmerインデックス)
        {
            var l_ペアエンドか = !string.IsNullOrWhiteSpace(p_引数.A_リード2のパス);
            Console.WriteLine(l_ペアエンドか ? "Loading File1" : "Loading File");

            if (p_引数.A_曖昧塩基を許容するか)
            {
                KmerCounting.V_読込_リードファイル_曖昧塩基あり(p_引数.A_リード1のパス, p_kmerインデックス);
            }
            else
            {
                KmerCounting.V_読込_リードファイル(p_引数.A_リード1のパス, p_kmerインデックス);
            }

            if (!l_ペアエンドか)
            {
                return;
            }

            Console.WriteLine("Loading File2");
            if (p_引数.A_曖昧塩基を許容するか)
            {
                KmerCounting.V_読込_リードファイル_曖昧塩基あり(p_引数.A_リード2のパス, p_kmerインデックス);
            }
            else
            {
                KmerCounting.V_読込_リードファイル(p_引数.A_リード2のパス, p_kmerインデックス);
            }
        }

        /// <summary>
        /// unitig を構築して FASTA へ書き出し、ID -> 配列 の対応を返す。
        /// 同じ配列を順鎖・逆鎖の両方で出さないよう既出集合で弾く。
        /// </summary>
        private static Dictionary<int, string> Get_ユニティグ(
            TrustedKmerIndex p_kmerインデックス, List<byte[]> p_開始kmer, string p_出力パス, out bool p_上限に達したか)
        {
            var l_ユニティグ構築 = new UnitigMaker(p_kmerインデックス);
            HashSet<string> l_既出 = [];
            Dictionary<int, string> l_ユニティグ配列 = [];
            var l_ID = 1;

            using (var l_書き込み = new FastaWriter(p_出力パス))
            {
                foreach (var l_kmer in p_開始kmer)
                {
                    var l_ユニティグ = l_ユニティグ構築.Get_ユニティグ(l_kmer);
                    if (l_既出.Contains(l_ユニティグ.A_配列) || l_既出.Contains(Util.V_逆相補(l_ユニティグ.A_配列)))
                    {
                        continue;
                    }
                    _ = l_既出.Add(l_ユニティグ.A_配列);
                    _ = l_既出.Add(Util.V_逆相補(l_ユニティグ.A_配列));
                    l_ユニティグ配列[l_ID] = l_ユニティグ.A_配列;
                    l_書き込み.V_書き込み(l_ID++, l_ユニティグ.A_配列);
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
