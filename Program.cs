using System.Runtime.InteropServices;
using Tsumiki.Common;
using Tsumiki.Core;
using Tsumiki.IO;
using Tsumiki.Utility;

namespace Tsumiki
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            try
            {
                V_実行(args);
            }
            catch (Exception ex)
            {
                Logger.V_出力_エラー("Unhandled Tsumiki's method", ex);
            }
            finally
            {
                var l_一時ディレクトリ = Path.Combine(
                    Environment.CurrentDirectory, ConfigurationManager.A_実行時引数.A_一時ディレクトリ);
                if (Directory.Exists(l_一時ディレクトリ))
                {
                    Directory.Delete(l_一時ディレクトリ, true);
                }
            }
        }

        private static void V_実行(string[] p_引数列)
        {
            if (p_引数列.Length == 0)
            {
                Console.WriteLine(Consts.概要テキスト);
                Environment.Exit(0);
            }

            var l_引数 = ArgumentsReader.Get_実行時引数(p_引数列);
            ConfigurationManager.A_実行時引数 = l_引数;

            if (l_引数.A_ヘルプモードか)
            {
                Console.WriteLine(Consts.ヘルプテキスト);
                Environment.Exit(0);
            }

            Console.WriteLine(l_引数);

            PhredSniffer.V_解決_Phredオフセット(l_引数, l_引数.A_リード1のパス, l_引数.A_リード2のパス);

            Logger.V_出力_タイムスタンプ();

            var l_一時ディレクトリ = Path.Combine(Environment.CurrentDirectory, l_引数.A_一時ディレクトリ);

            if (Path.Exists(l_一時ディレクトリ))
            {
                Console.WriteLine($"{l_引数.A_一時ディレクトリ} already exists");
                Console.WriteLine("Please check path!");
                Environment.Exit(0);
            }

            _ = Directory.CreateDirectory(l_一時ディレクトリ);

            if (l_引数.A_エラー訂正するか)
            {
                Console.WriteLine("Correcting reads before assembly");

                var l_訂正済み1 = Path.Combine(l_一時ディレクトリ, "corrected.1.fq");
                var l_リード2があるか = !string.IsNullOrWhiteSpace(l_引数.A_リード2のパス);
                var l_訂正済み2 = l_リード2があるか ? Path.Combine(l_一時ディレクトリ, "corrected.2.fq") : null;

                ErrorCorrector.V_訂正_リードファイル(
                    l_引数.A_リード1のパス,
                    l_リード2があるか ? l_引数.A_リード2のパス : null,
                    l_一時ディレクトリ,
                    l_訂正済み1,
                    l_訂正済み2);

                // 以降の全処理(k-merカウント・グラフ構築・リードの再マッピング)は
                // 訂正済みファイルを見るようにする。
                l_引数.A_リード1のパス = l_訂正済み1;
                if (l_リード2があるか)
                {
                    l_引数.A_リード2のパス = l_訂正済み2!;
                }

                Logger.V_出力_タイムスタンプ();
            }

            Console.WriteLine("Start construction k-mer index");
            var l_kmerインデックス = new TrustedKmerIndex(l_一時ディレクトリ);
            ConfigurationManager.A_kmerインデックス = l_kmerインデックス;

            if (string.IsNullOrWhiteSpace(l_引数.A_リード2のパス))
            {
                Console.WriteLine("Loading File");
            }
            else
            {
                Console.WriteLine("Loading File1");
            }

            if (l_引数.A_曖昧塩基を許容するか)
            {
                V_読込_リードファイル_曖昧塩基あり(l_引数.A_リード1のパス, l_kmerインデックス);
            }
            else
            {
                KmerCounting.V_読込_リードファイル(l_引数.A_リード1のパス, l_kmerインデックス);
            }

            if (!string.IsNullOrWhiteSpace(l_引数.A_リード2のパス))
            {
                Console.WriteLine("Loading File2");
                if (l_引数.A_曖昧塩基を許容するか)
                {
                    V_読込_リードファイル_曖昧塩基あり(l_引数.A_リード2のパス, l_kmerインデックス);
                }
                else
                {
                    KmerCounting.V_読込_リードファイル(l_引数.A_リード2のパス, l_kmerインデックス);
                }
            }

            Logger.V_出力_タイムスタンプ();

            Console.WriteLine("Applying k-mer cutoff");
            _ = l_kmerインデックス.V_カットオフ(l_引数.A_kmerカットオフ);

            Logger.V_出力_タイムスタンプ();

            Console.WriteLine("Clipping short tips");
            var l_開始kmer = GraphSimplifier.V_除去_tip(l_kmerインデックス, l_引数.A_k長);

            Logger.V_出力_タイムスタンプ();

            Console.WriteLine("Make unitigs");

            var l_ユニティグ構築 = new UnitigMaker(l_kmerインデックス);
            HashSet<string> l_既出 = [];
            Dictionary<int, string> l_ユニティグ配列 = [];
            var l_ID = 1;
            using (var l_書き込み = new FastaWriter(Consts.ユニティグファイル名))
            {
                foreach (var l_kmer in l_開始kmer)
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

            AssemblyStatsReporter.V_出力_統計("unitigs", Consts.ユニティグファイル名);

            // 各 unitig のカバレッジからコピー数を推定する。反復配列かどうかを
            // グラフの形ではなく量的な根拠で判定でき、後段の経路探索では
            // 「この unitig は何回まで使ってよいか」という予算になる。
            // k-mer インデックスがまだ生きているこの時点でしか計算できない。
            var l_ユニティグ長 = l_ユニティグ配列.ToDictionary(x => x.Key, x => x.Value.Length);
            var l_カバレッジ = CopyNumberEstimator.Get_カバレッジ(l_kmerインデックス, l_ユニティグ配列, l_引数.A_k長);
            var l_コピー数推定 = CopyNumberEstimator.Get_推定結果(l_カバレッジ, l_ユニティグ長);
            CopyNumberEstimator.V_出力_推定結果(l_コピー数推定, l_ユニティグ長);

            Logger.V_出力_タイムスタンプ();

            if (l_ID > Consts.ユニティグ数の上限)
            {
                Console.WriteLine("""

                    This genome is too complex to assembly...
                    Please adjust the parameters!

                    """);
                Environment.Exit(0);
            }

            Console.WriteLine("Map reads to unitigs");

            var l_コンティグ構築 = new ContigMaker(Consts.ユニティグファイル名);

            if (string.IsNullOrWhiteSpace(l_引数.A_リード2のパス))
            {
                Console.WriteLine(l_引数.A_リード1のパス);
                l_コンティグ構築.V_マッピング_リード(l_引数.A_リード1のパス);
            }
            else
            {
                // ペアエンドの場合、read1/read2 を同時に読み進めて
                // インサートサイズによる隣接検出も行う。
                Console.WriteLine(l_引数.A_リード1のパス);
                Console.WriteLine(l_引数.A_リード2のパス);
                l_コンティグ構築.V_マッピング_ペアリード(l_引数.A_リード1のパス, l_引数.A_リード2のパス);
            }

            Logger.V_出力_タイムスタンプ();

            Console.WriteLine("unite unitigs");

            l_コンティグ構築.V_結合_コンティグ(
                Consts.コンティグファイル名, l_引数.A_ペア結合閾値, l_引数.A_ペア支持数閾値, l_コピー数推定.A_コピー数);

            Console.WriteLine("Maked contigs");

            AssemblyStatsReporter.V_出力_統計("contigs", Consts.コンティグファイル名);

            Logger.V_出力_タイムスタンプ();

            // スキャフォールディングはペアエンド情報を前提とするため、
            // read2 が指定されている(=ペアエンドで実行された)場合のみ行う。
            if (!string.IsNullOrWhiteSpace(l_引数.A_リード2のパス))
            {
                Console.WriteLine("Scaffolding contigs");

                var l_スキャフォールド構築 = new Scaffolder(l_コンティグ構築, Consts.コンティグファイル名);
                l_スキャフォールド構築.V_実行(Consts.スキャフォールドファイル名);

                AssemblyStatsReporter.V_出力_統計("scaffolds", Consts.スキャフォールドファイル名);

                // スキャフォールドの N を、グラフ上で両端を繋ぐ経路を探して
                // 実配列に置き換える。contig が途切れたのは配列が無いからではなく
                // 分岐で決められなかったからであることが多く、その場合
                // ギャップを埋める配列はグラフ上に実在する。
                Console.WriteLine("Filling scaffold gaps");
                var l_ギャップ統計 = GapFiller.V_充填_ギャップ(
                    Consts.スキャフォールドファイル名, l_kmerインデックス, l_引数.A_k長);
                GapFiller.V_出力_充填統計(l_ギャップ統計);
                if (l_ギャップ統計.A_埋めたギャップ数 > 0)
                {
                    AssemblyStatsReporter.V_出力_統計("scaffolds (gaps filled)", Consts.スキャフォールドファイル名);
                }

                // 出来上がったアセンブリが、観測された k-mer とその出現回数に
                // 対して辻褄が合っているかを自己検査する(リファレンス不要)。
                AssemblyValidator.V_出力_検査結果(
                    "scaffolds",
                    AssemblyValidator.Get_検査結果(
                        Consts.スキャフォールドファイル名, l_kmerインデックス, l_引数.A_k長, l_コピー数推定.A_単一コピー基準値));

                Logger.V_出力_タイムスタンプ();
            }
            else
            {
                AssemblyValidator.V_出力_検査結果(
                    "contigs",
                    AssemblyValidator.Get_検査結果(
                        Consts.コンティグファイル名, l_kmerインデックス, l_引数.A_k長, l_コピー数推定.A_単一コピー基準値));

                Logger.V_出力_タイムスタンプ();
            }

            Console.WriteLine("開発中！");

            Logger.V_出力_タイムスタンプ();
        }

        /// <summary>
        /// 曖昧塩基を許容する経路。現状シングルスレッドのまま(ワーカー番号固定)で、
        /// 呼ばれる頻度が低い想定のため並列化の優先度を下げている。
        ///
        /// KmerCounting 側と同様、以前はここでもリードを逆相補にしてもう一度
        /// すべての k-mer を登録していた。登録側が正規化するようになった今は
        /// 完全に冗長で、すべてのカウントをちょうど2倍にしてしまうため削除した。
        /// </summary>
        private static void V_読込_リードファイル_曖昧塩基あり(string p_ファイルパス, TrustedKmerIndex p_kmerインデックス)
        {
            ulong l_件数 = 0;
            ulong l_ログ回数 = 0;
            var l_k長 = ConfigurationManager.A_実行時引数.A_k長;
            var l_Phredオフセット = ConfigurationManager.A_実行時引数.A_Phredオフセット;
            var l_クオリティカットオフ = ConfigurationManager.A_実行時引数.A_クオリティカットオフ;

            using var l_読み込み = new FastqReader(p_ファイルパス);
            while (l_読み込み.Get_続きがあるか())
            {
                var l_リード = l_読み込み.Get_次のリード();
                if (l_リード.A_塩基候補列!.Count < l_k長)
                {
                    continue;
                }
                var l_塩基候補 = CollectionsMarshal.AsSpan(l_リード.A_塩基候補列);
                for (var i = 0; i < l_リード.A_クオリティ.Length; i++)
                {
                    if (l_リード.A_クオリティ[i] - l_Phredオフセット - l_クオリティカットオフ < 0)
                    {
                        l_塩基候補[i] = [Consts.塩基ID.A, Consts.塩基ID.C, Consts.塩基ID.G, Consts.塩基ID.T];
                    }
                }
                p_kmerインデックス.V_登録_曖昧塩基あり(l_塩基候補[..l_k長], 0);
                for (var i = l_k長; i < l_リード.A_塩基候補列.Count; i++)
                {
                    p_kmerインデックス.V_登録_曖昧塩基あり(l_塩基候補.Slice(i - l_k長 + 1, l_k長), 0);
                }
                if (++l_件数 == Consts.進捗ログ間隔)
                {
                    Console.WriteLine((++l_ログ回数 * Consts.進捗ログ間隔) + " reads Loaded");
                    l_件数 = 0;
                }
            }
            Console.WriteLine($"Loaded {(l_ログ回数 * Consts.進捗ログ間隔) + l_件数} reads from {Path.GetFileName(p_ファイルパス)}");
        }
    }
}
