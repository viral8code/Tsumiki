using Tsumiki.Model;

namespace Tsumiki.Common
{
    /// <summary>
    /// 言語ごとの文言表。<see cref="Messages"/> から引く。
    ///
    /// 言語を増やすときは、その言語の辞書をここに足して <see cref="Get_辞書"/> に
    /// 繋ぐだけでよい。訳が無い ID は英語にそのまま落ちるので、部分的に訳した
    /// 状態でも表示が壊れない。
    /// </summary>
    internal static class MessageCatalog
    {
        /// <summary>p_ID の書式文字列。指定言語に無ければ英語を返す。</summary>
        public static string Get_書式(言語 p_言語, メッセージID p_ID)
        {
            var l_辞書 = Get_辞書(p_言語);
            if (l_辞書 is not null && l_辞書.TryGetValue(p_ID, out var l_書式))
            {
                return l_書式;
            }
            return _英語.TryGetValue(p_ID, out var l_英語) ? l_英語 : p_ID.ToString();
        }

        /// <summary>その言語にこの ID の訳があるか。訳の入れ忘れの検査に使う。</summary>
        public static bool Get_訳があるか(言語 p_言語, メッセージID p_ID)
        {
            return Get_辞書(p_言語)?.ContainsKey(p_ID) == true;
        }

        /// <summary>その言語の辞書。まだ用意していない言語は null。</summary>
        private static IReadOnlyDictionary<メッセージID, string>? Get_辞書(言語 p_言語)
        {
            return p_言語 switch
            {
                言語.英語 => _英語,
                言語.日本語 => _日本語,
                言語.中国語 => _中国語,
                _ => null,
            };
        }

        private static readonly Dictionary<メッセージID, string> _英語 = new()
        {
            [メッセージID.デブルーイングラフの要約] =
                "[Debug] Exact de Bruijn unitig graph: {0} directed edge(s), {1} branching vertex(es) out of {2}.",
            [メッセージID.先読みで解決した分岐数] =
                "[Debug] Beam-search lookahead resolved {0} further junction(s) that the single-step mutual-uniqueness rule could not decide.",
            [メッセージID.GFA出力完了] =
                "[Info] Wrote unitig graph to {0} (GFA1).",
            [メッセージID.分岐選択の重み内訳] =
                "[Debug] Branch-selection weights: {0} single-read adjacency pair(s) + {1} paired-end pair(s).",
            [メッセージID.単純化の収束] =
                "[Debug] Simplification converged after {0} round(s).",
            [メッセージID.単純化の打ち切り] =
                "[Debug] Simplification stopped at the round limit ({0}) without fully converging.",
            [メッセージID.バブル除去数] =
                "[Debug] Popped {0} simple bubble branch(es) total (kept as standalone contigs; only their graph edges were removed).",
            [メッセージID.反復解決数] =
                "[Debug] Repeat resolution: {0} short repeat(s) (<= {1}bp) total were duplicated and untangled using read pairs that span them.",
            [メッセージID.辺選択の内訳] =
                "[Debug] Edge selection: {0} vertex(es) had a single out-edge, {1} branch(es) resolved by read support, {2} branch(es) left unresolved because they leave a multi-copy repeat (reads inside a repeat cannot tell the copies apart).",
            [メッセージID.反復通り抜けで棄却した結合数] =
                "[Debug] {0} join(s) were refused because they would chain through a multi-copy repeat that has not been untangled (doing so skips whatever lies between the repeat's copies).",
            [メッセージID.相互一意で残った結合数] =
                "[Debug] {0} directed merge(s) survived the mutual-uniqueness check ({1} undirected join(s)).",
            [メッセージID.環状コンティグあり] =
                "[Info] {0} contig(s) closed into a circle ({1}). A closed circle means that replicon (chromosome or plasmid) was assembled end to end.",
            [メッセージID.環状コンティグなし] =
                "[Info] No contig closed into a circle; every replicon is still fragmented.",
            [メッセージID.確定辺標本数] =
                "[Info] InsertSize samples derived from resolved (actually-joined) unitig adjacency: {0}.",
            [メッセージID.確定辺標本の中央値] =
                "[Info] Resolved-edge sample median: {0} (from {1} samples; not subject to the same-unitig length bias).",
            [メッセージID.短すぎるユニティグの除外] =
                "[Warning] {0} unitig(s) shorter than k-mer length were skipped in mapping.",
            [メッセージID.曖昧なkmer登録] =
                "[Warning] {0} k-mer registration(s) were ambiguous (shared by multiple unitigs) and will be ignored during mapping.",
            [メッセージID.ペア隣接候補数] =
                "[Info] Paired-end adjacency candidates detected: {0} edges ({1} supporting pairs total).",
            [メッセージID.同一ユニティグのペア向き集計] =
                "[Info] Same-unitig pair orientation counts: same-orientation={0}, opposite-orientation={1}. Using '{2}' as the library's observed orientation for InsertSize estimation ({3} samples).",
            [メッセージID.同一ユニティグの断片長分布] =
                "[Info] Same-unitig fragment-length distribution: {0}.",
            [メッセージID.同一ユニティグの断片長中央値] =
                "[Info] Same-unitig fragment-length median: {0} (from {1} samples; read lengths added back to the inner distance, so this is a true fragment length. May still be biased short if unitigs are shorter than the true insert size).",
            [メッセージID.ペアリードIDの不一致] =
                "[Warning] Paired read IDs do not match at this position (\"{0}\" vs \"{1}\"). Paired-end adjacency detection may be unreliable for reads after this point; single-read adjacency detection is unaffected.",
            [メッセージID.統合できる接合点なし] =
                "[Merge] No junction in the backbone was spanned by another k; keeping it as-is.",
            [メッセージID.統合した接合点数] =
                "[Merge] Joined {0} backbone junction(s) using sequence from other k values.",
            [メッセージID.候補一覧の見出し] =
                "[Multi-k] Candidate assemblies (reference-free evaluation against a common anchor k-mer set):",
            [メッセージID.候補一覧の明細] =
                "[Multi-k]   k={0,3}: {1}{2}",
            [メッセージID.統計_ファイルなし] =
                "[Stats] {0}: (file not found: {1})",
            [メッセージID.統計] =
                "[Stats] {0}: {1}",
            [メッセージID.統計_長さで絞り込み] =
                "[Stats] {0} (>={1}bp, comparable to abyss-fac): {2}",
            [メッセージID.検査_対象外のk長] =
                "[Check] {0}: skipped (the self-check supports k <= 64 only).",
            [メッセージID.検査_取りこぼし] =
                "[Check] {0}: {1:N0} trusted k-mer(s); {2:N0} ({3:0.00}%) do not appear in the assembly at all (sequence that was trimmed away or never reached by any path).",
            [メッセージID.検査_出しすぎ] =
                "[Check] {0}: {1:N0} k-mer instance(s) in the assembly; {2:N0} ({3:0.00}%) are more copies than the coverage supports across {4:N0} distinct k-mer(s) -- this is the part of the total length that is inflated.",
            [メッセージID.kmerインデックス構築開始] =
                "Start construction k-mer index",
            [メッセージID.kmerカットオフ適用] =
                "Applying k-mer cutoff",
            [メッセージID.引き継ぎで追加したkmer数] =
                "[Carry-over] Added {0:N0} k-mer(s) from the previous k that this k did not observe.",
            [メッセージID.tip除去開始] =
                "Clipping short tips",
            [メッセージID.ユニティグ構築開始] =
                "Make unitigs",
            [メッセージID.グラフが複雑すぎる] =
                "[Warning] The graph is too complex to assemble at k={0} (unitig count exceeded {1}). Skipping this k.",
            [メッセージID.リードのマッピング開始] =
                "Map reads to unitigs",
            [メッセージID.ユニティグ結合開始] =
                "unite unitigs",
            [メッセージID.rMer検証の見送り] =
                "[Info] Repeat r-mer verification skipped for k={0}: r-mer length would be {1}bp, exceeding the 32bp packing limit.",
            [メッセージID.コンティグ構築完了] =
                "Maked contigs",
            [メッセージID.スキャフォールディング開始] =
                "Scaffolding contigs",
            [メッセージID.ギャップ充填開始] =
                "Filling scaffold gaps",
            [メッセージID.kが薄すぎて省略] =
                "[Multi-k] Skipping k={0}: predicted k-mer coverage {1:F1}x is below {2:F1}x, so this k cannot produce a usable graph.",
            [メッセージID.kの開始見出し] =
                "[Multi-k] ===== Assembling with k={0} =====",
            [メッセージID.kでアセンブリできず] =
                "[Multi-k] k={0} produced no assembly; skipping it.",
            [メッセージID.単一のkのみ成功] =
                "[Multi-k] Only one k produced an assembly; using it without comparison.",
            [メッセージID.アンカーkmer集合の構築] =
                "[Multi-k] Building the common anchor k-mer set (k={0}) for comparison",
            [メッセージID.アンカースペクトルが二峰でない] =
                "[Multi-k] The anchor k-mer spectrum is not bimodal; candidates cannot be compared.",
            [メッセージID.候補を評価できない] =
                "[Multi-k] Could not evaluate the candidates; falling back to the largest k.",
            [メッセージID.採用したk] =
                "[Multi-k] Selected k={0}.",
            [メッセージID.統合開始] =
                "[Merge] Merging the other k values into the selected assembly",
            [メッセージID.統合結果を評価できない] =
                "[Merge] Could not evaluate the merged assembly; keeping the selected one.",
            [メッセージID.統合前の評価] =
                "[Merge]   before: {0}",
            [メッセージID.統合後の評価] =
                "[Merge]   merged: {0}",
            [メッセージID.統合が骨格に勝てず] =
                "[Merge] The merged assembly did not beat the selected one; keeping the selected one.",
            [メッセージID.統合結果を採用] =
                "[Merge] Using the merged assembly.",
            [メッセージID.エラー訂正_スペクトル構築] =
                "[ErrorCorrection] Building k-mer spectrum...",
            [メッセージID.エラー訂正_訂正開始] =
                "[ErrorCorrection] Correcting reads...",
            [メッセージID.エラー訂正_ファイル別統計] =
                "[ErrorCorrection] {0}: {1}/{2} reads corrected ({3} base corrections total).",
            [メッセージID.リード読込完了] =
                "Loaded {0} reads from {1}",
            [メッセージID.リード2の読込開始] =
                "Loading File2",
            [メッセージID.前処理統計] =
                "[Preprocess] {0}/{1} pair(s) showed adapter read-through (trimmed to the overlapping fragment length); {2} base(s) mutually corrected.",
            [メッセージID.SuperRead統計] =
                "[SuperRead] {0}/{1} pair(s) bridged into a single synthetic read and carried into the next k.",
            [メッセージID.ギャップ充填_対象なし] =
                "[Info] Gap filling: no gaps to fill.",
            [メッセージID.ギャップ充填統計] =
                "[Info] Gap filling: {0}/{1} gap(s) closed with real sequence ({2:N0}bp of N replaced). {3} left as N because more than one path fits, {4} because no path through the graph connects the two sides.",
            [メッセージID.局所アセンブリ_対象なし] =
                "[Info] Local assembly: no remaining gap was eligible.",
            [メッセージID.局所アセンブリ統計] =
                "[Info] Local assembly: {0}/{1} remaining gap(s) closed with locally re-assembled reads ({2:N0}bp of N replaced). {3} had no local read at all, {4} had more than one local path, {5} had no local path connecting the two sides.",
            [メッセージID.スキャフォールディング省略_インサートサイズ不明] =
                "[Info] Scaffolding skipped: insert size was not specified and could not be estimated from mapped pairs.",
            [メッセージID.スキャフォールディング開始_インサートサイズ] =
                "[Info] Scaffolding with insert size = {0}",
            [メッセージID.スキャフォールディング省略_コンティグなし] =
                "[Info] Scaffolding skipped: no contigs were found.",
            [メッセージID.内部を指したペア候補] =
                "[Info] {0} pair-end candidate(s) pointed at unitigs interior to an already-joined contig and were skipped (endpoint already resolved by contig construction).",
            [メッセージID.未配置を指したペア候補] =
                "[Info] {0} pair-end candidate(s) referenced unitigs that were not placed into any contig (e.g. too short) and were skipped.",
            [メッセージID.閾値後のスキャフォールド辺] =
                "[Info] Scaffold edges resolved after thresholding: {0}; {1} rejected by the mutual-uniqueness check, {2} kept.",
            [メッセージID.スキャフォールド出力完了] =
                "[Info] Wrote {0} scaffold(s), total length {1}, to {2}",
            [メッセージID.インサートサイズ推定_同一ユニティグ] =
                "[Info] Insert size auto-estimated as {0} from {1} same-unitig sampled pairs (median; unitig N50 {2} is >= {3}x the estimate, so the short-fragment truncation bias does not apply).",
            [メッセージID.インサートサイズ推定_確定辺] =
                "[Info] Insert size auto-estimated as {0} from {1} resolved-edge sampled pairs (median, preferred over same-unitig samples because the unitigs are not long enough for same-unitig samples to be unbiased).",
            [メッセージID.インサートサイズ推定_標本不足] =
                "[Info] Insert size auto-estimation requires at least {0} samples; only {1} resolved-edge and {2} total samples were collected.",
            [メッセージID.インサートサイズ推定_全標本] =
                "[Info] Insert size auto-estimated as {0} from {1} sampled pairs (median; resolved-edge samples were too few ({2}), fell back to the full pool which may be biased short).",
            [メッセージID.単一コピー基準値] =
                "[Info] Single-copy coverage baseline estimated as {0:0.#} (length-weighted median).",
            [メッセージID.コピー数の要約] =
                "[Info] Estimated copy numbers -- {0}",
            [メッセージID.反復配列の割合] =
                "[Info] Multi-copy (repeat) content: {0:N0}bp of {1:N0}bp ({2:0.0}% of the assembly is sequence that occurs more than once).",
            [メッセージID.グラフ単純化の反復] =
                "[GraphSimplifier] Iteration {0}: examined {1} unitig(s) (tip threshold < {2}bp, coverage baseline {3:0.#}), removed {4} tip(s), trimmed {5} low-coverage k-mer(s) from {6} unitig edge(s).",
            [メッセージID.rMer検証による棄却] =
                "[Debug] r-mer verification vetoed {0} repeat duplication(s) whose winning pairing was not confirmed by reads actually crossing the junction (pair-count evidence alone would have duplicated them).",
            [メッセージID.例外を無視_見出し] =
                "[Warning] The following exception will be ignored.",
            [メッセージID.例外を無視_メソッド] =
                "          Method: {0}",
            [メッセージID.停止_見出し] =
                "[Error] Program was stopped.",
            [メッセージID.停止_メソッド] =
                "        Method: {0}",
            [メッセージID.タイムスタンプ] =
                "[Log] {0:yyyy/MM/dd HH:mm:ss}",
            [メッセージID.Phred_明示指定と不一致] =
                "[Warning] Quality strings look like Phred{0}, but -p {1} was given explicitly. Honouring the explicit value; re-run with -p {2} if the data really is Phred{3}.",
            [メッセージID.Phred_自動判定] =
                "[Info] Phred offset auto-detected as {0} from the quality strings (observed ASCII range [{1}, {2}]). Pass -p explicitly to override.",
            [メッセージID.Phred_検査の警告] =
                "[Warning] Phred encoding check for {0}: {1}",
            [メッセージID.kmerカットオフ_混合モデル] =
                "[Info] k-mer cutoff auto-selected as {0} from a 2-component mixture model (error + genomic; single-copy mean {1:0.#}, trusted from {2}, {3} EM iteration(s)). Pass -kc explicitly to override.",
            [メッセージID.kmerカットオフ_谷が不明] =
                "[Info] Could not identify a clear histogram valley (the spectrum may not be bimodal at this coverage); keeping -kc {0}.",
            [メッセージID.kmerカットオフ_スペクトル] =
                "[Info] k-mer cutoff auto-selected as {0} from the k-mer spectrum. Pass -kc explicitly to override.",
            [メッセージID.kmerヒストグラム] =
                "[Info] k-mer count histogram (count:#distinct kmers): {0}",
            [メッセージID.スペクトルの谷が不明] =
                "[Info] Could not identify a clear histogram valley (the spectrum may not be bimodal at this coverage).",
            [メッセージID.スペクトルの谷とピーク] =
                "[Info] k-mer spectrum: valley at count {0} ({1} distinct kmers), single-copy peak at count {2} ({3} distinct kmers)",
            [メッセージID.ゲノムサイズとカバレッジの推定] =
                "[Info] Estimated genome size: {0:N0} bp, estimated coverage: {1}",
            [メッセージID.k自動選択_リード長不明] =
                "[Info] Could not sample a read length; keeping -k {0}.",
            [メッセージID.k自動選択_kが長すぎる] =
                "[Warning] -k {0} is not shorter than the observed read length ({1} bp); no k-mer can be extracted from a read. Lower -k.",
            [メッセージID.k自動選択_リード長が短い] =
                "[Info] Observed read length ({0} bp) is too short to pick a k automatically; keeping -k {1}.",
            [メッセージID.k自動選択] =
                "[Info] k auto-selected as {0} from the observed read length ({1} bp). Pass -k explicitly to override.",
            [メッセージID.開始kmerの探索] =
                "Search First k-mer",
            [メッセージID.一時ディレクトリを残した] =
                "[Info] Per-k assemblies are kept in {0} (use {1} to delete it).",
            [メッセージID.リード長の観測値] =
                "[Info] Read length (median of sampled reads): {0} bp",
            [メッセージID.一時ディレクトリが既にある] =
                "{0} already exists",
            [メッセージID.パスの確認] =
                "Please check path!",
            [メッセージID.前処理省略_ペアなし] =
                "[Warning] -pp requires paired-end reads (read2 is not set). Skipping preprocessing.",
            [メッセージID.前処理開始] =
                "Preprocessing reads (adapter trim + pair correction)",
            [メッセージID.エラー訂正開始] =
                "Correcting reads before assembly",
            [メッセージID.開発中] =
                "開発中！",
            [メッセージID.コンティグ総延長] =
                "Total Length of contigs : {0}",
            [メッセージID.リードファイルのパス] =
                "{0}",
            [メッセージID.試すk一覧] =
                "[Multi-k] Trying k = {0}",
            [メッセージID.リード読込の進捗] =
                "{0} reads Loaded",
            [メッセージID.リード1の読込開始] =
                "Loading File1",
            [メッセージID.単一リードの読込開始] =
                "Loading File",
            [メッセージID.スキャフォールド候補辺数] =
                "[Info] Scaffold candidate edges (contig-level, before thresholding): {0}{1}",
            [メッセージID.理想本数モデルあり] =
                "; ideal-count model available",
            [メッセージID.理想本数モデルなし] =
                "; ideal-count model unavailable",
            [メッセージID.Phred_ファイル間で不一致] =
                "[Warning] Phred offset inference disagreed between the two read files (read1 -> {0}, read2 -> {1}). Keeping -p {2} as-is.",
            [メッセージID.Phred_未確定] =
                "undetermined",
            [メッセージID.kmer種類数] =
                "kmer count: {0}",
            [メッセージID.採用kmer数] =
                "good kmer: {0}",
            [メッセージID.アセンブリ不能] =
                "\nThis genome is too complex to assembly...\nPlease adjust the parameters!\n",
            [メッセージID.概要_説明] =
                "Tsumiki is a genome assembler.",
            [メッセージID.概要_作者] =
                "author: {0}",
            [メッセージID.概要_バージョン] =
                "version: {0}",
            [メッセージID.ヘルプ_使い方] =
                "Usage: tsumiki {0} <path> [{1} <path>] [options]",
            [メッセージID.ヘルプ節_入力] =
                "# Input",
            [メッセージID.ヘルプ節_kmerと品質] =
                "# k-mer / quality",
            [メッセージID.ヘルプ節_ペアエンド] =
                "# Paired-end / scaffolding",
            [メッセージID.ヘルプ節_前処理] =
                "# Preprocessing",
            [メッセージID.ヘルプ節_マルチk] =
                "# Multi-k assembly",
            [メッセージID.ヘルプ節_反復の安全策] =
                "# Repeat resolution safeguards",
            [メッセージID.ヘルプ節_出力とその他] =
                "# Output / misc",
            [メッセージID.ヘルプ_リード1] =
                "forward fastq(.gz) path (required; also use this for single-end reads)",
            [メッセージID.ヘルプ_リード2] =
                "reverse fastq(.gz) path",
            [メッセージID.ヘルプ_曖昧塩基] =
                "allow ambiguous bases (e.g. N) in reads (default: false)",
            [メッセージID.ヘルプ_k長] =
                "k-mer length; a comma-separated list (e.g. 31,63,95) tries each and keeps the best, like {0} (default: auto from read length, capped at {1})",
            [メッセージID.ヘルプ_kmerカットオフ] =
                "minimum k-mer count to trust (default: auto from the k-mer spectrum)",
            [メッセージID.ヘルプ_Phred] =
                "phred score base ({0}) (default: {1})",
            [メッセージID.ヘルプ_クオリティカットオフ] =
                "minimum base quality to trust (default: {0})",
            [メッセージID.ヘルプ_メモリ予算] =
                "memory budget for k-mer counting, e.g. 2G, 512M (default: {0})",
            [メッセージID.ヘルプ_インサートサイズ] =
                "expected insert size (default: auto-estimated from mapped pairs)",
            [メッセージID.ヘルプ_ペア結合閾値] =
                "minimum dominance ratio to accept a scaffold edge (default: {0})",
            [メッセージID.ヘルプ_ペア支持数閾値] =
                "minimum read-pair support to resolve a short repeat (default: {0})",
            [メッセージID.ヘルプ_積極性モード] =
                "preset for {0}/{1}: trades completeness for safety (default: {2})",
            [メッセージID.ヘルプ_前処理] =
                "trim adapter read-through and cross-correct low-quality bases via mate overlap (default: false)",
            [メッセージID.ヘルプ_エラー訂正] =
                "correct reads from the k-mer spectrum before assembly (default: false)",
            [メッセージID.ヘルプ_マルチk] =
                "assemble at several k values and keep the best (default: false; up to {0}x runtime)",
            [メッセージID.ヘルプ_引き継ぎなし] =
                "don't carry sequence from one k to the next (default: carry)",
            [メッセージID.ヘルプ_SuperRead] =
                "with paired-end reads, bridge pairs into synthetic long reads and carry those too (default: false)",
            [メッセージID.ヘルプ_マージ] =
                "splice sequence from other k values into junctions the best k left open (default: false; can raise misassemblies)",
            [メッセージID.ヘルプ_反復rMer検証] =
                "require raw-read r-mer support before duplicating a short repeat (default: false)",
            [メッセージID.ヘルプ_局所アセンブリ] =
                "reassemble unclosed scaffold gaps from reads mapped near their edges; safer than {0} (default: false)",
            [メッセージID.ヘルプ_GFA出力] =
                "also write {0} (the unitig graph, GFA1) for viewers like Bandage (default: false)",
            [メッセージID.ヘルプ_一時ディレクトリ] =
                "temp directory, keeps the per-k assemblies (default: {0})",
            [メッセージID.ヘルプ_一時ディレクトリ削除] =
                "delete the temp directory when finished (default: keep)",
            [メッセージID.ヘルプ_スレッド数] =
                "worker threads for read loading (default: number of logical processors)",
            [メッセージID.ヘルプ_言語] =
                "language of the messages (default: {0})",
            [メッセージID.ヘルプ_バージョン] =
                "show the version and exit",
            [メッセージID.ヘルプ_ヘルプ] =
                "show this text",
            [メッセージID.引き継ぎの統合開始] =
                "Merging {0:N0} carry-over sequence(s) from the previous k",
            [メッセージID.引き継ぎの統合進捗] =
                "[Carry-over] {0:N0}/{1:N0} sequence(s) merged",
            [メッセージID.引き継ぎの準備開始] =
                "Preparing carry-over sequences for the next k",
            [メッセージID.引き継ぎの準備完了] =
                "[Carry-over] Prepared {0:N0} sequence(s) for the next k.",
            [メッセージID.SuperRead橋渡し開始] =
                "[SuperRead] Bridging read pairs",
            [メッセージID.SuperRead橋渡し進捗] =
                "[SuperRead] {0:N0} pair(s) processed, {1:N0} bridged",
        };

        private static readonly Dictionary<メッセージID, string> _日本語 = new()
        {
            [メッセージID.デブルーイングラフの要約] =
                "[Debug] 厳密な de Bruijn ユニティググラフ: 有向辺 {0} 本、分岐頂点 {1} 個 / 全 {2} 頂点。",
            [メッセージID.先読みで解決した分岐数] =
                "[Debug] 先読み探索により、1手先の相互一意判定では決められなかった接合点をさらに {0} 箇所解決した。",
            [メッセージID.GFA出力完了] =
                "[Info] ユニティググラフを {0} に書き出した (GFA1)。",
            [メッセージID.分岐選択の重み内訳] =
                "[Debug] 分岐選択の重み: 単一リード由来の隣接 {0} 組 + ペアエンド由来 {1} 組。",
            [メッセージID.単純化の収束] =
                "[Debug] 単純化は {0} ラウンドで収束した。",
            [メッセージID.単純化の打ち切り] =
                "[Debug] 単純化はラウンド上限 ({0}) に達し、収束しきらずに打ち切った。",
            [メッセージID.バブル除去数] =
                "[Debug] 単純バブルの枝を合計 {0} 本除去した (配列自体は独立したコンティグとして残し、グラフ上の辺だけを外した)。",
            [メッセージID.反復解決数] =
                "[Debug] 反復解決: 短い反復 (<= {1}bp) を合計 {0} 箇所、跨いだリードペアを使って複製・分離した。",
            [メッセージID.辺選択の内訳] =
                "[Debug] 辺の選択: 出辺が1本だけの頂点 {0} 個、リード支持で解決した分岐 {1} 個、多コピー反復から出るため未解決のまま残した分岐 {2} 個 (反復内部のリードではコピーを区別できない)。",
            [メッセージID.反復通り抜けで棄却した結合数] =
                "[Debug] 未分離の多コピー反復を通り抜ける形になるため、{0} 件の結合を見送った (通すと反復のコピー間にある配列を飛ばしてしまう)。",
            [メッセージID.相互一意で残った結合数] =
                "[Debug] 相互一意の検査を通った有向マージは {0} 件 (無向では {1} 件)。",
            [メッセージID.環状コンティグあり] =
                "[Info] {0} 本のコンティグが環状に閉じた ({1})。環状に閉じたということは、そのレプリコン (染色体またはプラスミド) が端から端まで組み上がったということ。",
            [メッセージID.環状コンティグなし] =
                "[Info] 環状に閉じたコンティグは無い。どのレプリコンもまだ断片のまま。",
            [メッセージID.確定辺標本数] =
                "[Info] 確定した (実際に結合された) ユニティグ隣接から得たインサートサイズ標本: {0} 件。",
            [メッセージID.確定辺標本の中央値] =
                "[Info] 確定辺標本の中央値: {0} ({1} 件の標本。同一ユニティグ標本のような長さの偏りを受けない)。",
            [メッセージID.短すぎるユニティグの除外] =
                "[Warning] k-mer 長より短いユニティグ {0} 本を、マッピングの対象から外した。",
            [メッセージID.曖昧なkmer登録] =
                "[Warning] {0} 件の k-mer 登録が曖昧 (複数のユニティグで共有) だったため、マッピングでは無視する。",
            [メッセージID.ペア隣接候補数] =
                "[Info] ペアエンド由来の隣接候補: 辺 {0} 本 (支持したペアは延べ {1} 組)。",
            [メッセージID.同一ユニティグのペア向き集計] =
                "[Info] 同一ユニティグ内のペアの向き: 同方向={0}、逆方向={1}。インサートサイズ推定にはこのライブラリの向きとして '{2}' を採用する ({3} 件の標本)。",
            [メッセージID.同一ユニティグの断片長分布] =
                "[Info] 同一ユニティグ内の断片長分布: {0}。",
            [メッセージID.同一ユニティグの断片長中央値] =
                "[Info] 同一ユニティグ内の断片長中央値: {0} ({1} 件の標本。内側の距離にリード長を戻しているので、これは真の断片長にあたる。ただしユニティグが真のインサートサイズより短い場合は短めに偏りうる)。",
            [メッセージID.ペアリードIDの不一致] =
                "[Warning] この位置でペアのリードIDが一致しない (\"{0}\" と \"{1}\")。これ以降のリードではペアエンド由来の隣接検出が当てにならない可能性がある (単一リード由来の隣接検出には影響しない)。",
            [メッセージID.統合できる接合点なし] =
                "[Merge] 骨格の接合点を跨げた k は無かった。骨格をそのまま使う。",
            [メッセージID.統合した接合点数] =
                "[Merge] 他の k の配列を使って、骨格の接合点 {0} 箇所を繋いだ。",
            [メッセージID.候補一覧の見出し] =
                "[Multi-k] 候補アセンブリ (共通のアンカー k-mer 集合によるリファレンス無し評価):",
            [メッセージID.候補一覧の明細] =
                "[Multi-k]   k={0,3}: {1}{2}",
            [メッセージID.統計_ファイルなし] =
                "[Stats] {0}: (ファイルが見つからない: {1})",
            [メッセージID.統計] =
                "[Stats] {0}: {1}",
            [メッセージID.統計_長さで絞り込み] =
                "[Stats] {0} ({1}bp 以上、abyss-fac と比較できる条件): {2}",
            [メッセージID.検査_対象外のk長] =
                "[Check] {0}: 省略 (自己検査は k <= 64 のみ対応)。",
            [メッセージID.検査_取りこぼし] =
                "[Check] {0}: 信頼できる k-mer は {1:N0} 個。うち {2:N0} 個 ({3:0.00}%) はアセンブリにまったく現れない (削られたか、どの経路も通らなかった配列)。",
            [メッセージID.検査_出しすぎ] =
                "[Check] {0}: アセンブリ内の k-mer は延べ {1:N0} 個。うち {2:N0} 個 ({3:0.00}%) はカバレッジが許す以上のコピー数で、対象は {4:N0} 種類の k-mer -- ここが総延長の水増しにあたる。",
            [メッセージID.kmerインデックス構築開始] =
                "k-mer インデックスの構築を開始",
            [メッセージID.kmerカットオフ適用] =
                "k-mer カットオフを適用",
            [メッセージID.引き継ぎで追加したkmer数] =
                "[Carry-over] 前の k から、この k では観測されなかった k-mer を {0:N0} 個追加した。",
            [メッセージID.tip除去開始] =
                "短い tip を除去",
            [メッセージID.ユニティグ構築開始] =
                "ユニティグを構築",
            [メッセージID.グラフが複雑すぎる] =
                "[Warning] k={0} ではグラフが複雑すぎてアセンブリできない (ユニティグ数が {1} を超えた)。この k は見送る。",
            [メッセージID.リードのマッピング開始] =
                "リードをユニティグへマッピング",
            [メッセージID.ユニティグ結合開始] =
                "ユニティグを結合",
            [メッセージID.rMer検証の見送り] =
                "[Info] k={0} では反復の r-mer 検証を見送る: r-mer 長が {1}bp となり、32bp のパック上限を超えるため。",
            [メッセージID.コンティグ構築完了] =
                "コンティグの構築が完了",
            [メッセージID.スキャフォールディング開始] =
                "コンティグをスキャフォールディング",
            [メッセージID.ギャップ充填開始] =
                "スキャフォールドのギャップを充填",
            [メッセージID.kが薄すぎて省略] =
                "[Multi-k] k={0} は省略: 予測される k-mer カバレッジ {1:F1}x が {2:F1}x を下回り、使えるグラフにならない。",
            [メッセージID.kの開始見出し] =
                "[Multi-k] ===== k={0} でアセンブリ =====",
            [メッセージID.kでアセンブリできず] =
                "[Multi-k] k={0} ではアセンブリできなかった。この k は見送る。",
            [メッセージID.単一のkのみ成功] =
                "[Multi-k] アセンブリできた k が1つだけだったので、比較せずにそれを使う。",
            [メッセージID.アンカーkmer集合の構築] =
                "[Multi-k] 比較用の共通アンカー k-mer 集合 (k={0}) を構築",
            [メッセージID.アンカースペクトルが二峰でない] =
                "[Multi-k] アンカーの k-mer スペクトルが二峰にならないため、候補を比較できない。",
            [メッセージID.候補を評価できない] =
                "[Multi-k] 候補を評価できなかったため、最大の k を採用する。",
            [メッセージID.採用したk] =
                "[Multi-k] k={0} を採用した。",
            [メッセージID.統合開始] =
                "[Merge] 採用したアセンブリへ、他の k の配列を統合",
            [メッセージID.統合結果を評価できない] =
                "[Merge] 統合結果を評価できなかったため、採用済みのものを維持する。",
            [メッセージID.統合前の評価] =
                "[Merge]   統合前: {0}",
            [メッセージID.統合後の評価] =
                "[Merge]   統合後: {0}",
            [メッセージID.統合が骨格に勝てず] =
                "[Merge] 統合結果は採用済みのものに及ばなかったため、そのまま維持する。",
            [メッセージID.統合結果を採用] =
                "[Merge] 統合結果を採用する。",
            [メッセージID.エラー訂正_スペクトル構築] =
                "[ErrorCorrection] k-mer スペクトルを構築中...",
            [メッセージID.エラー訂正_訂正開始] =
                "[ErrorCorrection] リードを訂正中...",
            [メッセージID.エラー訂正_ファイル別統計] =
                "[ErrorCorrection] {0}: {1}/{2} リードを訂正 (訂正した塩基は延べ {3} 箇所)。",
            [メッセージID.リード読込完了] =
                "{1} から {0} リードを読み込んだ",
            [メッセージID.リード2の読込開始] =
                "ファイル2を読み込み",
            [メッセージID.前処理統計] =
                "[Preprocess] {1} 組のうち {0} 組でアダプタの読み抜けを検出 (重なった断片長まで切り詰めた)。相互に訂正した塩基は {2} 箇所。",
            [メッセージID.SuperRead統計] =
                "[SuperRead] {1} 組のうち {0} 組を1本の合成リードに繋ぎ、次の k へ引き継いだ。",
            [メッセージID.ギャップ充填_対象なし] =
                "[Info] ギャップ充填: 埋めるギャップが無い。",
            [メッセージID.ギャップ充填統計] =
                "[Info] ギャップ充填: {1} 箇所のうち {0} 箇所を実配列で埋めた ({2:N0}bp の N を置換)。{3} 箇所は経路が複数あるため、{4} 箇所は両側を繋ぐ経路がグラフ上に無いため N のまま。",
            [メッセージID.局所アセンブリ_対象なし] =
                "[Info] 局所アセンブリ: 対象となる残りギャップが無い。",
            [メッセージID.局所アセンブリ統計] =
                "[Info] 局所アセンブリ: 残り {1} 箇所のうち {0} 箇所を、局所的に組み直したリードで埋めた ({2:N0}bp の N を置換)。{3} 箇所は局所リードが集まらず、{4} 箇所は局所経路が複数あり、{5} 箇所は両側を繋ぐ局所経路が無かった。",
            [メッセージID.スキャフォールディング省略_インサートサイズ不明] =
                "[Info] スキャフォールディングを省略: インサートサイズが指定されておらず、マップされたペアからも推定できなかった。",
            [メッセージID.スキャフォールディング開始_インサートサイズ] =
                "[Info] インサートサイズ {0} でスキャフォールディング",
            [メッセージID.スキャフォールディング省略_コンティグなし] =
                "[Info] スキャフォールディングを省略: コンティグが見つからなかった。",
            [メッセージID.内部を指したペア候補] =
                "[Info] {0} 件のペアエンド候補は、既に結合済みのコンティグ内部にあるユニティグを指していたため見送った (端点はコンティグ構築の時点で解決済み)。",
            [メッセージID.未配置を指したペア候補] =
                "[Info] {0} 件のペアエンド候補は、どのコンティグにも配置されなかったユニティグ (短すぎる等) を指していたため見送った。",
            [メッセージID.閾値後のスキャフォールド辺] =
                "[Info] 閾値処理後に確定したスキャフォールド辺: {0} 本。うち {1} 本は相互一意の検査で棄却し、{2} 本を採用した。",
            [メッセージID.スキャフォールド出力完了] =
                "[Info] {0} 本のスキャフォールド (総延長 {1}) を {2} に書き出した",
            [メッセージID.インサートサイズ推定_同一ユニティグ] =
                "[Info] インサートサイズを {1} 組の同一ユニティグ標本から {0} と推定した (中央値。ユニティグ N50 {2} が推定値の {3} 倍以上あるため、短い断片への打ち切りバイアスは効かない)。",
            [メッセージID.インサートサイズ推定_確定辺] =
                "[Info] インサートサイズを {1} 組の確定辺標本から {0} と推定した (中央値。ユニティグが十分に長くなく同一ユニティグ標本は偏るため、こちらを優先した)。",
            [メッセージID.インサートサイズ推定_標本不足] =
                "[Info] インサートサイズの自動推定には最低 {0} 件の標本が必要だが、確定辺由来 {1} 件、全体でも {2} 件しか集まらなかった。",
            [メッセージID.インサートサイズ推定_全標本] =
                "[Info] インサートサイズを {1} 組の標本から {0} と推定した (中央値。確定辺由来の標本が {2} 件と少なかったため全標本に戻した。短めに偏る可能性がある)。",
            [メッセージID.単一コピー基準値] =
                "[Info] 単一コピーのカバレッジ基準値を {0:0.#} と推定した (長さで重み付けした中央値)。",
            [メッセージID.コピー数の要約] =
                "[Info] 推定したコピー数 -- {0}",
            [メッセージID.反復配列の割合] =
                "[Info] 多コピー (反復) の量: {1:N0}bp 中 {0:N0}bp (アセンブリの {2:0.0}% が2回以上現れる配列)。",
            [メッセージID.グラフ単純化の反復] =
                "[GraphSimplifier] 反復 {0}: ユニティグ {1} 本を検査し (tip 閾値 {2}bp 未満、カバレッジ基準値 {3:0.#})、tip を {4} 本除去、{6} 本のユニティグ端から低カバレッジの k-mer を {5} 個削った。",
            [メッセージID.rMer検証による棄却] =
                "[Debug] r-mer 検証により、接合点を実際に跨いだリードで裏付けが取れなかった反復の複製を {0} 件棄却した (ペア数の証拠だけなら複製していた)。",
            [メッセージID.例外を無視_見出し] =
                "[Warning] 以下の例外は無視して続行する。",
            [メッセージID.例外を無視_メソッド] =
                "          メソッド: {0}",
            [メッセージID.停止_見出し] =
                "[Error] 処理を中断した。",
            [メッセージID.停止_メソッド] =
                "        メソッド: {0}",
            [メッセージID.タイムスタンプ] =
                "[Log] {0:yyyy/MM/dd HH:mm:ss}",
            [メッセージID.Phred_明示指定と不一致] =
                "[Warning] クオリティ列は Phred{0} に見えるが、-p {1} が明示指定されている。指定値を優先する。本当に Phred{3} なら -p {2} で再実行すること。",
            [メッセージID.Phred_自動判定] =
                "[Info] クオリティ列から Phred オフセットを {0} と自動判定した (観測した ASCII 範囲 [{1}, {2}])。-p を明示すれば上書きできる。",
            [メッセージID.Phred_検査の警告] =
                "[Warning] {0} の Phred エンコーディング検査: {1}",
            [メッセージID.kmerカットオフ_混合モデル] =
                "[Info] k-mer カットオフを2成分混合モデル (エラー + ゲノム) から {0} と自動選択した (単一コピー平均 {1:0.#}、信頼下限 {2}、EM 反復 {3} 回)。-kc を明示すれば上書きできる。",
            [メッセージID.kmerカットオフ_谷が不明] =
                "[Info] ヒストグラムの谷をはっきり特定できなかった (このカバレッジではスペクトルが二峰にならない可能性がある)。-kc {0} のままにする。",
            [メッセージID.kmerカットオフ_スペクトル] =
                "[Info] k-mer カットオフを k-mer スペクトルから {0} と自動選択した。-kc を明示すれば上書きできる。",
            [メッセージID.kmerヒストグラム] =
                "[Info] k-mer 出現回数のヒストグラム (出現回数:異なる k-mer 数): {0}",
            [メッセージID.スペクトルの谷が不明] =
                "[Info] ヒストグラムの谷をはっきり特定できなかった (このカバレッジではスペクトルが二峰にならない可能性がある)。",
            [メッセージID.スペクトルの谷とピーク] =
                "[Info] k-mer スペクトル: 谷は出現回数 {0} ({1} 種類)、単一コピーのピークは出現回数 {2} ({3} 種類)",
            [メッセージID.ゲノムサイズとカバレッジの推定] =
                "[Info] 推定ゲノムサイズ: {0:N0} bp、推定カバレッジ: {1}",
            [メッセージID.k自動選択_リード長不明] =
                "[Info] リード長を標本抽出できなかった。-k {0} のままにする。",
            [メッセージID.k自動選択_kが長すぎる] =
                "[Warning] -k {0} は観測されたリード長 ({1} bp) より短くないため、リードから k-mer を取り出せない。-k を下げること。",
            [メッセージID.k自動選択_リード長が短い] =
                "[Info] 観測されたリード長 ({0} bp) は k を自動選択するには短すぎる。-k {1} のままにする。",
            [メッセージID.k自動選択] =
                "[Info] 観測されたリード長 ({1} bp) から k を {0} と自動選択した。-k を明示すれば上書きできる。",
            [メッセージID.開始kmerの探索] =
                "開始 k-mer を探索",
            [メッセージID.一時ディレクトリを残した] =
                "[Info] k ごとのアセンブリを {0} に残した ({1} を付けると削除する)。",
            [メッセージID.リード長の観測値] =
                "[Info] リード長 (標本の中央値): {0} bp",
            [メッセージID.一時ディレクトリが既にある] =
                "{0} は既に存在する",
            [メッセージID.パスの確認] =
                "パスを確認すること！",
            [メッセージID.前処理省略_ペアなし] =
                "[Warning] -pp はペアエンドのリードが必要 (read2 が指定されていない)。前処理を省略する。",
            [メッセージID.前処理開始] =
                "リードを前処理 (アダプタ除去 + ペア相互訂正)",
            [メッセージID.エラー訂正開始] =
                "アセンブリ前にリードを訂正",
            [メッセージID.開発中] =
                "開発中！",
            [メッセージID.コンティグ総延長] =
                "コンティグの総延長 : {0}",
            [メッセージID.リードファイルのパス] =
                "{0}",
            [メッセージID.試すk一覧] =
                "[Multi-k] 試す k = {0}",
            [メッセージID.リード読込の進捗] =
                "{0} リードを読み込み済み",
            [メッセージID.リード1の読込開始] =
                "ファイル1を読み込み",
            [メッセージID.単一リードの読込開始] =
                "ファイルを読み込み",
            [メッセージID.スキャフォールド候補辺数] =
                "[Info] スキャフォールドの候補辺 (コンティグ単位、閾値処理の前): {0}{1}",
            [メッセージID.理想本数モデルあり] =
                "。理想本数モデルは利用可能",
            [メッセージID.理想本数モデルなし] =
                "。理想本数モデルは利用不可",
            [メッセージID.Phred_ファイル間で不一致] =
                "[Warning] 2つのリードファイルで Phred オフセットの推定が食い違った (read1 -> {0}、read2 -> {1})。-p {2} をそのまま使う。",
            [メッセージID.Phred_未確定] =
                "不明",
            [メッセージID.kmer種類数] =
                "k-mer の種類数: {0}",
            [メッセージID.採用kmer数] =
                "採用した k-mer: {0}",
            [メッセージID.アセンブリ不能] =
                "\nこのゲノムは複雑すぎてアセンブリできない...\nパラメータを見直すこと！\n",
            [メッセージID.概要_説明] =
                "Tsumiki はゲノムアセンブラです。",
            [メッセージID.概要_作者] =
                "作者: {0}",
            [メッセージID.概要_バージョン] =
                "バージョン: {0}",
            [メッセージID.ヘルプ_使い方] =
                "使い方: tsumiki {0} <path> [{1} <path>] [options]",
            [メッセージID.ヘルプ節_入力] =
                "# 入力",
            [メッセージID.ヘルプ節_kmerと品質] =
                "# k-mer / 品質",
            [メッセージID.ヘルプ節_ペアエンド] =
                "# ペアエンド / スキャフォールディング",
            [メッセージID.ヘルプ節_前処理] =
                "# 前処理",
            [メッセージID.ヘルプ節_マルチk] =
                "# マルチ k アセンブリ",
            [メッセージID.ヘルプ節_反復の安全策] =
                "# 反復解決の安全策",
            [メッセージID.ヘルプ節_出力とその他] =
                "# 出力 / その他",
            [メッセージID.ヘルプ_リード1] =
                "forward 側 fastq(.gz) のパス (必須。シングルエンドの場合もこちらに指定する)",
            [メッセージID.ヘルプ_リード2] =
                "reverse 側 fastq(.gz) のパス",
            [メッセージID.ヘルプ_曖昧塩基] =
                "リード中の曖昧塩基 (N など) を許容する (既定: false)",
            [メッセージID.ヘルプ_k長] =
                "k-mer 長。カンマ区切り (例 31,63,95) にすると {0} と同じく全て試して最良を残す (既定: リード長からの自動選択。上限 {1})",
            [メッセージID.ヘルプ_kmerカットオフ] =
                "信頼する k-mer の最小出現回数 (既定: k-mer スペクトルからの自動選択)",
            [メッセージID.ヘルプ_Phred] =
                "phred スコアの基準値 ({0}) (既定: {1})",
            [メッセージID.ヘルプ_クオリティカットオフ] =
                "信頼する塩基クオリティの下限 (既定: {0})",
            [メッセージID.ヘルプ_メモリ予算] =
                "k-mer カウントに使うメモリ量。例: 2G, 512M (既定: {0})",
            [メッセージID.ヘルプ_インサートサイズ] =
                "期待するインサートサイズ (既定: マップされたペアからの自動推定)",
            [メッセージID.ヘルプ_ペア結合閾値] =
                "スキャフォールド辺を認める優勢比の下限 (既定: {0})",
            [メッセージID.ヘルプ_ペア支持数閾値] =
                "短い反復を解決するのに必要なリードペア支持数の下限 (既定: {0})",
            [メッセージID.ヘルプ_積極性モード] =
                "{0}/{1} のプリセット。完全性と正確性のどちらに倒すかを決める (既定: {2})",
            [メッセージID.ヘルプ_前処理] =
                "アダプタの読み抜けを切り、ペアの重なりを使って低品質の塩基を相互訂正する (既定: false)",
            [メッセージID.ヘルプ_エラー訂正] =
                "アセンブリ前に k-mer スペクトルからリードを訂正する (既定: false)",
            [メッセージID.ヘルプ_マルチk] =
                "複数の k でアセンブリし最良のものを残す (既定: false。実行時間は最大 {0} 倍)",
            [メッセージID.ヘルプ_引き継ぎなし] =
                "前の k の配列を次の k へ引き継がない (既定: 引き継ぐ)",
            [メッセージID.ヘルプ_SuperRead] =
                "ペアエンドの場合、ペアを橋渡しして合成ロングリードを作り、それも引き継ぐ (既定: false)",
            [メッセージID.ヘルプ_マージ] =
                "最良の k が繋げなかった接合点に、他の k の配列を継ぎ足す (既定: false。誤アセンブリが増えることがある)",
            [メッセージID.ヘルプ_反復rMer検証] =
                "短い反復を複製する前に、生リードの r-mer による裏付けを必須にする (既定: false)",
            [メッセージID.ヘルプ_局所アセンブリ] =
                "埋まらなかったスキャフォールドのギャップを、その両端に付いたリードだけで組み直す。{0} より安全 (既定: false)",
            [メッセージID.ヘルプ_GFA出力] =
                "Bandage などで見るために {0} (ユニティググラフ、GFA1) も書き出す (既定: false)",
            [メッセージID.ヘルプ_一時ディレクトリ] =
                "一時ディレクトリ。k ごとのアセンブリはここに残る (既定: {0})",
            [メッセージID.ヘルプ_一時ディレクトリ削除] =
                "実行後に一時ディレクトリを削除する (既定: 残す)",
            [メッセージID.ヘルプ_スレッド数] =
                "リード読み込みに使うワーカースレッド数 (既定: 論理プロセッサ数)",
            [メッセージID.ヘルプ_言語] =
                "メッセージの言語 (既定: {0})",
            [メッセージID.ヘルプ_バージョン] =
                "バージョンを表示して終了する",
            [メッセージID.ヘルプ_ヘルプ] =
                "このテキストを表示する",
            [メッセージID.引き継ぎの統合開始] =
                "前の k から引き継いだ配列 {0:N0} 本を統合中",
            [メッセージID.引き継ぎの統合進捗] =
                "[Carry-over] {1:N0} 本中 {0:N0} 本を統合済み",
            [メッセージID.引き継ぎの準備開始] =
                "次の k へ引き継ぐ配列を準備中",
            [メッセージID.引き継ぎの準備完了] =
                "[Carry-over] 次の k へ引き継ぐ配列を {0:N0} 本用意した。",
            [メッセージID.SuperRead橋渡し開始] =
                "[SuperRead] リードペアを橋渡し中",
            [メッセージID.SuperRead橋渡し進捗] =
                "[SuperRead] {0:N0} 組を処理、うち {1:N0} 組を橋渡し",
        };

        private static readonly Dictionary<メッセージID, string> _中国語 = new()
        {
            [メッセージID.デブルーイングラフの要約] =
                "[Debug] 精确的 de Bruijn unitig 图：有向边 {0} 条，分支顶点 {1} 个 / 共 {2} 个顶点。",
            [メッセージID.先読みで解決した分岐数] =
                "[Debug] 束搜索前瞻额外解决了 {0} 处单步互唯一性判定无法决定的连接点。",
            [メッセージID.GFA出力完了] =
                "[Info] 已将 unitig 图写入 {0}（GFA1）。",
            [メッセージID.分岐選択の重み内訳] =
                "[Debug] 分支选择权重：单条 read 相邻 {0} 对 + 双端 read {1} 对。",
            [メッセージID.単純化の収束] =
                "[Debug] 简化在第 {0} 轮收敛。",
            [メッセージID.単純化の打ち切り] =
                "[Debug] 简化达到轮数上限（{0}），未完全收敛即停止。",
            [メッセージID.バブル除去数] =
                "[Debug] 共移除 {0} 条简单气泡分支（序列本身作为独立 contig 保留，仅移除图中的边）。",
            [メッセージID.反復解決数] =
                "[Debug] 重复解析：共有 {0} 处短重复（<= {1}bp）借助跨越它们的 read 对完成复制与解缠。",
            [メッセージID.辺選択の内訳] =
                "[Debug] 边选择：出边唯一的顶点 {0} 个，凭 read 支持解决的分支 {1} 个，因离开多拷贝重复而未解决的分支 {2} 个（重复内部的 read 无法区分各拷贝）。",
            [メッセージID.反復通り抜けで棄却した結合数] =
                "[Debug] 有 {0} 处连接被拒绝，因为它们会穿过尚未解缠的多拷贝重复（这样会跳过重复各拷贝之间的序列）。",
            [メッセージID.相互一意で残った結合数] =
                "[Debug] 通过互唯一性检查的有向合并 {0} 处（无向 {1} 处）。",
            [メッセージID.環状コンティグあり] =
                "[Info] 有 {0} 条 contig 闭合成环（{1}）。闭合成环意味着该复制子（染色体或质粒）已被完整组装。",
            [メッセージID.環状コンティグなし] =
                "[Info] 没有 contig 闭合成环；所有复制子仍是片段状态。",
            [メッセージID.確定辺標本数] =
                "[Info] 来自已确定（实际连接）unitig 相邻关系的插入片段长度样本：{0} 个。",
            [メッセージID.確定辺標本の中央値] =
                "[Info] 已确定边样本的中位数：{0}（共 {1} 个样本；不受同一 unitig 样本的长度偏差影响）。",
            [メッセージID.短すぎるユニティグの除外] =
                "[Warning] 有 {0} 条短于 k-mer 长度的 unitig 在比对中被跳过。",
            [メッセージID.曖昧なkmer登録] =
                "[Warning] 有 {0} 处 k-mer 登记存在歧义（被多条 unitig 共享），比对时将忽略。",
            [メッセージID.ペア隣接候補数] =
                "[Info] 检测到双端相邻候选：{0} 条边（共 {1} 对支持的 read 对）。",
            [メッセージID.同一ユニティグのペア向き集計] =
                "[Info] 同一 unitig 内 read 对的方向统计：同向={0}，反向={1}。插入片段长度估计采用该文库观测到的方向 '{2}'（{3} 个样本）。",
            [メッセージID.同一ユニティグの断片長分布] =
                "[Info] 同一 unitig 内片段长度分布：{0}。",
            [メッセージID.同一ユニティグの断片長中央値] =
                "[Info] 同一 unitig 内片段长度中位数：{0}（共 {1} 个样本；已把 read 长度加回内部距离，因此这是真实的片段长度。若 unitig 短于真实插入片段，仍可能偏短）。",
            [メッセージID.ペアリードIDの不一致] =
                "[Warning] 此处配对 read 的 ID 不一致（\"{0}\" 与 \"{1}\"）。此后的 read 其双端相邻检测可能不可靠；单条 read 的相邻检测不受影响。",
            [メッセージID.統合できる接合点なし] =
                "[Merge] 没有其他 k 跨越骨架的连接点；保持原样。",
            [メッセージID.統合した接合点数] =
                "[Merge] 使用其他 k 的序列连接了骨架的 {0} 处连接点。",
            [メッセージID.候補一覧の見出し] =
                "[Multi-k] 候选组装（基于共同锚定 k-mer 集合的无参考评估）：",
            [メッセージID.候補一覧の明細] =
                "[Multi-k]   k={0,3}: {1}{2}",
            [メッセージID.統計_ファイルなし] =
                "[Stats] {0}：（找不到文件：{1}）",
            [メッセージID.統計] =
                "[Stats] {0}：{1}",
            [メッセージID.統計_長さで絞り込み] =
                "[Stats] {0}（>={1}bp，可与 abyss-fac 对比）：{2}",
            [メッセージID.検査_対象外のk長] =
                "[Check] {0}：已跳过（自检仅支持 k <= 64）。",
            [メッセージID.検査_取りこぼし] =
                "[Check] {0}：可信 k-mer 共 {1:N0} 个；其中 {2:N0} 个（{3:0.00}%）完全未出现在组装中（被修剪掉或没有任何路径到达的序列）。",
            [メッセージID.検査_出しすぎ] =
                "[Check] {0}：组装中 k-mer 实例共 {1:N0} 个；其中 {2:N0} 个（{3:0.00}%）超出覆盖度所能支持的拷贝数，涉及 {4:N0} 种不同的 k-mer —— 这部分即为总长度的虚增。",
            [メッセージID.kmerインデックス構築開始] =
                "开始构建 k-mer 索引",
            [メッセージID.kmerカットオフ適用] =
                "正在应用 k-mer 截断阈值",
            [メッセージID.引き継ぎで追加したkmer数] =
                "[Carry-over] 从上一个 k 补充了 {0:N0} 个当前 k 未观测到的 k-mer。",
            [メッセージID.tip除去開始] =
                "正在修剪短 tip",
            [メッセージID.ユニティグ構築開始] =
                "正在构建 unitig",
            [メッセージID.グラフが複雑すぎる] =
                "[Warning] k={0} 时图过于复杂，无法组装（unitig 数超过 {1}）。跳过该 k。",
            [メッセージID.リードのマッピング開始] =
                "正在将 read 比对到 unitig",
            [メッセージID.ユニティグ結合開始] =
                "正在合并 unitig",
            [メッセージID.rMer検証の見送り] =
                "[Info] k={0} 时跳过重复的 r-mer 验证：r-mer 长度将为 {1}bp，超过 32bp 的打包上限。",
            [メッセージID.コンティグ構築完了] =
                "contig 构建完成",
            [メッセージID.スキャフォールディング開始] =
                "正在对 contig 进行 scaffold",
            [メッセージID.ギャップ充填開始] =
                "正在填补 scaffold 的空缺",
            [メッセージID.kが薄すぎて省略] =
                "[Multi-k] 跳过 k={0}：预测的 k-mer 覆盖度 {1:F1}x 低于 {2:F1}x，无法得到可用的图。",
            [メッセージID.kの開始見出し] =
                "[Multi-k] ===== 使用 k={0} 进行组装 =====",
            [メッセージID.kでアセンブリできず] =
                "[Multi-k] k={0} 未产出组装结果；跳过。",
            [メッセージID.単一のkのみ成功] =
                "[Multi-k] 只有一个 k 产出了组装结果，直接采用，不再比较。",
            [メッセージID.アンカーkmer集合の構築] =
                "[Multi-k] 正在构建用于比较的共同锚定 k-mer 集合（k={0}）",
            [メッセージID.アンカースペクトルが二峰でない] =
                "[Multi-k] 锚定 k-mer 谱不是双峰，无法比较候选。",
            [メッセージID.候補を評価できない] =
                "[Multi-k] 无法评估候选，回退到最大的 k。",
            [メッセージID.採用したk] =
                "[Multi-k] 已选定 k={0}。",
            [メッセージID.統合開始] =
                "[Merge] 正在将其他 k 的序列合并到选定的组装中",
            [メッセージID.統合結果を評価できない] =
                "[Merge] 无法评估合并结果，保留已选定的组装。",
            [メッセージID.統合前の評価] =
                "[Merge]   合并前：{0}",
            [メッセージID.統合後の評価] =
                "[Merge]   合并后：{0}",
            [メッセージID.統合が骨格に勝てず] =
                "[Merge] 合并结果未优于已选定的组装，保持原样。",
            [メッセージID.統合結果を採用] =
                "[Merge] 采用合并结果。",
            [メッセージID.エラー訂正_スペクトル構築] =
                "[ErrorCorrection] 正在构建 k-mer 谱……",
            [メッセージID.エラー訂正_訂正開始] =
                "[ErrorCorrection] 正在校正 read……",
            [メッセージID.エラー訂正_ファイル別統計] =
                "[ErrorCorrection] {0}：已校正 {1}/{2} 条 read（共校正 {3} 个碱基）。",
            [メッセージID.リード読込完了] =
                "已从 {1} 读入 {0} 条 read",
            [メッセージID.リード2の読込開始] =
                "正在读取文件 2",
            [メッセージID.前処理統計] =
                "[Preprocess] {1} 对中有 {0} 对检测到接头读通（已修剪至重叠片段长度）；相互校正碱基 {2} 个。",
            [メッセージID.SuperRead統計] =
                "[SuperRead] {1} 对中有 {0} 对被桥接为单条合成 read 并传递到下一个 k。",
            [メッセージID.ギャップ充填_対象なし] =
                "[Info] 空缺填补：没有需要填补的空缺。",
            [メッセージID.ギャップ充填統計] =
                "[Info] 空缺填补：{1} 处空缺中有 {0} 处用真实序列填补（替换了 {2:N0}bp 的 N）。{3} 处因存在多条路径、{4} 处因图中没有连接两侧的路径而保留为 N。",
            [メッセージID.局所アセンブリ_対象なし] =
                "[Info] 局部组装：没有符合条件的剩余空缺。",
            [メッセージID.局所アセンブリ統計] =
                "[Info] 局部组装：剩余 {1} 处空缺中有 {0} 处用局部重组的 read 填补（替换了 {2:N0}bp 的 N）。{3} 处没有任何局部 read，{4} 处存在多条局部路径，{5} 处没有连接两侧的局部路径。",
            [メッセージID.スキャフォールディング省略_インサートサイズ不明] =
                "[Info] 跳过 scaffold：未指定插入片段长度，也无法从比对上的 read 对估计。",
            [メッセージID.スキャフォールディング開始_インサートサイズ] =
                "[Info] 以插入片段长度 {0} 进行 scaffold",
            [メッセージID.スキャフォールディング省略_コンティグなし] =
                "[Info] 跳过 scaffold：未找到 contig。",
            [メッセージID.内部を指したペア候補] =
                "[Info] 有 {0} 个双端候选指向已连接 contig 内部的 unitig，已跳过（端点在 contig 构建时已解决）。",
            [メッセージID.未配置を指したペア候補] =
                "[Info] 有 {0} 个双端候选指向未被放入任何 contig 的 unitig（如过短），已跳过。",
            [メッセージID.閾値後のスキャフォールド辺] =
                "[Info] 阈值处理后确定的 scaffold 边：{0} 条；其中 {1} 条被互唯一性检查拒绝，保留 {2} 条。",
            [メッセージID.スキャフォールド出力完了] =
                "[Info] 已将 {0} 条 scaffold（总长 {1}）写入 {2}",
            [メッセージID.インサートサイズ推定_同一ユニティグ] =
                "[Info] 依据 {1} 对同一 unitig 样本，将插入片段长度估计为 {0}（中位数；unitig N50 为 {2}，达到估计值的 {3} 倍以上，因此不受短片段截断偏差影响）。",
            [メッセージID.インサートサイズ推定_確定辺] =
                "[Info] 依据 {1} 对已确定边样本，将插入片段长度估计为 {0}（中位数；因 unitig 长度不足会使同一 unitig 样本产生偏差，故优先采用）。",
            [メッセージID.インサートサイズ推定_標本不足] =
                "[Info] 插入片段长度自动估计至少需要 {0} 个样本，但仅收集到已确定边样本 {1} 个、总计 {2} 个。",
            [メッセージID.インサートサイズ推定_全標本] =
                "[Info] 依据 {1} 对样本，将插入片段长度估计为 {0}（中位数；已确定边样本仅 {2} 个，过少，故回退到全部样本，可能偏短）。",
            [メッセージID.単一コピー基準値] =
                "[Info] 单拷贝覆盖度基线估计为 {0:0.#}（按长度加权的中位数）。",
            [メッセージID.コピー数の要約] =
                "[Info] 估计的拷贝数 —— {0}",
            [メッセージID.反復配列の割合] =
                "[Info] 多拷贝（重复）含量：{1:N0}bp 中的 {0:N0}bp（组装的 {2:0.0}% 是出现一次以上的序列）。",
            [メッセージID.グラフ単純化の反復] =
                "[GraphSimplifier] 第 {0} 轮：检查 unitig {1} 条（tip 阈值 < {2}bp，覆盖度基线 {3:0.#}），移除 tip {4} 条，从 {6} 条 unitig 的末端修剪低覆盖 k-mer {5} 个。",
            [メッセージID.rMer検証による棄却] =
                "[Debug] r-mer 验证否决了 {0} 处重复复制，其胜出配对未被真正跨越连接点的 read 证实（仅凭 read 对计数证据本会复制它们）。",
            [メッセージID.例外を無視_見出し] =
                "[Warning] 以下异常将被忽略。",
            [メッセージID.例外を無視_メソッド] =
                "          方法：{0}",
            [メッセージID.停止_見出し] =
                "[Error] 程序已停止。",
            [メッセージID.停止_メソッド] =
                "        方法：{0}",
            [メッセージID.タイムスタンプ] =
                "[Log] {0:yyyy/MM/dd HH:mm:ss}",
            [メッセージID.Phred_明示指定と不一致] =
                "[Warning] 质量字符串看起来是 Phred{0}，但显式指定了 -p {1}。将采用显式值；若数据确为 Phred{3}，请以 -p {2} 重新运行。",
            [メッセージID.Phred_自動判定] =
                "[Info] 依据质量字符串自动判定 Phred 偏移为 {0}（观测到的 ASCII 范围 [{1}, {2}]）。可用 -p 显式覆盖。",
            [メッセージID.Phred_検査の警告] =
                "[Warning] {0} 的 Phred 编码检查：{1}",
            [メッセージID.kmerカットオフ_混合モデル] =
                "[Info] 依据双成分混合模型（错误 + 基因组）自动选择 k-mer 截断阈值为 {0}（单拷贝均值 {1:0.#}，可信下限 {2}，EM 迭代 {3} 次）。可用 -kc 显式覆盖。",
            [メッセージID.kmerカットオフ_谷が不明] =
                "[Info] 无法明确识别直方图谷值（在此覆盖度下谱可能不是双峰）；保持 -kc {0}。",
            [メッセージID.kmerカットオフ_スペクトル] =
                "[Info] 依据 k-mer 谱自动选择 k-mer 截断阈值为 {0}。可用 -kc 显式覆盖。",
            [メッセージID.kmerヒストグラム] =
                "[Info] k-mer 计数直方图（计数:不同 k-mer 数）：{0}",
            [メッセージID.スペクトルの谷が不明] =
                "[Info] 无法明确识别直方图谷值（在此覆盖度下谱可能不是双峰）。",
            [メッセージID.スペクトルの谷とピーク] =
                "[Info] k-mer 谱：谷位于计数 {0}（{1} 种），单拷贝峰位于计数 {2}（{3} 种）",
            [メッセージID.ゲノムサイズとカバレッジの推定] =
                "[Info] 估计基因组大小：{0:N0} bp，估计覆盖度：{1}",
            [メッセージID.k自動選択_リード長不明] =
                "[Info] 无法采样得到 read 长度；保持 -k {0}。",
            [メッセージID.k自動選択_kが長すぎる] =
                "[Warning] -k {0} 不短于观测到的 read 长度（{1} bp），无法从 read 中提取 k-mer。请调小 -k。",
            [メッセージID.k自動選択_リード長が短い] =
                "[Info] 观测到的 read 长度（{0} bp）过短，无法自动选择 k；保持 -k {1}。",
            [メッセージID.k自動選択] =
                "[Info] 依据观测到的 read 长度（{1} bp）自动选择 k 为 {0}。可用 -k 显式覆盖。",
            [メッセージID.開始kmerの探索] =
                "正在搜索起始 k-mer",
            [メッセージID.一時ディレクトリを残した] =
                "[Info] 各 k 的组装结果保留在 {0}（加上 {1} 可删除）。",
            [メッセージID.リード長の観測値] =
                "[Info] read 长度（采样中位数）：{0} bp",
            [メッセージID.一時ディレクトリが既にある] =
                "{0} 已存在",
            [メッセージID.パスの確認] =
                "请检查路径！",
            [メッセージID.前処理省略_ペアなし] =
                "[Warning] -pp 需要双端 read（未设置 read2）。跳过预处理。",
            [メッセージID.前処理開始] =
                "正在预处理 read（去接头 + 配对互校正）",
            [メッセージID.エラー訂正開始] =
                "组装前正在校正 read",
            [メッセージID.開発中] =
                "开发中！",
            [メッセージID.コンティグ総延長] =
                "contig 总长度：{0}",
            [メッセージID.リードファイルのパス] =
                "{0}",
            [メッセージID.試すk一覧] =
                "[Multi-k] 将尝试 k = {0}",
            [メッセージID.リード読込の進捗] =
                "已读入 {0} 条 read",
            [メッセージID.リード1の読込開始] =
                "正在读取文件 1",
            [メッセージID.単一リードの読込開始] =
                "正在读取文件",
            [メッセージID.スキャフォールド候補辺数] =
                "[Info] scaffold 候选边（contig 级，阈值处理前）：{0}{1}",
            [メッセージID.理想本数モデルあり] =
                "；理想计数模型可用",
            [メッセージID.理想本数モデルなし] =
                "；理想计数模型不可用",
            [メッセージID.Phred_ファイル間で不一致] =
                "[Warning] 两个 read 文件推断出的 Phred 偏移不一致（read1 -> {0}，read2 -> {1}）。保持 -p {2} 不变。",
            [メッセージID.Phred_未確定] =
                "未确定",
            [メッセージID.kmer種類数] =
                "k-mer 种类数：{0}",
            [メッセージID.採用kmer数] =
                "采用的 k-mer：{0}",
            [メッセージID.アセンブリ不能] =
                "\n该基因组过于复杂，无法组装……\n请调整参数！\n",
            [メッセージID.概要_説明] =
                "Tsumiki 是一个基因组组装工具。",
            [メッセージID.概要_作者] =
                "作者：{0}",
            [メッセージID.概要_バージョン] =
                "版本：{0}",
            [メッセージID.ヘルプ_使い方] =
                "用法：tsumiki {0} <path> [{1} <path>] [options]",
            [メッセージID.ヘルプ節_入力] =
                "# 输入",
            [メッセージID.ヘルプ節_kmerと品質] =
                "# k-mer / 质量",
            [メッセージID.ヘルプ節_ペアエンド] =
                "# 双端 / scaffold",
            [メッセージID.ヘルプ節_前処理] =
                "# 预处理",
            [メッセージID.ヘルプ節_マルチk] =
                "# 多 k 组装",
            [メッセージID.ヘルプ節_反復の安全策] =
                "# 重复解析的安全措施",
            [メッセージID.ヘルプ節_出力とその他] =
                "# 输出 / 其他",
            [メッセージID.ヘルプ_リード1] =
                "forward 端 fastq(.gz) 路径（必填；单端 read 也使用此项）",
            [メッセージID.ヘルプ_リード2] =
                "reverse 端 fastq(.gz) 路径",
            [メッセージID.ヘルプ_曖昧塩基] =
                "允许 read 中的模糊碱基（如 N）（默认：false）",
            [メッセージID.ヘルプ_k長] =
                "k-mer 长度；用逗号分隔（如 31,63,95）时与 {0} 一样逐个尝试并保留最优（默认：依据 read 长度自动选择，上限 {1}）",
            [メッセージID.ヘルプ_kmerカットオフ] =
                "可信 k-mer 的最小计数（默认：依据 k-mer 谱自动选择）",
            [メッセージID.ヘルプ_Phred] =
                "phred 分值基准（{0}）（默认：{1}）",
            [メッセージID.ヘルプ_クオリティカットオフ] =
                "可信碱基质量的下限（默认：{0}）",
            [メッセージID.ヘルプ_メモリ予算] =
                "k-mer 计数的内存预算，如 2G、512M（默认：{0}）",
            [メッセージID.ヘルプ_インサートサイズ] =
                "期望的插入片段长度（默认：从比对上的 read 对自动估计）",
            [メッセージID.ヘルプ_ペア結合閾値] =
                "接受 scaffold 边所需的最小优势比（默认：{0}）",
            [メッセージID.ヘルプ_ペア支持数閾値] =
                "解决短重复所需的最小 read 对支持数（默认：{0}）",
            [メッセージID.ヘルプ_積極性モード] =
                "{0}/{1} 的预设：在完整性与准确性之间取舍（默认：{2}）",
            [メッセージID.ヘルプ_前処理] =
                "修剪接头读通，并借助配对重叠互相校正低质量碱基（默认：false）",
            [メッセージID.ヘルプ_エラー訂正] =
                "组装前依据 k-mer 谱校正 read（默认：false）",
            [メッセージID.ヘルプ_マルチk] =
                "使用多个 k 组装并保留最优结果（默认：false；运行时间最多 {0} 倍）",
            [メッセージID.ヘルプ_引き継ぎなし] =
                "不将上一个 k 的序列传递到下一个 k（默认：传递）",
            [メッセージID.ヘルプ_SuperRead] =
                "对于双端 read，将配对桥接为合成长 read 并一并传递（默认：false）",
            [メッセージID.ヘルプ_マージ] =
                "将其他 k 的序列拼接到最优 k 未能连接的接合点（默认：false；可能增加错误组装）",
            [メッセージID.ヘルプ_反復rMer検証] =
                "在复制短重复前要求原始 read 的 r-mer 支持（默认：false）",
            [メッセージID.ヘルプ_局所アセンブリ] =
                "用比对到空缺两端附近的 read 重新组装未闭合的 scaffold 空缺；比 {0} 更安全（默认：false）",
            [メッセージID.ヘルプ_GFA出力] =
                "同时写出 {0}（unitig 图，GFA1）以便用 Bandage 等查看（默认：false）",
            [メッセージID.ヘルプ_一時ディレクトリ] =
                "临时目录，保留各 k 的组装结果（默认：{0}）",
            [メッセージID.ヘルプ_一時ディレクトリ削除] =
                "运行结束后删除临时目录（默认：保留）",
            [メッセージID.ヘルプ_スレッド数] =
                "读取 read 所用的工作线程数（默认：逻辑处理器数）",
            [メッセージID.ヘルプ_言語] =
                "消息语言（默认：{0}）",
            [メッセージID.ヘルプ_バージョン] =
                "显示版本并退出",
            [メッセージID.ヘルプ_ヘルプ] =
                "显示本说明",
            [メッセージID.引き継ぎの統合開始] =
                "正在合并从上一个 k 传递来的 {0:N0} 条序列",
            [メッセージID.引き継ぎの統合進捗] =
                "[Carry-over] 已合并 {1:N0} 条中的 {0:N0} 条",
            [メッセージID.引き継ぎの準備開始] =
                "正在准备传递给下一个 k 的序列",
            [メッセージID.引き継ぎの準備完了] =
                "[Carry-over] 已为下一个 k 准备 {0:N0} 条序列。",
            [メッセージID.SuperRead橋渡し開始] =
                "[SuperRead] 正在桥接 read 对",
            [メッセージID.SuperRead橋渡し進捗] =
                "[SuperRead] 已处理 {0:N0} 对，其中 {1:N0} 对完成桥接",
        };
    }
}
