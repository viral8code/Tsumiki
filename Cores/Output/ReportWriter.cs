using System.Globalization;
using System.Text;
using Tsumiki.Commons;
using Tsumiki.Cores.Evaluation;
using Tsumiki.Models.Evaluation;
using Tsumiki.Models.Polishing;
using Tsumiki.Models.Reporting;

namespace Tsumiki.Cores.Output
{
    /// <summary>
    /// 完全長の判定と、その根拠になった数値をファイルへ書き出す
    /// </summary>
    internal static class ReportWriter
    {
        #region 定数

        /// <summary>
        /// JSON オブジェクト開始
        /// </summary>
        private const string C_JSONオブジェクト開始 = "{";

        /// <summary>
        /// JSON 書式 schema version
        /// </summary>
        private const string C_JSON書式_schema_version = "  \"schema_version\": {0},";

        /// <summary>
        /// JSON 書式 tsumiki version
        /// </summary>
        private const string C_JSON書式_tsumiki_version = "  \"tsumiki_version\": {0},";

        /// <summary>
        /// JSON 書式 complete
        /// </summary>
        private const string C_JSON書式_complete = "  \"complete\": {0},";

        /// <summary>
        /// JSON 真
        /// </summary>
        private const string C_JSON真 = "true";

        /// <summary>
        /// JSON 偽
        /// </summary>
        private const string C_JSON偽 = "false";

        /// <summary>
        /// JSON 書式 quality level
        /// </summary>
        private const string C_JSON書式_quality_level = "  \"quality_level\": {0},";

        /// <summary>
        /// 品質保証レベル接頭辞
        /// </summary>
        private const string C_品質保証レベル接頭辞 = "Q";

        /// <summary>
        /// JSON 書式 reason codes
        /// </summary>
        private const string C_JSON書式_reason_codes = "  \"reason_codes\": [{0}],";

        /// <summary>
        /// 項目区切り
        /// </summary>
        private const string C_項目区切り = ", ";

        /// <summary>
        /// JSON 書式 k
        /// </summary>
        private const string C_JSON書式_k = "  \"k\": {0},";

        /// <summary>
        /// JSON 書式 sequences
        /// </summary>
        private const string C_JSON書式_sequences = "  \"sequences\": {0},";

        /// <summary>
        /// JSON 書式 total length
        /// </summary>
        private const string C_JSON書式_total_length = "  \"total_length\": {0},";

        /// <summary>
        /// JSON 書式 largest
        /// </summary>
        private const string C_JSON書式_largest = "  \"largest\": {0},";

        /// <summary>
        /// JSON 書式 n50
        /// </summary>
        private const string C_JSON書式_n50 = "  \"n50\": {0},";

        /// <summary>
        /// JSON 書式 l50
        /// </summary>
        private const string C_JSON書式_l50 = "  \"l50\": {0},";

        /// <summary>
        /// JSON 書式 gc percent
        /// </summary>
        private const string C_JSON書式_gc_percent = "  \"gc_percent\": {0},";

        /// <summary>
        /// 項目 scaffold stats
        /// </summary>
        private const string C_項目_scaffold_stats = "scaffold_stats";

        /// <summary>
        /// 項目 n split contig stats
        /// </summary>
        private const string C_項目_n_split_contig_stats = "n_split_contig_stats";

        /// <summary>
        /// JSON 書式 assembly settings
        /// </summary>
        private const string C_JSON書式_assembly_settings = "  \"assembly_settings\": {";

        /// <summary>
        /// JSON 書式 copy number baseline requested
        /// </summary>
        private const string C_JSON書式_copy_number_baseline_requested = "    \"copy_number_baseline_requested\": {0},";

        /// <summary>
        /// JSON 書式 copy number baseline actual
        /// </summary>
        private const string C_JSON書式_copy_number_baseline_actual = "    \"copy_number_baseline_actual\": {0},";

        /// <summary>
        /// JSON 書式 trim low coverage ends
        /// </summary>
        private const string C_JSON書式_trim_low_coverage_ends = "    \"trim_low_coverage_ends\": {0}";

        /// <summary>
        /// JSON 空値
        /// </summary>
        private const string C_JSON空値 = "null";

        /// <summary>
        /// JSON 内側オブジェクト終了と区切り
        /// </summary>
        private const string C_JSON内側オブジェクト終了と区切り = "  },";

        /// <summary>
        /// JSON 書式 circular replicons
        /// </summary>
        private const string C_JSON書式_circular_replicons = "  \"circular_replicons\": {0},";

        /// <summary>
        /// JSON 書式 unresolved gaps
        /// </summary>
        private const string C_JSON書式_unresolved_gaps = "  \"unresolved_gaps\": {0},";

        /// <summary>
        /// JSON 書式 checks
        /// </summary>
        private const string C_JSON書式_checks = "  \"checks\": [";

        /// <summary>
        /// 区切り
        /// </summary>
        private const string C_区切り = ",";

        /// <summary>
        /// JSON 行書式 name
        /// </summary>
        private const string C_JSON行書式_name = "    {{\"name\": {0}, \"result\": {1}, \"detail\": {2}}}{3}";

        /// <summary>
        /// JSON 内側配列終了と区切り
        /// </summary>
        private const string C_JSON内側配列終了と区切り = "  ],";

        /// <summary>
        /// JSON 書式 ambiguous junctions
        /// </summary>
        private const string C_JSON書式_ambiguous_junctions = "  \"ambiguous_junctions\": {0},";

        /// <summary>
        /// JSON オブジェクト終了
        /// </summary>
        private const string C_JSONオブジェクト終了 = "}";

        /// <summary>
        /// 見出し Tsumiki assembly report
        /// </summary>
        private const string C_見出し_Tsumiki_assembly_report = "# Tsumiki assembly report";

        /// <summary>
        /// 見出し Summary
        /// </summary>
        private const string C_見出し_Summary = "## Summary";

        /// <summary>
        /// 項目 item   value
        /// </summary>
        private const string C_項目_item___value = "| item | value |";

        /// <summary>
        /// 項目 V 書き出し Markdown レポート
        /// </summary>
        private const string C_項目_V_書き出し_Markdownレポート = "|---|---|";

        /// <summary>
        /// 項目 Tsumiki version
        /// </summary>
        private const string C_項目_Tsumiki_version = "Tsumiki version";

        /// <summary>
        /// 項目 adopted k
        /// </summary>
        private const string C_項目_adopted_k = "adopted k";

        /// <summary>
        /// 項目 complete
        /// </summary>
        private const string C_項目_complete = "complete";

        /// <summary>
        /// 項目 yes
        /// </summary>
        private const string C_項目_yes = "yes";

        /// <summary>
        /// 項目 no
        /// </summary>
        private const string C_項目_no = "no";

        /// <summary>
        /// 項目 quality level
        /// </summary>
        private const string C_項目_quality_level = "quality level";

        /// <summary>
        /// 項目 reasons not complete
        /// </summary>
        private const string C_項目_reasons_not_complete = "reasons not complete";

        /// <summary>
        /// 未指定表示
        /// </summary>
        private const string C_未指定表示 = "-";

        /// <summary>
        /// 項目 circular replicons
        /// </summary>
        private const string C_項目_circular_replicons = "circular replicons";

        /// <summary>
        /// 項目 unresolved gaps
        /// </summary>
        private const string C_項目_unresolved_gaps = "unresolved gaps";

        /// <summary>
        /// 項目 ambiguous junctions
        /// </summary>
        private const string C_項目_ambiguous_junctions = "ambiguous junctions";

        /// <summary>
        /// 項目 copy number baseline  requested   actual
        /// </summary>
        private const string C_項目_copy_number_baseline__requested___actual = "copy-number baseline (requested / actual)";

        /// <summary>
        /// 見出し Sequence statistics
        /// </summary>
        private const string C_見出し_Sequence_statistics = "## Sequence statistics";

        /// <summary>
        /// 見出し Completeness checks
        /// </summary>
        private const string C_見出し_Completeness_checks = "## Completeness checks";

        /// <summary>
        /// 項目 check   result   detail
        /// </summary>
        private const string C_項目_check___result___detail = "| check | result | detail |";

        /// <summary>
        /// 検査項目の列書式
        /// </summary>
        private const string C_検査項目の列書式 = "|---|---|---|";

        /// <summary>
        /// 見出し Self check against trusted k mers
        /// </summary>
        private const string C_見出し_Self_check_against_trusted_k_mers = "## Self-check against trusted k-mers";

        /// <summary>
        /// 項目 trusted k mers
        /// </summary>
        private const string C_項目_trusted_k_mers = "trusted k-mers";

        /// <summary>
        /// 項目 N0
        /// </summary>
        private const string C_項目_N0 = "N0";

        /// <summary>
        /// 項目 missing k mers
        /// </summary>
        private const string C_項目_missing_k_mers = "missing k-mers";

        /// <summary>
        /// 項目 excess copies
        /// </summary>
        private const string C_項目_excess_copies = "excess copies";

        /// <summary>
        /// 項目 not measured
        /// </summary>
        private const string C_項目_not_measured = "not measured";

        /// <summary>
        /// 見出し Fixed anchor evaluation
        /// </summary>
        private const string C_見出し_Fixed_anchor_evaluation = "## Fixed-anchor evaluation";

        /// <summary>
        /// 項目 completeness
        /// </summary>
        private const string C_項目_completeness = "completeness";

        /// <summary>
        /// 項目 accuracy
        /// </summary>
        private const string C_項目_accuracy = "accuracy";

        /// <summary>
        /// 項目 NG50
        /// </summary>
        private const string C_項目_NG50 = "NG50";

        /// <summary>
        /// 項目 sequences
        /// </summary>
        private const string C_項目_sequences = "sequences";

        /// <summary>
        /// 見出し Polishing
        /// </summary>
        private const string C_見出し_Polishing = "## Polishing";

        /// <summary>
        /// 項目 mapped reads
        /// </summary>
        private const string C_項目_mapped_reads = "mapped reads";

        /// <summary>
        /// 項目 corrected bases
        /// </summary>
        private const string C_項目_corrected_bases = "corrected bases";

        /// <summary>
        /// 項目 median depth
        /// </summary>
        private const string C_項目_median_depth = "median depth";

        /// <summary>
        /// 項目 low depth positions
        /// </summary>
        private const string C_項目_low_depth_positions = "low-depth positions";

        /// <summary>
        /// 見出し Read support
        /// </summary>
        private const string C_見出し_Read_support = "## Read support";

        /// <summary>
        /// 項目 r mer length
        /// </summary>
        private const string C_項目_r_mer_length = "r-mer length";

        /// <summary>
        /// 項目 unsupported positions
        /// </summary>
        private const string C_項目_unsupported_positions = "unsupported positions";

        /// <summary>
        /// 項目 unsupported intervals
        /// </summary>
        private const string C_項目_unsupported_intervals = "unsupported intervals";

        /// <summary>
        /// 見出し Circular closure
        /// </summary>
        private const string C_見出し_Circular_closure = "## Circular closure";

        /// <summary>
        /// 項目 sequence   length   span ired   supported
        /// </summary>
        private const string C_項目_sequence___length___span_ired___supported = "| sequence | length | spanning reads | required | supported |";

        /// <summary>
        /// 支持検査の列書式
        /// </summary>
        private const string C_支持検査の列書式 = "|---|---:|---:|---:|---|";

        /// <summary>
        /// 項目 no circular sequences
        /// </summary>
        private const string C_項目_no_circular_sequences = "no circular sequences";

        /// <summary>
        /// 見出し Stage timings
        /// </summary>
        private const string C_見出し_Stage_timings = "## Stage timings";

        /// <summary>
        /// 項目 stage   elapsed s   CPU s   peak working set MB
        /// </summary>
        private const string C_項目_stage___elapsed_s___CPU_s___peak_working_set_MB = "| stage | elapsed s | CPU s | peak working set MB |";

        /// <summary>
        /// フェーズ計測の列書式
        /// </summary>
        private const string C_フェーズ計測の列書式 = "|---|---:|---:|---:|";

        /// <summary>
        /// 項目 sequence
        /// </summary>
        private const string C_項目_sequence = "sequence";

        /// <summary>
        /// 項目 start
        /// </summary>
        private const string C_項目_start = "start";

        /// <summary>
        /// 項目 end
        /// </summary>
        private const string C_項目_end = "end";

        /// <summary>
        /// 項目 length
        /// </summary>
        private const string C_項目_length = "length";

        /// <summary>
        /// 項目 r
        /// </summary>
        private const string C_項目_r = "r";

        /// <summary>
        /// 項目 k
        /// </summary>
        private const string C_項目_k = "k";

        /// <summary>
        /// 項目 type
        /// </summary>
        private const string C_項目_type = "type";

        /// <summary>
        /// 項目 location
        /// </summary>
        private const string C_項目_location = "location";

        /// <summary>
        /// 項目 stable id
        /// </summary>
        private const string C_項目_stable_id = "stable_id";

        /// <summary>
        /// 項目 top support
        /// </summary>
        private const string C_項目_top_support = "top_support";

        /// <summary>
        /// 項目 second support
        /// </summary>
        private const string C_項目_second_support = "second_support";

        /// <summary>
        /// 項目 margin
        /// </summary>
        private const string C_項目_margin = "margin";

        /// <summary>
        /// 項目 raw support
        /// </summary>
        private const string C_項目_raw_support = "raw_support";

        /// <summary>
        /// 項目 confidence
        /// </summary>
        private const string C_項目_confidence = "confidence";

        /// <summary>
        /// 項目 repeat length
        /// </summary>
        private const string C_項目_repeat_length = "repeat_length";

        /// <summary>
        /// 項目 gap
        /// </summary>
        private const string C_項目_gap = "gap";

        /// <summary>
        /// 項目 copy number
        /// </summary>
        private const string C_項目_copy_number = "copy_number";

        /// <summary>
        /// 項目 spanning reads
        /// </summary>
        private const string C_項目_spanning_reads = "spanning_reads";

        /// <summary>
        /// 項目 spanning pairs
        /// </summary>
        private const string C_項目_spanning_pairs = "spanning_pairs";

        /// <summary>
        /// 項目 expected pairs
        /// </summary>
        private const string C_項目_expected_pairs = "expected_pairs";

        /// <summary>
        /// 項目 pair ratio
        /// </summary>
        private const string C_項目_pair_ratio = "pair_ratio";

        /// <summary>
        /// 項目 discordant anchors
        /// </summary>
        private const string C_項目_discordant_anchors = "discordant_anchors";

        /// <summary>
        /// 項目 discordant fraction
        /// </summary>
        private const string C_項目_discordant_fraction = "discordant_fraction";

        /// <summary>
        /// 項目 clipped
        /// </summary>
        private const string C_項目_clipped = "clipped";

        /// <summary>
        /// 項目 left depth ratio
        /// </summary>
        private const string C_項目_left_depth_ratio = "left_depth_ratio";

        /// <summary>
        /// 項目 right depth ratio
        /// </summary>
        private const string C_項目_right_depth_ratio = "right_depth_ratio";

        /// <summary>
        /// 項目 class
        /// </summary>
        private const string C_項目_class = "class";

        /// <summary>
        /// 項目 misjoin probability
        /// </summary>
        private const string C_項目_misjoin_probability = "misjoin_probability";

        /// <summary>
        /// ファイル名 0
        /// </summary>
        private const string C_ファイル名_0 = "0.######";

        /// <summary>
        /// 項目 span
        /// </summary>
        private const string C_項目_span = "span";

        /// <summary>
        /// 項目 nospan
        /// </summary>
        private const string C_項目_nospan = "nospan";

        /// <summary>
        /// 表区切り
        /// </summary>
        private const string C_表区切り = "|";

        /// <summary>
        /// 表区切りエスケープ
        /// </summary>
        private const string C_表区切りエスケープ = "\\|";

        /// <summary>
        /// 書式
        /// </summary>
        private const string C_書式 = "\r";

        /// <summary>
        /// 空白
        /// </summary>
        private const string C_空白 = " ";

        /// <summary>
        /// 改行
        /// </summary>
        private const string C_改行 = "\n";

        /// <summary>
        /// 項目 0
        /// </summary>
        private const string C_項目_0 = "  \"{0}\": {{";

        /// <summary>
        /// JSON 書式 minimum sequence length
        /// </summary>
        private const string C_JSON書式_minimum_sequence_length = "    \"minimum_sequence_length\": {0},";

        /// <summary>
        /// JSON 書式 sequences 区切りあり
        /// </summary>
        private const string C_JSON書式_sequences_区切りあり = "    \"sequences\": {0},";

        /// <summary>
        /// JSON 書式 total length 区切りあり
        /// </summary>
        private const string C_JSON書式_total_length_区切りあり = "    \"total_length\": {0},";

        /// <summary>
        /// JSON 書式 largest 区切りあり
        /// </summary>
        private const string C_JSON書式_largest_区切りあり = "    \"largest\": {0},";

        /// <summary>
        /// JSON 書式 smallest
        /// </summary>
        private const string C_JSON書式_smallest = "    \"smallest\": {0},";

        /// <summary>
        /// JSON 書式 n50 区切りあり
        /// </summary>
        private const string C_JSON書式_n50_区切りあり = "    \"n50\": {0},";

        /// <summary>
        /// JSON 書式 l50 区切りあり
        /// </summary>
        private const string C_JSON書式_l50_区切りあり = "    \"l50\": {0},";

        /// <summary>
        /// JSON 書式 gc percent 区切りなし
        /// </summary>
        private const string C_JSON書式_gc_percent_区切りなし = "    \"gc_percent\": {0}";

        /// <summary>
        /// 項目 no support
        /// </summary>
        private const string C_項目_no_support = "no-support";

        /// <summary>
        /// 項目 insufficient margin
        /// </summary>
        private const string C_項目_insufficient_margin = "insufficient-margin";

        /// <summary>
        /// 項目 multiple paths
        /// </summary>
        private const string C_項目_multiple_paths = "multiple-paths";

        /// <summary>
        /// 項目 unreachable
        /// </summary>
        private const string C_項目_unreachable = "unreachable";

        /// <summary>
        /// 項目 anchor missing
        /// </summary>
        private const string C_項目_anchor_missing = "anchor-missing";

        /// <summary>
        /// 項目 no local reads
        /// </summary>
        private const string C_項目_no_local_reads = "no-local-reads";

        /// <summary>
        /// 項目 search cutoff
        /// </summary>
        private const string C_項目_search_cutoff = "search-cutoff";

        /// <summary>
        /// 項目 inside repeat
        /// </summary>
        private const string C_項目_inside_repeat = "inside-repeat";

        /// <summary>
        /// JSON 書式 phase timings
        /// </summary>
        private const string C_JSON書式_phase_timings = "  \"phase_timings\": [";

        /// <summary>
        /// JSON 行書式 stage
        /// </summary>
        private const string C_JSON行書式_stage = "    {{\"stage\": {0}, \"elapsed_s\": {1}, \"cpu_s\": {2}, \"allocated_mb\": {3}, \"working_set_mb\": {4}, \"peak_working_set_mb\": {5}, \"gen2_gc\": {6}}}{7}";

        /// <summary>
        /// 項目 V 追加 フェーズ計測
        /// </summary>
        private const string C_項目_V_追加_フェーズ計測 = "  ]";

        /// <summary>
        /// JSON 書式 self check
        /// </summary>
        private const string C_JSON書式_self_check = "  \"self_check\": null,";

        /// <summary>
        /// JSON 書式 self check 区切りなし
        /// </summary>
        private const string C_JSON書式_self_check_区切りなし = "  \"self_check\": {";

        /// <summary>
        /// JSON 書式 trusted kmers
        /// </summary>
        private const string C_JSON書式_trusted_kmers = "    \"trusted_kmers\": {0},";

        /// <summary>
        /// JSON 書式 missing kmers
        /// </summary>
        private const string C_JSON書式_missing_kmers = "    \"missing_kmers\": {0},";

        /// <summary>
        /// JSON 書式 missing percent
        /// </summary>
        private const string C_JSON書式_missing_percent = "    \"missing_percent\": {0},";

        /// <summary>
        /// JSON 書式 excess percent
        /// </summary>
        private const string C_JSON書式_excess_percent = "    \"excess_percent\": {0}";

        /// <summary>
        /// JSON 書式 anchor evaluation
        /// </summary>
        private const string C_JSON書式_anchor_evaluation = "  \"anchor_evaluation\": null,";

        /// <summary>
        /// JSON 書式 anchor evaluation 区切りなし
        /// </summary>
        private const string C_JSON書式_anchor_evaluation_区切りなし = "  \"anchor_evaluation\": {";

        /// <summary>
        /// JSON 書式 expected copies
        /// </summary>
        private const string C_JSON書式_expected_copies = "    \"expected_copies\": {0},";

        /// <summary>
        /// JSON 書式 missing copies
        /// </summary>
        private const string C_JSON書式_missing_copies = "    \"missing_copies\": {0},";

        /// <summary>
        /// JSON 書式 excess copies
        /// </summary>
        private const string C_JSON書式_excess_copies = "    \"excess_copies\": {0},";

        /// <summary>
        /// JSON 書式 completeness
        /// </summary>
        private const string C_JSON書式_completeness = "    \"completeness\": {0},";

        /// <summary>
        /// JSON 書式 accuracy
        /// </summary>
        private const string C_JSON書式_accuracy = "    \"accuracy\": {0},";

        /// <summary>
        /// JSON 書式 ng50
        /// </summary>
        private const string C_JSON書式_ng50 = "    \"ng50\": {0},";

        /// <summary>
        /// JSON 書式 circular replicons 区切りあり
        /// </summary>
        private const string C_JSON書式_circular_replicons_区切りあり = "    \"circular_replicons\": {0},";

        /// <summary>
        /// JSON 書式 circularized fraction
        /// </summary>
        private const string C_JSON書式_circularized_fraction = "    \"circularized_fraction\": {0}";

        /// <summary>
        /// JSON 書式 polish
        /// </summary>
        private const string C_JSON書式_polish = "  \"polish\": null,";

        /// <summary>
        /// JSON 書式 polish 区切りなし
        /// </summary>
        private const string C_JSON書式_polish_区切りなし = "  \"polish\": {";

        /// <summary>
        /// JSON 書式 mapped reads
        /// </summary>
        private const string C_JSON書式_mapped_reads = "    \"mapped_reads\": {0},";

        /// <summary>
        /// JSON 書式 rejected reads
        /// </summary>
        private const string C_JSON書式_rejected_reads = "    \"rejected_reads\": {0},";

        /// <summary>
        /// JSON 書式 corrected bases
        /// </summary>
        private const string C_JSON書式_corrected_bases = "    \"corrected_bases\": {0},";

        /// <summary>
        /// JSON 書式 median depth
        /// </summary>
        private const string C_JSON書式_median_depth = "    \"median_depth\": {0},";

        /// <summary>
        /// JSON 書式 low depth fraction
        /// </summary>
        private const string C_JSON書式_low_depth_fraction = "    \"low_depth_fraction\": {0}";

        /// <summary>
        /// JSON 書式 circular closure
        /// </summary>
        private const string C_JSON書式_circular_closure = "  \"circular_closure\": null,";

        /// <summary>
        /// JSON 書式 circular closure 区切りなし
        /// </summary>
        private const string C_JSON書式_circular_closure_区切りなし = "  \"circular_closure\": [";

        /// <summary>
        /// JSON 行書式 id
        /// </summary>
        private const string C_JSON行書式_id = "    {{\"id\": {0}, \"length\": {1}, \"spanning_reads\": {2}, \"required\": {3}, \"supported\": {4}}}{5}";

        /// <summary>
        /// ファイル名 0 Get 数値
        /// </summary>
        private const string C_ファイル名_0_Get_数値 = "0.####";

        /// <summary>
        /// 引用符
        /// </summary>
        private const string C_引用符 = "\"";

        /// <summary>
        /// JSON 引用符エスケープ
        /// </summary>
        private const string C_JSON引用符エスケープ = "\\\"";

        /// <summary>
        /// JSON 逆斜線エスケープ
        /// </summary>
        private const string C_JSON逆斜線エスケープ = "\\\\";

        /// <summary>
        /// JSON 改行エスケープ
        /// </summary>
        private const string C_JSON改行エスケープ = "\\n";

        /// <summary>
        /// JSON 復帰文字エスケープ
        /// </summary>
        private const string C_JSON復帰文字エスケープ = "\\r";

        /// <summary>
        /// JSON タブエスケープ
        /// </summary>
        private const string C_JSONタブエスケープ = "\\t";

        #endregion

        #region 公開メソッド

        /// <summary>
        /// 完全性レポートを JSON で書き出す
        /// </summary>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_k長">採用した k の長さ</param>
        /// <param name="p_統計">アセンブリ統計</param>
        /// <param name="p_未解決ギャップ数">埋まらなかったギャップの数</param>
        /// <param name="p_環状本数">環状に閉じた配列の本数</param>
        /// <param name="p_判定">完全長の判定結果</param>
        /// <param name="p_整合性">自己検査の結果</param>
        /// <param name="p_閉鎖検証">環状閉鎖の検証結果</param>
        /// <param name="p_ポリッシュ">ポリッシュの結果</param>
        /// <param name="p_曖昧箇所">決めきれなかった箇所</param>
        /// <param name="p_N分割統計"></param>
        /// <param name="p_統計の最小長"></param>
        /// <param name="p_要求コピー数基準"></param>
        /// <param name="p_実際のコピー数基準"></param>
        /// <param name="p_Is低カバレッジ端トリミング"></param>
        /// <param name="p_scaffold比較統計"></param>
        /// <param name="p_固定アンカー評価">採用した k・コピー数基準に依らない固定アンカーでの独立評価</param>
        /// <param name="p_フェーズ計測">工程ごとの資源使用量</param>
        public static void V_書き出し_レポート(string p_出力パス, int p_k長, アセンブリ統計 p_統計, int p_未解決ギャップ数, int p_環状本数, 完全性判定結果 p_判定, 整合性検査結果? p_整合性, IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証, ポリッシュ統計? p_ポリッシュ, IReadOnlyList<曖昧箇所> p_曖昧箇所, アセンブリ統計? p_N分割統計 = null, int p_統計の最小長 = 500, string? p_要求コピー数基準 = null, string? p_実際のコピー数基準 = null, bool? p_Is低カバレッジ端トリミング = null, アセンブリ統計? p_scaffold比較統計 = null, アセンブリ評価? p_固定アンカー評価 = null, IReadOnlyList<フェーズ計測>? p_フェーズ計測 = null)
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine(C_JSONオブジェクト開始);
            V_追加(l_文, C_JSON書式_schema_version, 1);
            V_追加(l_文, C_JSON書式_tsumiki_version, Get_文字列(Consts.バージョン));
            V_追加(l_文, C_JSON書式_complete, p_判定.A_Is完全長 ? C_JSON真 : C_JSON偽);
            V_追加(l_文, C_JSON書式_quality_level, Get_文字列(C_品質保証レベル接頭辞 + (int)p_判定.A_品質保証レベル));
            V_追加(l_文, C_JSON書式_reason_codes, string.Join(C_項目区切り, p_判定.A_未達理由.Select(x => Get_文字列(CompletenessValidator.Get_理由コード(x)))));
            V_追加(l_文, C_JSON書式_k, p_k長);
            V_追加(l_文, C_JSON書式_sequences, p_統計.A_配列数);
            V_追加(l_文, C_JSON書式_total_length, p_統計.A_総延長);
            V_追加(l_文, C_JSON書式_largest, p_統計.A_最大長);
            V_追加(l_文, C_JSON書式_n50, p_統計.A_N50);
            V_追加(l_文, C_JSON書式_l50, p_統計.A_L50);
            V_追加(l_文, C_JSON書式_gc_percent, Get_数値(p_統計.A_GC率));
            V_追加_統計(l_文, C_項目_scaffold_stats, p_scaffold比較統計 ?? p_統計, p_統計の最小長);
            V_追加_統計(l_文, C_項目_n_split_contig_stats, p_N分割統計 ?? p_統計, p_統計の最小長);
            _ = l_文.AppendLine(C_JSON書式_assembly_settings);
            V_追加(l_文, C_JSON書式_copy_number_baseline_requested, Get_文字列(p_要求コピー数基準));
            V_追加(l_文, C_JSON書式_copy_number_baseline_actual, Get_文字列(p_実際のコピー数基準));
            V_追加(l_文, C_JSON書式_trim_low_coverage_ends, p_Is低カバレッジ端トリミング is { } l_trim ? (l_trim ? C_JSON真 : C_JSON偽) : C_JSON空値);
            _ = l_文.AppendLine(C_JSON内側オブジェクト終了と区切り);
            V_追加(l_文, C_JSON書式_circular_replicons, p_環状本数);
            V_追加(l_文, C_JSON書式_unresolved_gaps, p_未解決ギャップ数);
            _ = l_文.AppendLine(C_JSON書式_checks);
            for (var i = 0; i < p_判定.A_検査項目.Count; i++)
            {
                var l_項目 = p_判定.A_検査項目[i];
                var l_末尾 = i == p_判定.A_検査項目.Count - 1 ? string.Empty : C_区切り;
                V_追加(l_文, C_JSON行書式_name, Get_文字列(l_項目.A_キー), Get_文字列(CompletenessValidator.Get_判定コード(l_項目.A_判定)), Get_文字列(l_項目.A_内訳), l_末尾);
            }

            _ = l_文.AppendLine(C_JSON内側配列終了と区切り);
            V_追加_自己検査(l_文, p_整合性);
            V_追加_アンカー評価(l_文, p_固定アンカー評価);
            V_追加_ポリッシュ(l_文, p_ポリッシュ);
            V_追加_閉鎖検証(l_文, p_閉鎖検証);
            V_追加(l_文, C_JSON書式_ambiguous_junctions, p_曖昧箇所.Count);
            V_追加_フェーズ計測(l_文, p_フェーズ計測 ?? []);
            _ = l_文.AppendLine(C_JSONオブジェクト終了);
            File.WriteAllText(p_出力パス, l_文.ToString());
        }

        /// <summary>
        /// 最終成果物の統計と検査・計測の結果を、人が読む Markdown で書き出す
        /// </summary>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_k長">採用した k の長さ</param>
        /// <param name="p_統計表"><see cref="AssemblyStatsReporter.Get_統計表"/> の行</param>
        /// <param name="p_判定">完全長の判定結果</param>
        /// <param name="p_未解決ギャップ数">埋まらなかったギャップの数</param>
        /// <param name="p_環状本数">環状に閉じた配列の本数</param>
        /// <param name="p_曖昧箇所数">決めきれなかった箇所の数</param>
        /// <param name="p_要求コピー数基準">要求したコピー数基準</param>
        /// <param name="p_実際のコピー数基準">実際に使ったコピー数基準</param>
        /// <param name="p_整合性">自己検査の結果</param>
        /// <param name="p_固定アンカー評価">固定アンカーでの独立評価</param>
        /// <param name="p_ポリッシュ">ポリッシュの結果</param>
        /// <param name="p_支持検査">リード支持の検査結果</param>
        /// <param name="p_閉鎖検証">環状閉鎖の検証結果</param>
        /// <param name="p_フェーズ計測">工程ごとの資源使用量</param>
        public static void V_書き出し_Markdownレポート(string p_出力パス, int p_k長, IReadOnlyList<string> p_統計表, 完全性判定結果 p_判定, int p_未解決ギャップ数, int p_環状本数, int p_曖昧箇所数, string? p_要求コピー数基準, string? p_実際のコピー数基準, 整合性検査結果? p_整合性, アセンブリ評価? p_固定アンカー評価, ポリッシュ統計? p_ポリッシュ, 支持検査結果? p_支持検査, IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証, IReadOnlyList<フェーズ計測> p_フェーズ計測)
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine(C_見出し_Tsumiki_assembly_report).AppendLine();
            _ = l_文.AppendLine(C_見出し_Summary).AppendLine().AppendLine(C_項目_item___value).AppendLine(C_項目_V_書き出し_Markdownレポート);
            V_追加_表の行(l_文, C_項目_Tsumiki_version, Consts.バージョン);
            V_追加_表の行(l_文, C_項目_adopted_k, p_k長);
            V_追加_表の行(l_文, C_項目_complete, p_判定.A_Is完全長 ? C_項目_yes : C_項目_no);
            V_追加_表の行(l_文, C_項目_quality_level, C_品質保証レベル接頭辞 + (int)p_判定.A_品質保証レベル);
            V_追加_表の行(l_文, C_項目_reasons_not_complete, p_判定.A_未達理由.Count == 0 ? C_未指定表示 : string.Join(C_項目区切り, p_判定.A_未達理由.Select(CompletenessValidator.Get_理由コード)));
            V_追加_表の行(l_文, C_項目_circular_replicons, p_環状本数);
            V_追加_表の行(l_文, C_項目_unresolved_gaps, p_未解決ギャップ数);
            V_追加_表の行(l_文, C_項目_ambiguous_junctions, p_曖昧箇所数);
            V_追加_表の行(l_文, C_項目_copy_number_baseline__requested___actual, $"{p_要求コピー数基準 ?? C_未指定表示} / {p_実際のコピー数基準 ?? C_未指定表示}");
            _ = l_文.AppendLine();
            _ = l_文.AppendLine(C_見出し_Sequence_statistics).AppendLine();
            foreach (var l_行 in p_統計表)
            {
                _ = l_文.AppendLine(l_行);
            }

            _ = l_文.AppendLine();
            _ = l_文.AppendLine(C_見出し_Completeness_checks).AppendLine().AppendLine(C_項目_check___result___detail).AppendLine(C_検査項目の列書式);
            foreach (var l_項目 in p_判定.A_検査項目)
            {
                _ = l_文.AppendLine($"| {Get_表の値(l_項目.A_キー)} | {CompletenessValidator.Get_判定コード(l_項目.A_判定)} | {Get_表の値(l_項目.A_内訳)} |");
            }

            _ = l_文.AppendLine();
            _ = l_文.AppendLine(C_見出し_Self_check_against_trusted_k_mers).AppendLine();
            if (p_整合性 is { } l_整合性)
            {
                _ = l_文.AppendLine(C_項目_item___value).AppendLine(C_項目_V_書き出し_Markdownレポート);
                V_追加_表の行(l_文, C_項目_trusted_k_mers, l_整合性.A_信頼kmer数.ToString(C_項目_N0, CultureInfo.InvariantCulture));
                V_追加_表の行(l_文, C_項目_missing_k_mers, string.Create(CultureInfo.InvariantCulture, $"{l_整合性.A_取りこぼし数:N0} ({l_整合性.A_取りこぼし率:F3}%)"));
                V_追加_表の行(l_文, C_項目_excess_copies, string.Create(CultureInfo.InvariantCulture, $"{l_整合性.A_出しすぎ率:F3}%"));
            }
            else
            {
                _ = l_文.AppendLine(C_項目_not_measured);
            }

            _ = l_文.AppendLine();
            _ = l_文.AppendLine(C_見出し_Fixed_anchor_evaluation).AppendLine();
            if (p_固定アンカー評価 is { } l_評価)
            {
                _ = l_文.AppendLine(C_項目_item___value).AppendLine(C_項目_V_書き出し_Markdownレポート);
                V_追加_表の行(l_文, C_項目_completeness, string.Create(CultureInfo.InvariantCulture, $"{l_評価.A_完全性 * 100D:F2}%"));
                V_追加_表の行(l_文, C_項目_accuracy, string.Create(CultureInfo.InvariantCulture, $"{l_評価.A_正確性 * 100D:F2}%"));
                V_追加_表の行(l_文, C_項目_NG50, l_評価.A_NG50.ToString(C_項目_N0, CultureInfo.InvariantCulture));
                V_追加_表の行(l_文, C_項目_sequences, l_評価.A_本数);
                V_追加_表の行(l_文, C_項目_circular_replicons, string.Create(CultureInfo.InvariantCulture, $"{l_評価.A_環状本数} ({l_評価.A_環状化率 * 100D:F1}% of genome)"));
            }
            else
            {
                _ = l_文.AppendLine(C_項目_not_measured);
            }

            _ = l_文.AppendLine();
            _ = l_文.AppendLine(C_見出し_Polishing).AppendLine();
            if (p_ポリッシュ is { } l_ポリッシュ)
            {
                _ = l_文.AppendLine(C_項目_item___value).AppendLine(C_項目_V_書き出し_Markdownレポート);
                V_追加_表の行(l_文, C_項目_mapped_reads, string.Create(CultureInfo.InvariantCulture, $"{l_ポリッシュ.A_マップされたリード数:N0} (rejected {l_ポリッシュ.A_棄却されたリード数:N0})"));
                V_追加_表の行(l_文, C_項目_corrected_bases, string.Create(CultureInfo.InvariantCulture, $"{l_ポリッシュ.A_訂正した塩基数:N0} / {l_ポリッシュ.A_総延長:N0}"));
                V_追加_表の行(l_文, C_項目_median_depth, string.Create(CultureInfo.InvariantCulture, $"{l_ポリッシュ.A_深度の中央値:F1}x"));
                V_追加_表の行(l_文, C_項目_low_depth_positions, string.Create(CultureInfo.InvariantCulture, $"{l_ポリッシュ.A_深度不足の位置数:N0} / {l_ポリッシュ.A_評価できた位置数:N0} ({l_ポリッシュ.A_深度不足率 * 100D:F2}%)"));
            }
            else
            {
                _ = l_文.AppendLine(C_項目_not_measured);
            }

            _ = l_文.AppendLine();
            _ = l_文.AppendLine(C_見出し_Read_support).AppendLine();
            if (p_支持検査 is { } l_支持検査)
            {
                _ = l_文.AppendLine(C_項目_item___value).AppendLine(C_項目_V_書き出し_Markdownレポート);
                V_追加_表の行(l_文, C_項目_r_mer_length, l_支持検査.A_r長);
                V_追加_表の行(l_文, C_項目_unsupported_positions, string.Create(CultureInfo.InvariantCulture, $"{l_支持検査.A_支持のない位置数:N0} / {l_支持検査.A_調べた位置数:N0} ({l_支持検査.A_支持のない率:F3}%)"));
                V_追加_表の行(l_文, C_項目_unsupported_intervals, l_支持検査.A_区間.Count.ToString(C_項目_N0, CultureInfo.InvariantCulture));
            }
            else
            {
                _ = l_文.AppendLine(C_項目_not_measured);
            }

            _ = l_文.AppendLine();
            _ = l_文.AppendLine(C_見出し_Circular_closure).AppendLine();
            if (p_閉鎖検証 is { Count: > 0 } l_閉鎖検証)
            {
                _ = l_文.AppendLine(C_項目_sequence___length___span_ired___supported).AppendLine(C_支持検査の列書式);
                foreach (var l_検証 in l_閉鎖検証)
                {
                    _ = l_文.AppendLine(string.Create(CultureInfo.InvariantCulture, $"| {Get_表の値(l_検証.A_配列ID)} | {l_検証.A_長さ:N0} | {l_検証.A_跨いだリード数:N0} | {l_検証.A_必要本数:N0} | {(l_検証.A_Has支持 ? C_項目_yes : C_項目_no)} |"));
                }
            }
            else
            {
                _ = l_文.AppendLine(p_閉鎖検証 is null ? C_項目_not_measured : C_項目_no_circular_sequences);
            }

            _ = l_文.AppendLine();
            _ = l_文.AppendLine(C_見出し_Stage_timings).AppendLine().AppendLine(C_項目_stage___elapsed_s___CPU_s___peak_working_set_MB).AppendLine(C_フェーズ計測の列書式);
            foreach (var l_計測 in p_フェーズ計測)
            {
                _ = l_文.AppendLine(string.Create(CultureInfo.InvariantCulture, $"| {Get_表の値(l_計測.A_工程)} | {l_計測.A_経過秒:F1} | {l_計測.A_CPU秒:F1} | {l_計測.A_ピークワーキングセットMB:N0} |"));
            }

            File.WriteAllText(p_出力パス, l_文.ToString());
        }

        /// <summary>
        /// リードに裏付けの無い区間を TSV で書き出す
        /// </summary>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_区間">支持のない区間</param>
        /// <param name="p_r長">支持を問うた r-mer の長さ</param>
        public static void V_書き出し_未支持箇所(string p_出力パス, IReadOnlyList<支持のない区間> p_区間, int p_r長)
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine(string.Join('\t', C_項目_sequence, C_項目_start, C_項目_end, C_項目_length, C_項目_r));
            foreach (var l_区間 in p_区間)
            {
                _ = l_文.AppendLine(string.Join('\t', l_区間.A_配列ID, l_区間.A_開始, l_区間.A_終了, l_区間.A_長さ, p_r長));
            }

            File.WriteAllText(p_出力パス, l_文.ToString());
        }

        /// <summary>
        /// 決めきれなかった箇所を TSV で書き出す
        /// </summary>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_曖昧箇所">決めきれなかった箇所</param>
        public static void V_書き出し_曖昧箇所(string p_出力パス, IReadOnlyList<曖昧箇所> p_曖昧箇所)
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine(string.Join('\t', C_項目_k, C_項目_type, C_項目_location, C_項目_stable_id, C_項目_top_support, C_項目_second_support, C_項目_margin, C_項目_raw_support, C_項目_confidence));
            foreach (var l_箇所 in p_曖昧箇所)
            {
                _ = l_文.AppendLine(string.Join('\t', l_箇所.A_k長, Get_種別コード(l_箇所.A_種別), l_箇所.A_場所, l_箇所.A_安定ID, Get_数値(l_箇所.A_首位の支持), Get_数値(l_箇所.A_次点の支持), Get_数値(l_箇所.A_余裕), l_箇所.A_首位の生支持数, Get_数値(l_箇所.A_確信度)));
            }

            File.WriteAllText(p_出力パス, l_文.ToString());
        }

        /// <summary>
        /// 継ぎ目ごとの証拠と、繋ぐ相手を誤っている確率を TSV で書き出す (位置は 0 始まり・終わりを含まない)
        /// </summary>
        /// <param name="p_出力パス">書き出し先</param>
        /// <param name="p_評価群">継ぎ目の評価</param>
        public static void V_書き出し_継ぎ目(string p_出力パス, IReadOnlyList<継ぎ目の評価> p_評価群)
        {
            var l_文 = new StringBuilder();
            _ = l_文.AppendLine(string.Join('\t', C_項目_sequence, C_項目_start, C_項目_end, C_項目_repeat_length, C_項目_gap, C_項目_copy_number, C_項目_spanning_reads, C_項目_spanning_pairs, C_項目_expected_pairs, C_項目_pair_ratio, C_項目_discordant_anchors, C_項目_discordant_fraction, C_項目_clipped, C_項目_left_depth_ratio, C_項目_right_depth_ratio, C_項目_class, C_項目_misjoin_probability));
            foreach (var (l_候補, l_組, l_確率) in p_評価群)
            {
                _ = l_文.AppendLine(string.Join('\t', l_候補.A_配列名, l_候補.A_開始, l_候補.A_終了, l_候補.A_反復長, l_候補.A_Isギャップ ? 1 : 0, Get_数値(l_候補.A_コピー数), l_候補.A_跨ぐ読み, l_候補.A_跨ぐ組, Get_数値(l_候補.A_期待の組), Get_数値(l_候補.A_組の比), l_候補.A_外れ錨, Get_数値(l_候補.A_外れ割合), Get_数値(l_候補.A_切れ端), Get_数値(l_候補.A_左の深さ比), Get_数値(l_候補.A_右の深さ比), Get_組コード(l_組), l_確率.ToString(C_ファイル名_0, CultureInfo.InvariantCulture)));
            }

            File.WriteAllText(p_出力パス, l_文.ToString());
        }

        #endregion

        #region 内部メソッド

        /// <summary>
        /// 継ぎ目の組を TSV に書く短い名前
        /// </summary>
        /// <param name="p_組"></param>
        /// <returns></returns>
        private static string Get_組コード(継ぎ目の組 p_組)
        {
            return p_組 switch
            {
                継ぎ目の組.ギャップ => C_項目_gap,
                継ぎ目の組.組で跨げる => C_項目_span,
                _ => C_項目_nospan,
            };
        }

        /// <summary>
        /// Markdown の 2 列の表に 1 行足す
        /// </summary>
        /// <param name="p_文"></param>
        /// <param name="p_項目"></param>
        /// <param name="p_値"></param>
        private static void V_追加_表の行(StringBuilder p_文, string p_項目, object p_値)
        {
            _ = p_文.AppendLine($"| {p_項目} | {Get_表の値(Convert.ToString(p_値, CultureInfo.InvariantCulture))} |");
        }

        /// <summary>
        /// Markdown の表のセルを壊さない文字列にする
        /// </summary>
        /// <param name="p_値"></param>
        /// <returns></returns>
        private static string Get_表の値(string? p_値)
        {
            return string.IsNullOrEmpty(p_値) ? C_未指定表示 : p_値.Replace(C_表区切り, C_表区切りエスケープ).Replace(C_書式, C_空白).Replace(C_改行, C_空白);
        }

        /// <summary>
        /// アセンブリ統計を JSON の本文へ追加する
        /// </summary>
        /// <param name="p_文"></param>
        /// <param name="p_名前"></param>
        /// <param name="p_統計"></param>
        /// <param name="p_最小長"></param>
        private static void V_追加_統計(StringBuilder p_文, string p_名前, アセンブリ統計 p_統計, int p_最小長)
        {
            V_追加(p_文, C_項目_0, p_名前);
            V_追加(p_文, C_JSON書式_minimum_sequence_length, p_最小長);
            V_追加(p_文, C_JSON書式_sequences_区切りあり, p_統計.A_配列数);
            V_追加(p_文, C_JSON書式_total_length_区切りあり, p_統計.A_総延長);
            V_追加(p_文, C_JSON書式_largest_区切りあり, p_統計.A_最大長);
            V_追加(p_文, C_JSON書式_smallest, p_統計.A_最小長);
            V_追加(p_文, C_JSON書式_n50_区切りあり, p_統計.A_N50);
            V_追加(p_文, C_JSON書式_l50_区切りあり, p_統計.A_L50);
            V_追加(p_文, C_JSON書式_gc_percent_区切りなし, Get_数値(p_統計.A_GC率));
            _ = p_文.AppendLine(C_JSON内側オブジェクト終了と区切り);
        }

        /// <summary>
        /// TSV に出す固定の種別名
        /// </summary>
        /// <param name="p_種別"></param>
        /// <returns></returns>
        private static string Get_種別コード(曖昧箇所の種別 p_種別)
        {
            return p_種別 switch
            {
                曖昧箇所の種別.支持なし => C_項目_no_support,
                曖昧箇所の種別.僅差 => C_項目_insufficient_margin,
                曖昧箇所の種別.経路が一意でない => C_項目_multiple_paths,
                曖昧箇所の種別.到達不能 => C_項目_unreachable,
                曖昧箇所の種別.アンカー不足 => C_項目_anchor_missing,
                曖昧箇所の種別.リード無し => C_項目_no_local_reads,
                曖昧箇所の種別.探索打切り => C_項目_search_cutoff,
                _ => C_項目_inside_repeat,
            };
        }

        /// <summary>
        /// 工程ごとの資源使用量をレポートへ足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_フェーズ計測">工程ごとの資源使用量</param>
        private static void V_追加_フェーズ計測(StringBuilder p_文, IReadOnlyList<フェーズ計測> p_フェーズ計測)
        {
            _ = p_文.AppendLine(C_JSON書式_phase_timings);
            for (var i = 0; i < p_フェーズ計測.Count; i++)
            {
                var l_計測 = p_フェーズ計測[i];
                var l_末尾 = i == p_フェーズ計測.Count - 1 ? string.Empty : C_区切り;
                V_追加(p_文, C_JSON行書式_stage, Get_文字列(l_計測.A_工程), Get_数値(l_計測.A_経過秒), Get_数値(l_計測.A_CPU秒), Get_数値(l_計測.A_確保MB), Get_数値(l_計測.A_ワーキングセットMB), Get_数値(l_計測.A_ピークワーキングセットMB), l_計測.A_世代2回収回数, l_末尾);
            }

            _ = p_文.AppendLine(C_項目_V_追加_フェーズ計測);
        }

        /// <summary>
        /// 自己検査の結果をレポートへ足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_整合性">自己検査の結果</param>
        private static void V_追加_自己検査(StringBuilder p_文, 整合性検査結果? p_整合性)
        {
            if (p_整合性 is not { } l_整合性)
            {
                _ = p_文.AppendLine(C_JSON書式_self_check);
                return;
            }

            _ = p_文.AppendLine(C_JSON書式_self_check_区切りなし);
            V_追加(p_文, C_JSON書式_trusted_kmers, l_整合性.A_信頼kmer数);
            V_追加(p_文, C_JSON書式_missing_kmers, l_整合性.A_取りこぼし数);
            V_追加(p_文, C_JSON書式_missing_percent, Get_数値(l_整合性.A_取りこぼし率));
            V_追加(p_文, C_JSON書式_excess_percent, Get_数値(l_整合性.A_出しすぎ率));
            _ = p_文.AppendLine(C_JSON内側オブジェクト終了と区切り);
        }

        /// <summary>
        /// 固定アンカーでの独立評価をレポートへ足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_評価">固定アンカーでの評価</param>
        private static void V_追加_アンカー評価(StringBuilder p_文, アセンブリ評価? p_評価)
        {
            if (p_評価 is not { } l_評価)
            {
                _ = p_文.AppendLine(C_JSON書式_anchor_evaluation);
                return;
            }

            _ = p_文.AppendLine(C_JSON書式_anchor_evaluation_区切りなし);
            V_追加(p_文, C_JSON書式_expected_copies, l_評価.A_期待延べ数);
            V_追加(p_文, C_JSON書式_missing_copies, l_評価.A_欠損延べ数);
            V_追加(p_文, C_JSON書式_excess_copies, l_評価.A_過剰延べ数);
            V_追加(p_文, C_JSON書式_completeness, Get_数値(l_評価.A_完全性));
            V_追加(p_文, C_JSON書式_accuracy, Get_数値(l_評価.A_正確性));
            V_追加(p_文, C_JSON書式_ng50, l_評価.A_NG50);
            V_追加(p_文, C_JSON書式_circular_replicons_区切りあり, l_評価.A_環状本数);
            V_追加(p_文, C_JSON書式_circularized_fraction, Get_数値(l_評価.A_環状化率));
            _ = p_文.AppendLine(C_JSON内側オブジェクト終了と区切り);
        }

        /// <summary>
        /// ポリッシュの結果をレポートへ足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_ポリッシュ">ポリッシュの結果</param>
        private static void V_追加_ポリッシュ(StringBuilder p_文, ポリッシュ統計? p_ポリッシュ)
        {
            if (p_ポリッシュ is not { } l_ポリッシュ)
            {
                _ = p_文.AppendLine(C_JSON書式_polish);
                return;
            }

            _ = p_文.AppendLine(C_JSON書式_polish_区切りなし);
            V_追加(p_文, C_JSON書式_mapped_reads, l_ポリッシュ.A_マップされたリード数);
            V_追加(p_文, C_JSON書式_rejected_reads, l_ポリッシュ.A_棄却されたリード数);
            V_追加(p_文, C_JSON書式_corrected_bases, l_ポリッシュ.A_訂正した塩基数);
            V_追加(p_文, C_JSON書式_median_depth, Get_数値(l_ポリッシュ.A_深度の中央値));
            V_追加(p_文, C_JSON書式_low_depth_fraction, Get_数値(l_ポリッシュ.A_深度不足率));
            _ = p_文.AppendLine(C_JSON内側オブジェクト終了と区切り);
        }

        /// <summary>
        /// 環状閉鎖の検証結果をレポートへ足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_閉鎖検証">環状閉鎖の検証結果</param>
        private static void V_追加_閉鎖検証(StringBuilder p_文, IReadOnlyList<環状閉鎖検証結果>? p_閉鎖検証)
        {
            if (p_閉鎖検証 is null)
            {
                _ = p_文.AppendLine(C_JSON書式_circular_closure);
                return;
            }

            _ = p_文.AppendLine(C_JSON書式_circular_closure_区切りなし);
            for (var i = 0; i < p_閉鎖検証.Count; i++)
            {
                var l_検証 = p_閉鎖検証[i];
                var l_末尾 = i == p_閉鎖検証.Count - 1 ? string.Empty : C_区切り;
                V_追加(p_文, C_JSON行書式_id, Get_文字列(l_検証.A_配列ID), l_検証.A_長さ, l_検証.A_跨いだリード数, l_検証.A_必要本数, l_検証.A_Has支持 ? C_JSON真 : C_JSON偽, l_末尾);
            }

            _ = p_文.AppendLine(C_JSON内側配列終了と区切り);
        }

        /// <summary>
        /// レポートへ 1 行足す
        /// </summary>
        /// <param name="p_文">組み立て中のレポート</param>
        /// <param name="p_書式">行の書式</param>
        /// <param name="p_引数">書式へ埋める値</param>
        private static void V_追加(StringBuilder p_文, string p_書式, params object?[] p_引数)
        {
            _ = p_文.AppendLine(string.Format(CultureInfo.InvariantCulture, p_書式, p_引数));
        }

        /// <summary>
        /// JSON へ書ける形の数値にして返す
        /// </summary>
        /// <param name="p_値">元の値</param>
        /// <returns>JSON へ書ける文字列</returns>
        private static string Get_数値(double p_値)
        {
            return p_値.ToString(C_ファイル名_0_Get_数値, CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// JSON の文字列リテラル
        /// </summary>
        /// <param name="p_値"></param>
        /// <returns></returns>
        private static string Get_文字列(string? p_値)
        {
            if (p_値 is null)
            {
                return C_JSON空値;
            }

            var l_文 = new StringBuilder(C_引用符);
            foreach (var l_文字 in p_値)
            {
                _ = l_文字 switch
                {
                    '"' => l_文.Append(C_JSON引用符エスケープ),
                    '\\' => l_文.Append(C_JSON逆斜線エスケープ),
                    '\n' => l_文.Append(C_JSON改行エスケープ),
                    '\r' => l_文.Append(C_JSON復帰文字エスケープ),
                    '\t' => l_文.Append(C_JSONタブエスケープ),
                    _ => l_文字 < ' ' ? l_文.Append(CultureInfo.InvariantCulture, $"\\u{(int)l_文字:x4}") : l_文.Append(l_文字),
                };
            }

            return l_文.Append('"').ToString();
        }

        #endregion
    }
}
