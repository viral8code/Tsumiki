namespace Tsumiki.Common
{
    internal class Consts
    {
        public const string バージョン = "1.0";

        public static readonly List<string> 作者一覧 = [
            "viral",
            ];

        public static readonly string 概要テキスト = $"""

            Tsumiki is a genome assembler.
            author: {string.Join(", ", 作者一覧)}
            version: {バージョン}

            """;

        /// <summary>
        /// コマンドライン引数のキー。定数をまとめるための入れ子であり、
        /// 値を持つ型(エンティティ)ではないため Model 配下へは展開しない。
        /// </summary>
        public static class 引数キー
        {
            public const string リード1のパス = "-1";

            public const string リード2のパス = "-2";

            public const string k長 = "-k";

            public const string kmerカットオフ = "-kc";

            public const string Phredオフセット = "-p";

            public const string クオリティカットオフ = "-q";

            public const string インサートサイズ = "-i";

            public const string 曖昧塩基を許容 = "-ab";

            public const string ヘルプ = "-h";

            public const string 一時ディレクトリ = "-t";

            public const string スレッド数 = "-th";

            public const string ペア結合閾値 = "-pu";

            public const string ペア支持数閾値 = "-pc";

            public const string エラー訂正 = "-ec";

            public const string メモリ予算 = "-mem";

            public const string マルチk = "-mk";
        }

        public const string インサートサイズ未指定表示 = "unspecified";

        /// <summary>
        /// -k も、リード長の標本抽出も当てにできなかった場合の最後の拠り所。
        /// 通常は実際のリード長から自動選択されるため、この値は使われない。
        /// </summary>
        public const int k長の既定値 = 31;

        /// <summary>
        /// 自動選択する k 長の、リード長に対する比。大きいほど反復を跨げるが、
        /// 1リードから取れる k-mer の本数が減ってカバレッジが痩せる。
        /// </summary>
        public const double 自動k長のリード長比 = 0.6;

        /// <summary>
        /// 自動選択する k 長の上限。k が 64 を超えると 2bit パックが
        /// UInt128 に収まらず高速経路から外れるため、そこで頭を抑える
        /// (偶数を避けるので実際に選ばれる最大値は 63)。
        /// </summary>
        public const int 自動k長の上限 = 63;

        /// <summary>
        /// k 長の自動選択を行うために最低限必要なリード長。
        /// これより短いリードでは k を十分に取れず、自動選択しても意味がない。
        /// </summary>
        public const int 自動k長に必要な最小リード長 = 32;

        /// <summary>
        /// -mk で試す k の個数。増やすほど実行時間が線形に伸びる。
        /// ヘルプに実行時間の目安として出るため、ここに置いている。
        /// </summary>
        public const int マルチkで試す個数 = 3;

        /// <summary>
        /// -kc も、k-mer スペクトルの解析も当てにできなかった場合の最後の拠り所。
        /// 通常は k-mer スペクトルから自動選択される。
        /// </summary>
        public const int kmerカットオフの既定値 = 2;

        /// <summary>
        /// k-mer カウント時にメモリ上へ保持するカウントの総量(バイト)。
        /// 大きいほどフラッシュ回数が減って I/O が軽くなる代わりメモリを使う。
        /// </summary>
        public const long メモリ予算の既定値 = 768L * 1024 * 1024;

        public const int Phredオフセットの既定値 = 33;

        public const int クオリティカットオフの既定値 = 1;

        public const string 一時ディレクトリの既定値 = "temp";

        public static readonly int[] 許容Phredオフセット = [33, 64];

        public const decimal ペア結合閾値の既定値 = 0.8m;

        public const ulong ペア支持数閾値の既定値 = 10;

        public static readonly string ヘルプテキスト = $"""
            {概要テキスト}

            # Arguments
            {引数キー.リード1のパス} [path] : forward fastq(.gz) path (required) (when using single reads, set the path using this argument)
            {引数キー.リード2のパス} [path] : backward fastq(.gz) path
            {引数キー.k長} [integer] : length of k-mer (default : auto-selected from the observed read length, capped at {自動k長の上限}; falls back to {k長の既定値})
            {引数キー.kmerカットオフ} [integer] : threshold of k-mer count (use kmers with this value or higher) (default : auto-selected from the k-mer count spectrum; falls back to {kmerカットオフの既定値})
            {引数キー.Phredオフセット} [integer] : base of phred score ({string.Join(" or ", 許容Phredオフセット)}) (default : {Phredオフセットの既定値})
            {引数キー.クオリティカットオフ} [integer] : threshold of base quality (use kmers with this value or higher) (default : {クオリティカットオフの既定値})
            {引数キー.メモリ予算} [decimal] : memory budget for k-mer counting (e.g. 2G, 512M; a bare number means MB). Raise it to reduce disk I/O, lower it to fit a smaller machine (default : {Util.Get_表示用メモリサイズ(メモリ予算の既定値)})
            {引数キー.インサートサイズ} : excepted insert size of pair-end reads (default : {インサートサイズ未指定表示}, auto-estimated from mapped pairs when possible)
            {引数キー.一時ディレクトリ} [path] : temp directory (default : {一時ディレクトリの既定値})
            {引数キー.スレッド数} [integer] : number of worker threads used for loading reads (default : number of logical processors)
            {引数キー.ペア結合閾値} [decimal] : minimum ratio of the best-supported pair-end scaffold edge among all candidates for a node (default : {ペア結合閾値の既定値})
            {引数キー.ペア支持数閾値} [integer] : minimum read-pair support required for a pair-end scaffold edge (default : {ペア支持数閾値の既定値})
            {引数キー.マルチk} : assemble at several k and keep the best one, judged without a reference. The best k depends on how repetitive the genome is, which cannot be known from the reads alone, so the only way to find it is to try (costs roughly {マルチkで試す個数 + 1}x the runtime) (default : false)
            {引数キー.エラー訂正} : run k-mer-spectrum-based read error correction before assembly (default : false)
            {引数キー.ヘルプ} : output this text (default : false)

            """;

        public const string ユニティグファイル名 = "unitigs.fasta";

        public const string コンティグファイル名 = "contigs.fasta";

        public const string スキャフォールドファイル名 = "scaffolds.fasta";

        /// <summary>
        /// 塩基の内部表現。A/C/G/T は塩基記号そのものなので英字のまま残す
        /// (日本語にするとかえって読みにくいため)。
        /// </summary>
        public static class 塩基ID
        {
            public const byte A = 1;
            public const byte C = 2;
            public const byte G = 3;
            public const byte T = 4;
        }

        public const ulong 進捗ログ間隔 = 100_000;

        public const byte 無効な塩基 = 5;

        public const int ユニティグ数の上限 = 100_000;

        /// <summary>
        /// インサートサイズが未指定の場合の自動推定に使う、単一unitigへ両リードが
        /// マップされたペアの最小標本数。これに満たない場合は推定を諦め、
        /// ペアエンド由来のスキャフォールディングをスキップする。
        /// </summary>
        public const int インサートサイズ標本数の下限 = 30;

        /// <summary>
        /// ギャップ長が推定上0以下になった場合に最低限挿入するNの数。
        /// </summary>
        public const int ギャップ長の下限 = 1;
    }
}
