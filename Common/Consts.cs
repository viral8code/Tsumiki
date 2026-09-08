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

            public const string 前処理 = "-pp";

            public const string メモリ予算 = "-mem";

            public const string マルチk = "-mk";

            public const string マージ = "-mg";

            public const string 引き継ぎなし = "-nc";

            public const string SuperRead = "-sr";

            public const string 反復r_mer検証 = "-rv";

            public const string 局所アセンブリ = "-la";

            public const string 積極性モード = "-mode";

            public const string GFA出力 = "-gfa";
        }

        /// <summary>
        /// -mode が束ねる値。利用者が本当に決めたいのは「完全性と正確性の
        /// どちらに倒すか」の1軸であり、-pu(優勢閾値)と-pc(支持数閾値)を
        /// 個別の数値として露出するより、Unicycler の
        /// --mode {conservative, normal, bold} のようなプリセットのほうが
        /// 意図を素直に表せる。normal は既定値そのもの。
        /// </summary>
        public static class 積極性モード名
        {
            public const string 保守的 = "conservative";

            public const string 標準 = "normal";

            public const string 積極的 = "bold";
        }

        public const decimal 保守的モードのペア結合閾値 = 0.9m;

        public const ulong 保守的モードのペア支持数閾値 = 15;

        public const decimal 積極的モードのペア結合閾値 = 0.65m;

        public const ulong 積極的モードのペア支持数閾値 = 5;

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
        public const int マルチkで試す個数 = 6;

        /// <summary>
        /// 行き止まりの短い unitig を tip(エラー由来)とみなすカバレッジの上限。
        /// グラフ全体のカバレッジ基準値に対する比。
        /// </summary>
        public const double tipとみなすカバレッジ比 = 0.5;

        /// <summary>
        /// 自動で試す k の上限(リード長に対する比)。
        /// 単一 k の自動選択より大きく取る。カバレッジが十分あれば
        /// リード長に近い k のほうが良いことがあり、その領域は
        /// 実際に試さないと分からない。
        /// </summary>
        public const double マルチk上限のリード長比 = 0.9;

        /// <summary>自動で試す k の下限。</summary>
        public const int マルチkの下限 = 21;

        /// <summary>
        /// 試す価値があるとみなす k-mer カバレッジの下限。
        /// k を上げると1リードから取れる k-mer が減るため、カバレッジの薄い
        /// データで高い k を試しても時間を捨てるだけになる。最初の k の結果から
        /// 各 k のカバレッジを予測し、これを下回るものは実行前に捨てる。
        /// </summary>
        public const double マルチkの最小kmerカバレッジ = 10.0;

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

        /// <summary>フラグメント長の経験分布を刻むビン幅。</summary>
        public const int フラグメント長のビン幅 = 5;

        /// <summary>
        /// 同一のギャップ長から出たとみなす既知長のばらつき幅の下限。
        /// ライブラリが極端に狭いときに窓が潰れないようにする。
        /// </summary>
        public const int 既知長のばらつき幅の下限 = 25;

        /// <summary>
        /// スキャフォールド辺を認めるのに必要な、距離が揃っているペアの本数。
        /// 反復解決が使う -pc とは数える対象が違うので別に持つ。
        /// </summary>
        public const ulong スキャフォールド支持数の下限 = 3;

        /// <summary>
        /// 短い反復解決の拒否権(-rv)で使う r-mer 長を、そのkでのアセンブリの
        /// k 長にこれだけ足して決める(r = k + この値)。
        ///
        /// head と repeat、repeat と tail は de Bruijn グラフの辺である以上、
        /// 必ず k-1 塩基を共有しており、その共有区間は repeat 自身の配列にも
        /// そのまま現れる。r <= k だと、接合点を跨ぐと判定した窓もこの共有区間の
        /// 内側に収まってしまい、head だけ・repeat だけを読んだリードでも
        /// 真になる(=対応付けの正しさを何も検定できない)。r を k より
        /// 確実に長く取ることで、共有区間の外側まで踏み込んだリードでなければ
        /// 真になり得ない窓だけを見られるようにする。
        /// ulong に 2bit パックする都合上 r は 32 が上限で、k + この値が
        /// それを超えるkでは検証自体をスキップする(高い k では反復自体が
        /// 少なく、他の判定で十分間に合っていることが多い)。
        /// </summary>
        public const int rMer長のk超過分の既定値 = 10;

        /// <summary>
        /// 短い反復解決の拒否権で、経路の接合点が「実際にリードに読まれている」と
        /// 認めるのに必要な、接合点を跨ぐ r-mer の最小本数。
        /// </summary>
        public const int r_mer接合点支持の閾値の既定値 = 4;

        public static readonly string ヘルプテキスト = $"""
            {概要テキスト}
            Usage: tsumiki {引数キー.リード1のパス} <path> [{引数キー.リード2のパス} <path>] [options]

            # Input
            {$"{引数キー.リード1のパス} <path>",-20}forward fastq(.gz) path (required; also use this for single-end reads)
            {$"{引数キー.リード2のパス} <path>",-20}reverse fastq(.gz) path
            {$"{引数キー.曖昧塩基を許容}",-20}allow ambiguous bases (e.g. N) in reads (default: false)

            # k-mer / quality
            {$"{引数キー.k長} <int[,int...]>",-20}k-mer length; a comma-separated list (e.g. 31,63,95) tries each and keeps the best, like {引数キー.マルチk} (default: auto from read length, capped at {自動k長の上限})
            {$"{引数キー.kmerカットオフ} <int>",-20}minimum k-mer count to trust (default: auto from the k-mer spectrum)
            {$"{引数キー.Phredオフセット} <int>",-20}phred score base ({string.Join(" or ", 許容Phredオフセット)}) (default: {Phredオフセットの既定値})
            {$"{引数キー.クオリティカットオフ} <int>",-20}minimum base quality to trust (default: {クオリティカットオフの既定値})
            {$"{引数キー.メモリ予算} <size>",-20}memory budget for k-mer counting, e.g. 2G, 512M (default: {Util.Get_表示用メモリサイズ(メモリ予算の既定値)})

            # Paired-end / scaffolding
            {$"{引数キー.インサートサイズ} <int>",-20}expected insert size (default: auto-estimated from mapped pairs)
            {$"{引数キー.ペア結合閾値} <decimal>",-20}minimum dominance ratio to accept a scaffold edge (default: {ペア結合閾値の既定値})
            {$"{引数キー.ペア支持数閾値} <int>",-20}minimum read-pair support to resolve a short repeat (default: {ペア支持数閾値の既定値})
            {引数キー.積極性モード} <{積極性モード名.保守的}|{積極性モード名.標準}|{積極性モード名.積極的}>
            {"",20}preset for {引数キー.ペア結合閾値}/{引数キー.ペア支持数閾値}: trades completeness for safety (default: {積極性モード名.標準})

            # Preprocessing
            {$"{引数キー.前処理}",-20}trim adapter read-through and cross-correct low-quality bases via mate overlap (default: false)
            {$"{引数キー.エラー訂正}",-20}correct reads from the k-mer spectrum before assembly (default: false)

            # Multi-k assembly
            {$"{引数キー.マルチk}",-20}assemble at several k values and keep the best (default: false; up to {マルチkで試す個数 + 1}x runtime)
            {$"{引数キー.引き継ぎなし}",-20}don't carry sequence from one k to the next (default: carry)
            {$"{引数キー.SuperRead}",-20}with paired-end reads, bridge pairs into synthetic long reads and carry those too (default: false)
            {$"{引数キー.マージ}",-20}splice sequence from other k values into junctions the best k left open (default: false; can raise misassemblies)

            # Repeat resolution safeguards
            {$"{引数キー.反復r_mer検証}",-20}require raw-read r-mer support before duplicating a short repeat (default: false)
            {$"{引数キー.局所アセンブリ}",-20}re-assemble unclosed scaffold gaps from reads mapped near their edges; safer than {引数キー.マージ} (default: false)

            # Output / misc
            {$"{引数キー.GFA出力}",-20}also write {GFAファイル名} (the unitig graph, GFA1) for viewers like Bandage (default: false)
            {$"{引数キー.一時ディレクトリ} <path>",-20}temp directory (default: {一時ディレクトリの既定値})
            {$"{引数キー.スレッド数} <int>",-20}worker threads for read loading (default: number of logical processors)
            {$"{引数キー.ヘルプ}",-20}show this text

            """;

        public const string ユニティグファイル名 = "unitigs.fasta";

        public const string コンティグファイル名 = "contigs.fasta";

        public const string スキャフォールドファイル名 = "scaffolds.fasta";

        public const string GFAファイル名 = "assembly.gfa";

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

        /// <summary>
        /// ギャップ充填(GapFiller・LocalAssembler共通)で、推定ギャップ長に
        /// 対して許容する誤差(塩基)。インサートサイズ推定のばらつきが
        /// そのままギャップ長推定のばらつきになるため、ぴったりの長さだけを
        /// 探すと現実にはまず当たらない。
        /// </summary>
        public const int ギャップ充填の長さの余裕幅 = 30;

        /// <summary>
        /// ギャップ充填(GapFiller・LocalAssembler共通)の対象とするギャップ長の
        /// 上限。これより長いギャップは探索空間が広すぎるうえ、推定長の誤差も
        /// 大きく一意に定まる見込みが薄いため対象外とする。
        /// </summary>
        public const int ギャップ充填のギャップ長上限 = 500;
    }
}
